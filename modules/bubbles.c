/* Gerris - The GNU Flow Solver
 * Copyright (C) 2010-2012 Daniel Fuster/CNRS
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.  
 */

#include "particulatecommon.h"
#include "ftt.h"

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

/* GfsBubbleParams: Header */
typedef struct _GfsBubbleParams GfsBubbleParams;

typedef enum {
  GFS_RP,
  GFS_RPKM,
  GFS_NORP
} GfsRPType;

typedef void (* GfsRayleighPlessetFunc)       (double t, const double y[], double f[], void * params);

struct _GfsBubbleParams {
  GfsEvent parent;

  gdouble RKeps, cl,dtfactor, sigma, visc ;
  GfsVariable * p;
  const gsl_odeiv_step_type * odeiv_type;
  GfsRPType RPType;
  GfsRayleighPlessetFunc RKeq;
};

#define GFS_BUBBLE_PARAMS(obj)            GTS_OBJECT_CAST (obj,\
					         GfsBubbleParams,\
					         gfs_bubble_params_class ())
#define GFS_IS_BUBBLE_PARAMS(obj)         (gts_object_is_from_class (obj,\
						 gfs_bubble_params_class ()))

GfsEventClass * gfs_bubble_params_class  (void);


/* GfsBubble: Header */

typedef struct _GfsBubble GfsBubble;

struct _GfsBubble {
  /*< private >*/
  GfsParticulate parent;
  gdouble rliq;
  GfsBubbleParams * params;
  
  /*< public >*/
  gdouble velR, p0, R0, vol_liq, dpotdt;
};

#define GFS_BUBBLE(obj)            GTS_OBJECT_CAST (obj, GfsBubble, gfs_bubble_class ())
#define GFS_IS_BUBBLE(obj)         (gts_object_is_from_class (obj, gfs_bubble_class ()))

static GfsEventClass * gfs_bubble_class  (void);


/* GfsBubble: Object */
/* The radius of each bubble varies according to the Rayleigh-Plesset equation */

typedef struct {
  gdouble liqpres, liqdens;
  GfsBubble * bubble;
  GfsDomain * domain;
} RPData;

gdouble static p_state_ec (gdouble p0, gdouble R0, gdouble rb) {
    /* check negative radius: what to do in this case? */
    if (rb <= 1.e-3*R0 ) 
      rb = 1.e-2*R0;

    return p0*pow (R0/rb, 3.*1.4);
}

static void RPeq (double t, const double y[], double f[], void * params)
{
  RPData * rp = (RPData *) params;
  gdouble pbubble = p_state_ec (rp->bubble->p0, rp->bubble->R0, y[0]);

  f[1] = ((pbubble - 2.*rp->bubble->params->sigma/y[0] + 4.*rp->bubble->params->visc*y[1]/y[0] - rp->liqpres)/rp->liqdens - 3./2.*y[1]*y[1])/y[0];
}

static void RPKMeq (double t, const double y[], double f[], void * params)
{
  RPData * rp = (RPData *) params;
  gdouble pbubble = p_state_ec (rp->bubble->p0, rp->bubble->R0, y[0]);
  f[1] = (pbubble - 2.*rp->bubble->params->sigma/y[0] + 4.*rp->bubble->params->visc*y[1]/y[0] - rp->liqpres)/rp->liqdens;
  f[1] *= ( 1. + y[1]/rp->bubble->params->cl);
  f[1] -= 3./2.*y[1]*y[1]*( 1. - y[1]/rp->bubble->params->cl/3.);
  f[1] /= (y[0]*( 1. - y[1]/rp->bubble->params->cl));
}

static void NORPeq (double t, const double y[], double f[], void * params)
{
   f[1] = 0.;
}

static int func (double t, const double y[], double f[], void * params)
{
  RPData * rp = (RPData *) params;

  GfsParticulate * particulate = GFS_PARTICULATE (rp->bubble);
  GfsParticle * p = GFS_PARTICLE (particulate);

  p->pos.x = y[2];
  p->pos.y = y[3];
  p->pos.z = y[4];
  particulate->vel.x = y[5];
  particulate->vel.y = y[6];
  particulate->vel.z = y[7];
  FttComponent c;
  for (c = 0; c < 3; c++)
    (&particulate->force.x)[c] = 0.;

  f[0] = y[1];

  /* interface acceleration- incompressible RP equation */
  (* rp->bubble->params->RKeq) (t, y, f, params);

  f[2] = y[5];
  f[3] = y[6];
  f[4] = y[7];
  gts_container_foreach (GTS_CONTAINER (particulate->forces),
      (GtsFunc) compute_forces, particulate);
  f[5] = particulate->force.x/particulate->mass;
  f[6] = particulate->force.y/particulate->mass;
#if !FTT_2D
  f[7] = particulate->force.z/particulate->mass;
#else
  f[7] = 0.;
#endif

  return GSL_SUCCESS;
}

/* jacobian matrix */
int static jac (double t, const double y[], double *dfdy, 
    double dfdt[], void *params)
{
  RPData * rp = (RPData *) params;
  gdouble pbubble = p_state_ec (rp->bubble->p0, rp->bubble->R0, y[0]);
  gdouble term10  = 2.*rp->bubble->params->sigma + rp->liqpres*y[0] - 8.*rp->bubble->params->visc*y[1] - 0.5*rp->liqdens*y[0]*pow(y[1],2);
  term10 /= y[0];
  term10 -= 4.*pbubble;
  term10 *= y[1];
  term10 += rp->bubble->params->cl*(2.*rp->bubble->params->sigma + rp->liqpres*y[0] - 8.*rp->bubble->params->visc*y[1] - 1.5*rp->liqdens*y[0]*pow(y[1],2))/y[0];
  term10 -= 4.*rp->bubble->params->cl*pbubble;
  term10 /= (rp->liqdens*pow(y[0],2)*(rp->bubble->params->cl-y[1]));

  gdouble term11 = pow(rp->bubble->params->cl,2)*(4.*rp->bubble->params->visc-3.*rp->liqdens*y[0]*y[1])/y[0];
  term11 += pow(y[1],2)*(-4.*rp->bubble->params->visc-rp->liqdens*y[0]*y[1])/y[0];
  term11 += rp->bubble->params->cl*2.*pbubble;
  term11 += rp->bubble->params->cl*(-4.*rp->bubble->params->sigma-2.*rp->liqpres*y[0] + 8.*rp->bubble->params->visc*y[1] + 3.*rp->liqdens*y[0]*pow(y[1],2))/y[0];
  term11 /= (rp->liqdens*y[0]*pow(rp->bubble->params->cl-y[1],2));
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
  gsl_matrix * m = &dfdy_mat.matrix; 
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 1, 0, term10);
  gsl_matrix_set (m, 1, 1, term11);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

static gboolean gfs_bubble_event (GfsEvent * event, 
    GfsSimulation * sim)
{

  GfsParticle * p = GFS_PARTICLE (event);
  GfsParticulate * particulate = GFS_PARTICULATE (event);

  GfsBubble * bubble = GFS_BUBBLE (event);
  GfsDomain * domain = GFS_DOMAIN (sim);
  if (bubble->params == NULL) {
    GSList * i = sim->events->items;
    gboolean bubbleparams = FALSE;
    while (i) {
      if (GFS_IS_BUBBLE_PARAMS (i->data)) {
        bubbleparams = TRUE;
        bubble->params = GFS_BUBBLE_PARAMS (i->data);
      }
      i = i->next;
    }
    if (!bubbleparams) {
      bubble->params = GFS_BUBBLE_PARAMS (gts_object_new (GTS_OBJECT_CLASS (gfs_bubble_params_class ())));
      bubble->params->p = gfs_variable_from_name (domain->variables, "P");
    }

  }

  GfsVariable * liqpres = bubble->params->p;

  FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
  if (cell == NULL) 
    return TRUE;
  gdouble liq_rho = sim->physical_params.alpha ? 1./
    gfs_function_value (sim->physical_params.alpha, cell) : 1.;

  FttVector pos = p->pos;
  gfs_simulation_map (sim, &pos);

  if (cell) 
  	bubble->dpotdt = gfs_interpolate (cell, p->pos, liqpres);
  else
	bubble->dpotdt = bubble->p0;

  gdouble point_pres = bubble->dpotdt;

  RPData rp = { point_pres, liq_rho, bubble, domain };

  gsl_odeiv_step * s    = gsl_odeiv_step_alloc (bubble->params->odeiv_type, 8);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new (bubble->params->RKeps, bubble->params->RKeps);
  gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (8);

  gsl_odeiv_system sys = {func, NULL, 8, &rp};

  gdouble t = sim->time.t;
  gdouble t1 = t + sim->advection_params.dt;
  gdouble h = sim->advection_params.dt*0.1; /* better criterion?? */
  /* variables R, dot{R} */
  gdouble y[8] = { pow(3./(4.*M_PI)*particulate->volume,1./3.) , bubble->velR, p->pos.x, p->pos.y, p->pos.z, 
    particulate->vel.x, particulate->vel.y, particulate->vel.z };

  while (t < t1) {
    gdouble told = t;
    int status = gsl_odeiv_evolve_apply (e, c, s,
        &sys, &t, t1, &h, y);
    if (bubble->params->dtfactor > 0.)
      sim->time.dtmax  = (t - told)*bubble->params->dtfactor;
    gdouble pbubble = p_state_ec (bubble->p0, bubble->R0, y[0]);
    if (status != GSL_SUCCESS) 
      g_error ("Error in the RK method");
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);

  bubble->velR = y[1];
  if (y[0] <= 0.001*bubble->R0 ) 
      y[0] = 1.e-2*bubble->R0;
  particulate->volume = 4./3.*M_PI*y[0]*y[0]*y[0];

  p->pos.x = y[2];
  p->pos.y = y[3];
  p->pos.z = y[4];
  particulate->vel.x = y[5];
  particulate->vel.y = y[6];
  particulate->vel.z = y[7];

  /* trick to impose the pressure with/without interactions. Improve?? */
  //    bubble->dpotdt = gfs_interpolate (cell, p->pos, liqpres);
  return TRUE;
} 

static void gfs_bubble_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_bubble_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_bubble_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  GfsBubble * p = GFS_BUBBLE (*o);
  GfsParticulate * part = GFS_PARTICULATE (*o);
  gdouble L = gfs_object_simulation (*o)->physical_params.L;
    
  p->vol_liq = 0;
  p->dpotdt  = 101300.;
  p->R0 = pow (part->volume*3./(4.*M_PI), 1./3.);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (radial velocity)");
    return;
  }
  p->velR = atof (fp->token->str)/L;
  gts_file_next_token (fp);
  
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (reference pressure)");
    return;
  }
  p->p0 = atof (fp->token->str)*L;
  gts_file_next_token (fp);
}

static void gfs_bubble_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_bubble_class ())->parent_class->write) (o, fp); 
  GfsBubble * p = GFS_BUBBLE (o);
  gdouble L = gfs_object_simulation (o)->physical_params.L;
  fprintf (fp, " %g %g", p->velR*L, p->p0/L);
}

static void gfs_bubble_class_init (GfsEventClass * klass)
{
  klass->event = gfs_bubble_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_bubble_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_bubble_write;
}


static GfsEventClass * gfs_bubble_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bubble_info = {
      "GfsBubble",
      sizeof (GfsBubble),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_bubble_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particulate_class ()),
				  &gfs_bubble_info);
  }
  return klass;
}

/* GfsBubbleParams: Object */

static void gfs_bubble_params_init (GfsBubbleParams * bparams)
{
  GFS_EVENT (bparams)->istep = 1;

  bparams->RKeps = 1.e-6;
  bparams->cl = 1500.;
  bparams->sigma = 0.;
  bparams->visc  = 0.;
  bparams->dtfactor = 0.;
  bparams->RPType = GFS_RP;
  bparams->RKeq = (GfsRayleighPlessetFunc) RPeq;

  bparams->odeiv_type = gsl_odeiv_step_rk8pd;
}


static void bubble_params_write (GtsObject * o, FILE * fp)
{
  
  (* GTS_OBJECT_CLASS (gfs_bubble_params_class ())->parent_class->write) (o, fp);

  GfsBubbleParams * bparams = GFS_BUBBLE_PARAMS (o);
  fprintf (fp,
           "{\n"
	   "  eps      = %g\n"
	   "  cl       = %g\n"
	   "  sigma    = %g\n"
	   "  visc     = %g\n"
	   "  dtfactor = %g\n"
	   "  p = %s\n"
	   "  RK = %s\n",
	  bparams->RKeps,
 	  bparams->cl,
 	  bparams->sigma,
 	  bparams->visc,
          bparams->dtfactor,
          GFS_VARIABLE(bparams->p)->name,
	  bparams->odeiv_type == gsl_odeiv_step_rk2 ?
	  "rk2" : 
          bparams->odeiv_type == gsl_odeiv_step_rk4 ?
          "rk4" :
          bparams->odeiv_type == gsl_odeiv_step_rkf45 ? 
          "rkf45" :
          bparams->odeiv_type == gsl_odeiv_step_rkck ?
          "rkck" :
          bparams->odeiv_type == gsl_odeiv_step_rk8pd ?
          "rk8pd" :
          bparams->odeiv_type == gsl_odeiv_step_rk2imp ?
          "rk2imp" :
          bparams->odeiv_type == gsl_odeiv_step_rk4imp ?
          "rk4imp" :
          bparams->odeiv_type == gsl_odeiv_step_bsimp ?
          "bsimp" :
          bparams->odeiv_type == gsl_odeiv_step_gear1 ?
          "gear1" :
          bparams->odeiv_type == gsl_odeiv_step_gear2 ?
          "gear2" :
          "none"
 	  );
  switch (bparams->RPType) {
  case GFS_RP: fputs ("  RP   = RP\n", fp); break;
  case GFS_RPKM: fputs ("  RP   = RPKM\n", fp); break;
  case GFS_NORP: fputs ("  RP   = none\n", fp); break;
  }

  fputc ('}', fp);
  
}

static void bubble_params_read (GtsObject ** o, GtsFile * fp)
{
  GfsBubbleParams * bparams = GFS_BUBBLE_PARAMS (*o);

  (* GTS_OBJECT_CLASS (gfs_bubble_params_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  gfs_bubble_params_init (bparams);
  gchar * RK = NULL, * pname = NULL, * RPtype = NULL;
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "eps", TRUE, &bparams->RKeps},
    {GTS_DOUBLE, "cl",  TRUE, &bparams->cl},
    {GTS_DOUBLE, "sigma",  TRUE, &bparams->sigma},
    {GTS_DOUBLE, "visc",  TRUE, &bparams->visc},
    {GTS_DOUBLE, "dtfactor",  TRUE, &bparams->dtfactor},
    {GTS_STRING, "p",  TRUE, &pname},
    {GTS_STRING, "RK",  TRUE, &RK},
    {GTS_STRING, "RP",  TRUE, &RPtype},
    {GTS_NONE}
  };

  gts_file_assign_variables (fp, var);

  if (RK) {
    if (!strcmp (RK, "rk2"))
      bparams->odeiv_type = gsl_odeiv_step_rk2;
    else if (!strcmp (RK, "rk4"))
      bparams->odeiv_type = gsl_odeiv_step_rk4;
    else if (!strcmp (RK, "rkf45"))
      bparams->odeiv_type = gsl_odeiv_step_rkf45;
    else if (!strcmp (RK, "rkck"))
      bparams->odeiv_type = gsl_odeiv_step_rkck;
    else if (!strcmp (RK, "rk8pd"))
      bparams->odeiv_type = gsl_odeiv_step_rk8pd;
    else if (!strcmp (RK, "rk2imp"))
      bparams->odeiv_type = gsl_odeiv_step_rk2imp;
    else if (!strcmp (RK, "rk4imp"))
      bparams->odeiv_type = gsl_odeiv_step_rk4imp;
    else if (!strcmp (RK, "bsimp"))
      bparams->odeiv_type = gsl_odeiv_step_bsimp;
    else if (!strcmp (RK, "gear1"))
      bparams->odeiv_type = gsl_odeiv_step_gear1;
    else if (!strcmp (RK, "gear2"))
      bparams->odeiv_type = gsl_odeiv_step_gear2;
    else if (fp->type != GTS_ERROR)
      gts_file_variable_error (fp, var, "Runge Kutta",
			       "unknown RK parameter `%s'", RK);
    g_free (RK);
  }

  if (RPtype) {
    if (!strcmp (RPtype, "RP")) {
      bparams->RPType = GFS_RP;
      bparams->RKeq   = RPeq;
    }
    else if (!strcmp (RPtype, "RPKM")) {
      bparams->RPType = GFS_RPKM;
      bparams->RKeq   = RPKMeq;
    }
    else if (!strcmp (RPtype, "none")) {
      bparams->RPType = GFS_NORP;
      bparams->RKeq = NORPeq;
    }
    else if (fp->type != GTS_ERROR)
      gts_file_variable_error (fp, var, "equation",
			       "unknown Rayleigh-Plesset type equation `%s'", RPtype);
    g_free (RPtype);
  }

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (bparams));
  if (pname)
	bparams->p = gfs_variable_from_name (domain->variables, pname);
  else
	bparams->p = gfs_variable_from_name (domain->variables, "P");

}

static void gfs_bubble_params_class_init (GtsObjectClass * klass)
{
  klass->read =    bubble_params_read;
  klass->write =   bubble_params_write;
}

GfsEventClass * gfs_bubble_params_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bubble_params_info = {
      "GfsBubbleParams",
      sizeof (GfsBubbleParams),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_bubble_params_class_init,
      (GtsObjectInitFunc) gfs_bubble_params_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_bubble_params_info);
  }

  return klass;
}

/* GfsBubbleFraction: header */

typedef struct _GfsBubbleFraction                GfsBubbleFraction;

struct _GfsBubbleFraction {
  /*< private >*/
  GfsParticulateField parent;

  /*< public >*/
  gdouble rliq;
  GfsFunction * kernel_function;
};

#define GFS_BUBBLE_FRACTION(obj)            GTS_OBJECT_CAST (obj,	\
							     GfsBubbleFraction, \
							     gfs_bubble_fraction_class ())
#define GFS_IS_BUBBLE_FRACTION(obj)         (gts_object_is_from_class (obj, \
              					             gfs_bubble_fraction_class ()))

GfsVariableClass * gfs_bubble_fraction_class  (void);

typedef struct {
  gdouble correction;
  GfsBubble * bubble;
  GfsVariable * v;
  GfsBubbleFraction * bf;
} BubbleData;

typedef struct {
  FttVector * pos;
  gdouble distance; 
} CondData;

/** \beginobject{GfsBubbleFraction} */
/* do it more general? */

static void voidfraction_from_bubbles (FttCell * cell, BubbleData * p)
{
  FttVector pos;
  ftt_cell_pos (cell, &pos);
  distance_normalization (&pos, GFS_PARTICULATE (p->bubble));
  GFS_VALUE (cell, p->v) += GFS_PARTICULATE (p->bubble)->volume*
    gfs_function_spatial_value (p->bf->kernel_function, &pos)/p->correction;
}

static void kernel_volume (FttCell * cell, BubbleData * p)
{
  gdouble cellvol = gfs_cell_volume (cell, p->v->domain);

  p->bubble->vol_liq += cellvol;

  /* correction term to make a discretely conservative kernel */
  FttVector pos;
  ftt_cell_pos (cell, &pos);
  distance_normalization (&pos, GFS_PARTICULATE (p->bubble));
  p->correction += gfs_function_spatial_value (p->bf->kernel_function, &pos)*cellvol;
}

static gboolean cond_bubble (FttCell * cell, gpointer data)
{
  CondData * p = data;
  FttVector pos;
  ftt_cell_pos (cell, &pos);
  gdouble radeq;
  gdouble size = ftt_cell_size (cell)/2.;

#if FTT_2D
  radeq = size*sqrt(2.);
#else
  radeq = size*sqrt(3.);
#endif /* 3D */


  if (ftt_vector_distance (&pos, p->pos) - radeq <= p->distance) 
    return TRUE;

  /* Check also if the bubble is inside the cell*/
  if (p->pos->x > pos.x + size || p->pos->x < pos.x - size ||
      p->pos->y > pos.y + size || p->pos->y < pos.y - size 
#if !FTT_2D
      || p->pos->z > pos.z + size || p->pos->z < pos.z - size 
#endif
     )
    return FALSE;

  return TRUE;
}

static gboolean bubble_fraction_event (GfsEvent * event, 
    GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_particulate_field_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsVariable * v = GFS_VARIABLE (event);
    GfsParticulateField * pfield = GFS_PARTICULATE_FIELD (v);
    GfsBubbleFraction * bf = GFS_BUBBLE_FRACTION (event);

    /* Reset variable */
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
        (FttCellTraverseFunc) gfs_cell_reset, v);
    /* Loop over the list of particles in the selected object */
    GSList * i = GFS_EVENT_LIST (pfield->plist)->list->items;
    while (i) {
      GfsBubble * bubble = GFS_BUBBLE (i->data);
      bubble->vol_liq = 0;
      bubble->rliq = pow (GFS_PARTICULATE(i->data)->volume*3./(4.*M_PI), 1./3.)*bf->rliq;
      BubbleData p = { 0, bubble, v, bf };
      CondData cd = { &GFS_PARTICLE (i->data)->pos, bubble->rliq };
      gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
          (FttCellTraverseFunc) kernel_volume, &p,
          cond_bubble, &cd);          
      gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
          (FttCellTraverseFunc) pfield->voidfraction_func, &p,
          cond_bubble, &cd);          
      i = i->next;
    }
    return TRUE;
  }
  return FALSE;
}

static void bubble_fraction_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_BUBBLE_FRACTION (o)->kernel_function));

  (* GTS_OBJECT_CLASS (gfs_bubble_fraction_class ())->parent_class->destroy) (o); 
}

static void bubble_fraction_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_bubble_fraction_class ())->parent_class->read) (o, fp); 
  if (fp->type == GTS_ERROR)
    return;

  GfsBubbleFraction * b = GFS_BUBBLE_FRACTION(*o);

  if (fp->type != '{') {
    gts_file_error (fp, "expecting an opening brace");
    return;
  }
  fp->scope_max++;
  gts_file_next_token (fp);

  while (fp->type != GTS_ERROR && fp->type != '}') {
    if (fp->type == '\n') {
      gts_file_next_token (fp);
      continue;
    }
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a keyword");
      return;
    }  
    else if (!strcmp (fp->token->str, "rkernel")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting '='");
	return;
      }    
      gts_file_next_token (fp);
      b->rliq = atof (fp->token->str);
      gts_file_next_token (fp);
    }
    else if (!strcmp (fp->token->str, "kernel")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting '='");
	return;
      }
      gts_file_next_token (fp);
      gfs_function_read (b->kernel_function, gfs_object_simulation (*o), fp);
    }
    else {
      gts_file_error (fp, "unknown keyword `%s'", fp->token->str);
      return;
    }
  }
  if (fp->type == GTS_ERROR)
    return;
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return;
  }
  fp->scope_max--;
  gts_file_next_token (fp);
}

static void bubble_fraction_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_bubble_fraction_class ())->parent_class->write) (o, fp); 
  fprintf (fp, " { rkernel = %g ", GFS_BUBBLE_FRACTION (o)->rliq);
  fputs (" kernel =", fp);
  gfs_function_write (GFS_BUBBLE_FRACTION (o)->kernel_function, fp);
  fputc ('}', fp);
}

static void bubble_fraction_init (GfsVariable * v)
{
  v->units = 0.;
  GFS_PARTICULATE_FIELD (v)->voidfraction_func = voidfraction_from_bubbles;
  GFS_BUBBLE_FRACTION (v)->kernel_function = gfs_function_new (gfs_function_spatial_class (), 0.);
}

static void bubble_fraction_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = bubble_fraction_event;
  klass->destroy = bubble_fraction_destroy;
  klass->read =  bubble_fraction_read;
  klass->write = bubble_fraction_write;
}

GfsVariableClass * gfs_bubble_fraction_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bubble_fraction_info = {
      "GfsBubbleFraction",
      sizeof (GfsBubbleFraction),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) bubble_fraction_class_init,
      (GtsObjectInitFunc) bubble_fraction_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particulate_field_class ()),
        &gfs_bubble_fraction_info);
  }

  return klass;
}

/** \endobject{GfsBubbleFraction} */

/* GfsBubbleFractionDt: header */

#define GFS_IS_BUBBLE_FRACTION_DT(obj)         (gts_object_is_from_class (obj, \
						     gfs_bubble_fraction_dt_class ()))

GfsVariableClass * gfs_bubble_fraction_dt_class (void);

/** \beginobject{GfsBubbleFractionDt} */

static void dVpdt_from_particles (FttCell * cell, BubbleData * p)
{
  FttVector pos;
  ftt_cell_pos (cell, &pos);
  distance_normalization (&pos, GFS_PARTICULATE (p->bubble));
  gdouble rad = pow (3.0*GFS_PARTICULATE (p->bubble)->volume/(4.0*M_PI), 1./3.);
  GFS_VALUE (cell, p->v) += 4.*M_PI*pow(rad,2)*p->bubble->velR*gfs_function_spatial_value (p->bf->kernel_function, &pos)/p->correction;
}

static void bubble_fraction_dt_init (GtsObject * o)
{
  GFS_PARTICULATE_FIELD (o)->voidfraction_func = dVpdt_from_particles;
}

GfsVariableClass * gfs_bubble_fraction_dt_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bubble_fraction_dt_info = {
      "GfsBubbleFractionDt",
      sizeof (GfsBubbleFraction),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) bubble_fraction_dt_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS ( gfs_bubble_fraction_class ()),
				  &gfs_bubble_fraction_dt_info);
  }
  return klass;
}

/** \endobject{GfsBubbleFractionDt} */

/** \beginobject{GfsBubbleInteractions} */

/* GfsBubbleInteractions: header */

typedef struct _GfsBubbleInteractions   GfsBubbleInteractions;

struct _GfsBubbleInteractions {
  /*< private >*/
  GfsEvent parent;

  /*< public >*/
  GfsParticleList * plist;
};

#define GFS_BUBBLE_INTERACTIONS(obj)            GTS_OBJECT_CAST (obj, GfsBubbleInteractions, gfs_bubble_interactions_class ())
#define GFS_IS_BUBBLE_INTERACTIONS(obj)         (gts_object_is_from_class (obj, gfs_bubble_interactions_class ()))

GfsVariableClass * gfs_bubble_interactions_class  (void);

/* GfsBubbleInteractions: object */

typedef struct {
  GfsVariable * v, * p;
  gdouble Rb3,Rb,sumRi2v,sumRi2;
  GtsSListContainer * listp;
  GfsBubble * b;
} InteractionData;


typedef struct {
  FttVector * pos;
  gdouble distance,pcell;
  GfsVariable * p;
} PcellData;


static void averaged_pressure (FttCell * cell, gpointer data)
{
  PcellData * pcell = data;
  pcell->pcell += ftt_cell_volume (cell)*GFS_VALUE(cell,pcell->p);
  return;
}

static gdouble get_pcell (FttCell * cell, InteractionData * data) {

  /*ensure that it executes after GfsBubbleFraction to obtain bubble->rliq*/
  PcellData pcell = {  &GFS_PARTICLE (data->b)->pos,pow(data->b->vol_liq*3./(4.*M_PI), 1./3.),0, data->p };
  CondData cd = { &GFS_PARTICLE (data->b)->pos,data->b->rliq};
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (data->p));

  gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
      (FttCellTraverseFunc) averaged_pressure, &pcell,
      cond_bubble, &cd);          
  return pcell.pcell/data->b->vol_liq;
}
static void solve_one_bubble (FttCell * cell, InteractionData * data) {

  GSList * i = data->listp->items;
  GfsBubble * b = GFS_BUBBLE(i->data);
  data->b = b; //Why did I do that ?? check
  gdouble Ri= pow(3.*GFS_PARTICULATE(i->data)->volume/4./M_PI,1./3.);
  gdouble Rb3 = data->Rb3;
  gdouble Rb  = data->Rb;
  gdouble sumRi2v  = data->sumRi2v;
  gdouble sumRi2  = data->sumRi2;

  GfsSimulation * sim = gfs_object_simulation(data->p); 
  gdouble liq_rho = sim->physical_params.alpha ? 1./
    gfs_function_value (sim->physical_params.alpha, cell) : 1.;

  gdouble prescell = get_pcell(cell,data); 

  gdouble radeq = pow(3.*b->vol_liq/(4.*M_PI),1./3.);
  gdouble aux1 = 3./2.*Rb3/(pow(radeq,3)-Rb3)*pow(sumRi2v/sumRi2,2)*(1.-Rb/radeq);
  gdouble aux2 = 3./2.*Rb*(pow(radeq,2)-pow(Rb,2))/(pow(radeq,3)-Rb3);

  gdouble dpotdt = aux1 - 0.5*pow(b->velR,2) + (prescell-p_state_ec(b->p0, b->R0, Ri))/liq_rho;
  dpotdt = dpotdt/(1. - aux2);
  b->dpotdt = aux1 +  prescell/liq_rho + aux2*dpotdt;
  b->dpotdt = b->dpotdt*liq_rho;
}

static void solve_cluster (FttCell * cell, InteractionData * data, gint nbubbles) {

  gdouble ** m,dij,Ri,Rj,aux,prescell;
  gint nb1,nb2;
  GfsBubble * bi, * bj;
  gdouble radeq;

  gdouble Rb3 = data->Rb3;
  gdouble Rb  = data->Rb;
  gdouble sumRi2v  = data->sumRi2v;
  gdouble sumRi2  = data->sumRi2;

  GfsSimulation * sim = gfs_object_simulation (data->p); 
  gdouble liq_rho = sim->physical_params.alpha ? 1./
    gfs_function_value (sim->physical_params.alpha, cell) : 1.;

  GSList * i,* j;
  i = data->listp->items;

  m = gfs_matrix_new (nbubbles, nbubbles, sizeof(gdouble));
  gdouble * rhs = g_malloc0(nbubbles*sizeof(gdouble));
  nb1 = 0;
  while (i) {
    bi = GFS_BUBBLE(i->data);
    data->b = bi;
    Ri= pow(3.*GFS_PARTICULATE (i->data)->volume/4./M_PI,1./3.);
    prescell = get_pcell(cell,data); 
    nb2 = 0;
    j = data->listp->items;
    while (j) {
      Rj = pow(3.*GFS_PARTICULATE (j->data)->volume/4./M_PI,1./3.);
      bj = GFS_BUBBLE(j->data);
      if (nb1 == nb2) aux=1.;
      else {
        dij = ftt_vector_distance (&(GFS_PARTICLE(i->data)->pos), &(GFS_PARTICLE(j->data)->pos)); 
        aux = MIN(Rj/dij,0.1); //restriction on the bubble interaction term
      }
      radeq = pow(3.*bj->vol_liq/(4.*M_PI),1./3.);
      m[nb1][nb2]=aux-3./2.*Rb*(pow(radeq,2)-pow(Rb,2))/(pow(radeq,3)-Rb3)*pow(Rj,2)/sumRi2;
      nb2 += 1;
      j = j->next;
    }
    radeq = pow(3.*bi->vol_liq/(4.*M_PI),1./3.);
    rhs[nb1] = 3./2.*Rb3/(pow(radeq,3)-Rb3)*pow(sumRi2v/sumRi2,2)*(1.-Rb/radeq) - 0.5*pow(bi->velR,2) \
               + (prescell - p_state_ec(bi->p0, bi->R0 ,Ri) )/liq_rho ; 
    nb1 += 1;
    i = i->next;
  }

  /* solve matrix*/
  i = data->listp->items;
  gdouble * dpotidt = g_malloc(nbubbles*sizeof(gdouble));

  if (gfs_matrix_inverse (m, nbubbles, 1e-10)) {
    for (nb1 = 0; nb1 < nbubbles; nb1++) { 
      dpotidt[nb1] = 0.; 
      for (nb2 = 0; nb2 < nbubbles; nb2++) 
        dpotidt[nb1] += m[nb1][nb2]*rhs[nb2];
    }
  }
  else { /* this may be a degenerate/isolated interface fragment */
    g_warning ("singular matrix in bubble interactions");
    for (nb1 = 0; nb1 < nbubbles; nb1++)  
      dpotidt[nb1] = 0.;
  }

  /*term from finite scale of the cell*/
  gdouble aux2 = 0.;

  i = data->listp->items;
  nb1=0.;
  while (i) {
    Ri= pow(3.*GFS_PARTICULATE (i->data)->volume/4./M_PI,1./3.);
    aux2 += dpotidt[nb1]*pow(Ri,2);
    i = i->next;
    nb1++;
  }

  i = data->listp->items;
  while (i) {
    bi = GFS_BUBBLE(i->data);        
    radeq = pow(3.*bi->vol_liq/(4.*M_PI),1./3.);
    bi->dpotdt = aux2*3./2.*Rb*(pow(radeq,2)-pow(Rb,2))/(pow(radeq,3)-Rb3)/sumRi2;
    bi->dpotdt += 3./2.*Rb3/(pow(radeq,3)-Rb3)*pow(sumRi2v/sumRi2,2)*(1.-Rb/radeq);
  }

  /*bubble-bubble interaction term*/
  i = data->listp->items;
  while (i) {
    bi = GFS_BUBBLE(i->data);
    j = data->listp->items;
    nb2 = 0;
    while (j) {
      if (j != i) {
        dij = ftt_vector_distance (&(GFS_PARTICLE(i->data)->pos), &(GFS_PARTICLE(j->data)->pos)); 
        aux = MIN(Rj/dij,0.1); //restriction on the bubble interaction term
        GFS_BUBBLE(i->data)->dpotdt -= dpotidt[nb2]*aux;
      }
      nb2++;
      j=j->next;
    }
  }

  g_free(dpotidt);
  g_free(rhs);
  gfs_matrix_free (m);
  return;
}

static void interaction_bubbles(FttCell * cell, InteractionData * data)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (data != NULL);

  gdouble Rb,sumRi2,sumRi2v,aux;

  GtsSListContainer * listp = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell,data->v));
  if (listp){
    gint nbubbles = gts_container_size (GTS_CONTAINER (listp));

    if (nbubbles > 10) {
      GSList * j = listp->items;
      while (j) {
        GFS_BUBBLE(j->data)->dpotdt = 0.;
        j = j->next;
      }
      return;
    }

    //Relevant quantities
    Rb = 0; sumRi2 = 0; sumRi2v = 0;
    GSList * j = listp->items;
    while (j) {
      aux = GFS_PARTICULATE (j->data)->volume;
      Rb += aux;
      aux = pow(3.*aux/4./M_PI,2./3.);
      sumRi2 += aux;
      sumRi2v += aux*GFS_BUBBLE(j->data)->velR;
      j = j->next;
    }
    data->Rb3 = 3.*Rb/4./M_PI;
    data->Rb = pow(data->Rb3,1./3.);
    data->sumRi2v=sumRi2v;
    data->sumRi2=sumRi2;
    data->listp = listp;

    if (nbubbles == 1) solve_one_bubble (cell,data);
    else solve_cluster (cell,data,nbubbles);
  }
}


static void destroy_list (FttCell * cell, GfsVariable * v)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (v != NULL);

  GtsSListContainer * listp = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell,v));
  if (listp){
    gts_object_destroy (GTS_OBJECT (listp)); //do I need to destry something else?
  }
}

static gboolean bubble_interactions_event (GfsEvent * event, 
    GfsSimulation * sim) 
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_bubble_interactions_class ())->parent_class)->event)
      (event, sim)) {

    GfsDomain * domain = GFS_DOMAIN (sim); 
    GtsSListContainer * listp ;
    GfsBubbleInteractions * bint = GFS_BUBBLE_INTERACTIONS(event);
    GfsVariable * sv = NULL;
    GfsVariable * p = gfs_variable_from_name (domain->variables, "P");
    sv = gfs_temporary_variable (domain);

    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
        (FttCellTraverseFunc) gfs_cell_reset, sv);
    /* build lists
       I think that the "build list" and "destroy list" should be a more general subroutine... 
       probably in src/particle.c (discuss)*/
    GSList * i = GFS_EVENT_LIST (bint->plist)->list->items; 
    while (i) {
      GfsParticle * part = GFS_PARTICLE (i->data);
      FttCell * cellpart = gfs_domain_locate (domain, part->pos, -1, NULL);
      if (cellpart) {
        if (GFS_DOUBLE_TO_POINTER (GFS_VALUE (cellpart,sv))) listp = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cellpart,sv));
        else { 
          listp = GTS_SLIST_CONTAINER (gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class())));
          GFS_DOUBLE_TO_POINTER (GFS_VALUE (cellpart,sv)) = listp;
        }
        gts_container_add (GTS_CONTAINER (listp), (GtsContainee *) part); // check if it is correct
      }
      i = i->next;
    } 

    /*solve for bubble potentials*/
    InteractionData data = { sv, p }; 
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
        (FttCellTraverseFunc) interaction_bubbles, &data);

    /*destroy lists*/
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
        (FttCellTraverseFunc) destroy_list, sv);

    gts_object_destroy (GTS_OBJECT (sv));
    return TRUE;
  }
  return FALSE;
}

static void bubble_interactions_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_bubble_interactions_class ())->parent_class->read) (o, fp); 
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (object name)");
    return;
  }
  GfsBubbleInteractions * pfield = GFS_BUBBLE_INTERACTIONS (*o);
  GtsObject * object = gfs_object_from_name (GFS_DOMAIN (gfs_object_simulation (*o)), 
      fp->token->str);
  if (object == NULL) {
    gts_file_error (fp, "unknown object '%s'", fp->token->str);
    return;
  }
  if (!GFS_IS_PARTICLE_LIST (object)) {
    gts_file_error (fp, "object '%s' is not a GfsParticleList", fp->token->str);
    return;  
  }
  pfield->plist = GFS_PARTICLE_LIST (object);
  gts_file_next_token (fp);

}

static void bubble_interactions_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_bubble_interactions_class ())->parent_class->write) (o, fp); 
  fprintf (fp, " %s", GFS_EVENT (GFS_BUBBLE_INTERACTIONS (o)->plist)->name);
}

static void bubble_interactions_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = bubble_interactions_event;
  klass->read =  bubble_interactions_read;
  klass->write = bubble_interactions_write;
}

GfsVariableClass * gfs_bubble_interactions_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bubble_interactions_info = {
      "GfsBubbleInteractions",
      sizeof (GfsBubbleInteractions),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) bubble_interactions_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS ( gfs_event_class ()),
        &gfs_bubble_interactions_info);
  }

  return klass;
}

/** \endobject{GfsBubbleInteractions} */

/* GfsFeedBubble: header */

#define GFS_IS_FEED_BUBBLE(obj)         (gts_object_is_from_class (obj,\
					         gfs_feed_bubble_class ()))

GfsEventClass * gfs_feed_bubble_class  (void);

/** \beginobject{GfsFeedBubble} */

static void add_bubble (GfsBubble data_part,
                             GfsParticleList * plist)
{
  guint c;
  GfsSimulation * sim = gfs_object_simulation (plist); 
  GfsEventList * l = GFS_EVENT_LIST (plist); 
  GtsObjectClass * klass = l->klass;
  if (klass == NULL) {
    gfs_error (0, "Unknown particle class\n");
    return;
  }

  /*fixme: better criterion?*/
  if ( data_part.parent.mass < 1.e-20 ) return;
  GtsObject * object = gts_object_new (klass);
  gfs_object_simulation_set (object, sim);
  l->list->items = g_slist_reverse (l->list->items);	
  gts_container_add (GTS_CONTAINER (l->list), GTS_CONTAINEE (object));
  l->list->items = g_slist_reverse (l->list->items);
  GfsEvent * list = GFS_EVENT (l);	
  gfs_event_set (GFS_EVENT (object), 
                 list->start, list->end, list->step, list->istart, list->iend, list->istep);
  GfsParticulate * part = GFS_PARTICULATE (object);
  GfsParticle * p = GFS_PARTICLE (part);
  
  part->vel = data_part.parent.vel;
  p->pos = data_part.parent.parent.pos;
  part->volume = data_part.parent.volume;

  p->id = particle_id(plist);
  part->mass = data_part.parent.mass;
  for (c = 0; c < 3; c++)
    (&part->force.x)[c] = 0.;
  assign_forces ( part , plist->forces);
  
  part->force.x = 0.;
  part->force.y = 0.;
  part->force.z = 0.;

  GfsBubble * bubble = GFS_BUBBLE (object);
  bubble->p0 = data_part.p0;
  bubble->velR = data_part.velR;
  bubble->vol_liq = data_part.vol_liq;
  bubble->dpotdt = data_part.dpotdt;
  bubble->R0 = data_part.R0;
  bubble->params = NULL;
}



static void feed_bubble (GfsDomain * domain, 
	    GfsFeedParticle * feedp, GfsVariable * p, gdouble sigma)
{
  GfsBubble newbubble; 
  newbubble.parent.parent.pos.x = gfs_function_value (feedp->posx, NULL); 
  newbubble.parent.parent.pos.y = gfs_function_value (feedp->posy, NULL); 
  newbubble.parent.parent.pos.z = gfs_function_value (feedp->posz, NULL); 

  FttCell * cell = gfs_domain_locate (domain, newbubble.parent.parent.pos, -1, NULL);    
  if (cell) {
    GfsVariable ** u = gfs_domain_velocity (domain);
    newbubble.parent.vel.x  = gfs_interpolate (cell, newbubble.parent.parent.pos, u[0]);
    newbubble.parent.vel.y  = gfs_interpolate (cell, newbubble.parent.parent.pos, u[1]);
#if FTT_2D
    newbubble.parent.vel.z  = 0.;
#else
    newbubble.parent.vel.z  = gfs_interpolate (cell, newbubble.parent.parent.pos, u[2]);
#endif
    newbubble.parent.volume = gfs_function_value (feedp->vol, cell);
    newbubble.R0 = pow(newbubble.parent.volume*3./(4.*M_PI), 1./3.);
    newbubble.parent.mass = gfs_function_value (feedp->mass, cell);
    newbubble.p0 = gfs_interpolate (cell, newbubble.parent.parent.pos, p) + 2.*sigma/newbubble.R0;
    newbubble.velR = 0.;
    newbubble.vol_liq = 0;
    newbubble.dpotdt  = newbubble.p0;

    add_bubble (newbubble,feedp->plist);
    
  }       
}

static gboolean gfs_feed_bubble_event (GfsEvent * event, GfsSimulation * sim)
{ 
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_feed_particle_class ())->parent_class)->event)
      (event, sim)) {  
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsFeedParticle * feedp = GFS_FEED_PARTICLE (event);
    gint j;
    guint np = gfs_function_value (feedp->np, NULL);

    GSList * i = sim->events->items;
    GfsVariable * p;
    gdouble sigma = 0.;
    gboolean bubbleparams = FALSE;
    while (i) {
	    if (GFS_IS_BUBBLE_PARAMS (i->data)) {
		    bubbleparams = TRUE;
		    p = GFS_BUBBLE_PARAMS (i->data)->p;
		    sigma = GFS_BUBBLE_PARAMS (i->data)->sigma;
	    }
	    i = i->next;
    }
    if (!bubbleparams) 
	p = gfs_variable_from_name (domain->variables, "P");

    for (j = 0; j < np; j++)
      feed_bubble (domain, feedp, p, sigma);
    return TRUE;
  }
  return FALSE;
}

static void gfs_feed_bubble_class_init (GfsEventClass * klass)
{
  GFS_EVENT_CLASS (klass)->event    = gfs_feed_bubble_event;
}

GfsEventClass * gfs_feed_bubble_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_feed_bubble_info = {
      "GfsFeedBubble",
      sizeof (GfsFeedParticle),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_feed_bubble_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_feed_particle_class ()),
				  &gfs_feed_bubble_info);
  }
  return klass;
}

/** \endobject{GfsFeedBubble} */


/* Initialize module */

const gchar gfs_module_name[] = "bubbles";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{ 
  gfs_bubble_params_class ();
  gfs_bubble_class ();
  gfs_bubble_fraction_class ();
  gfs_bubble_fraction_dt_class ();

  gfs_bubble_interactions_class ();
  gfs_feed_bubble_class ();
  return NULL; 
}
