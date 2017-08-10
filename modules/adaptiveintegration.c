/* Gerris - The GNU Flow Solver
 * Copyright (C) 2011-2012 Daniel Fuster
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

#include "simulation.h"
#include "source.h"
#include "adaptive.h"
#include "output.h"
#include "init.h"

/* GfsAdaptiveIntegration: Header */

/* GfsOutputAdaptiveNorm: Header */

typedef struct _GfsOutputAdaptiveNorm        GfsOutputAdaptiveNorm;

struct _GfsOutputAdaptiveNorm {
  /*< private >*/
  GfsOutputScalar parent;
  GfsVariable * v;
  
  /*< public >*/
  GfsFunction * s;
  GfsFunction * w;
};

#define GFS_OUTPUT_ADAPTIVE_NORM(obj)            GTS_OBJECT_CAST (obj,\
					         GfsOutputAdaptiveNorm,\
					         gfs_output_adaptive_norm_class ())
#define GFS_IS_OUTPUT_ADAPTIVE_NORM(obj)         (gts_object_is_from_class (obj,\
						 gfs_output_adaptive_norm_class ()))

GfsOutputClass * gfs_output_adaptive_norm_class  (void);

/**
 * Computing differences to a reference solution.
 * \beginobject{GfsOutputAdaptiveNorm}
 */

static gboolean cell_condition (FttCell * cell, gpointer condition)
{
  return gfs_function_value (condition, cell);
}

static void output_scalar_traverse (GfsOutputScalar * output, 
				    FttTraverseType order,
				    FttTraverseFlags flags,
				    gint max_depth,
				    FttCellTraverseFunc func,
				    gpointer data)
{
  if (output->condition)
    gfs_domain_cell_traverse_condition (GFS_DOMAIN (gfs_object_simulation (output)),
					order, flags, max_depth, 
					func, data,
					cell_condition, output->condition);
  else
    gfs_domain_cell_traverse (GFS_DOMAIN (gfs_object_simulation (output)),
			      order, flags, max_depth, 
			      func, data);
}

static void output_adaptive_norm_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_OUTPUT_ADAPTIVE_NORM (o)->s));
  gts_object_destroy (GTS_OBJECT (GFS_OUTPUT_ADAPTIVE_NORM (o)->w));

  (* GTS_OBJECT_CLASS (gfs_output_adaptive_norm_class ())->parent_class->destroy) (o);
}

static void output_adaptive_norm_read (GtsObject ** o, GtsFile * fp)
{
  GfsOutputAdaptiveNorm * n;

  if (GTS_OBJECT_CLASS (gfs_output_adaptive_norm_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_output_adaptive_norm_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  n = GFS_OUTPUT_ADAPTIVE_NORM (*o);
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
      gts_file_error (fp, "expecting a parameter");
      return;
    }
    else if (!strcmp (fp->token->str, "s")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting `='");
	return;
      }
      gts_file_next_token (fp);
      gfs_function_read (n->s, gfs_object_simulation (*o), fp);
      if (fp->type == GTS_ERROR)
	return;
    }
    else if (!strcmp (fp->token->str, "w")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting `='");
	return;
      }
      gts_file_next_token (fp);
      gfs_function_read (n->w, gfs_object_simulation (*o), fp);
      if (fp->type == GTS_ERROR)
	return;
    }
    else if (!strcmp (fp->token->str, "v")) {
      GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));

      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting `='");
	return;
      }
      gts_file_next_token (fp);
      if (fp->type != GTS_STRING) {
	gts_file_error (fp, "expecting a variable name");
	return;
      }
      if (!(n->v = gfs_domain_get_or_add_variable (domain, fp->token->str, "Error field"))) {
	gts_file_error (fp, "`%s' is a reserved keyword", fp->token->str);
	return;
      }
      gts_file_next_token (fp);
    }
    else {
      gts_file_error (fp, "unknown identifier `%s'", fp->token->str);
      return;
    }
  }
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return;
  }
  fp->scope_max--;
  gts_file_next_token (fp);
}

static void output_adaptive_norm_write (GtsObject * o, FILE * fp)
{
  GfsOutputAdaptiveNorm * n = GFS_OUTPUT_ADAPTIVE_NORM (o);

  if (GTS_OBJECT_CLASS (gfs_output_adaptive_norm_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_output_adaptive_norm_class ())->parent_class->write) 
      (o, fp);
  fputs (" { s = ", fp);
  gfs_function_write (n->s, fp);
  fputs (" w = ", fp);
  gfs_function_write (n->w, fp);
  if (n->v)
    fprintf (fp, " v = %s }", n->v->name);
  else
    fputs (" }", fp);
}

typedef struct {
  GfsSimulation * sim;
  GfsOutputScalar * o;
  GSList * virtualcells;
  GfsNorm * n;
  GfsFunction * w;
  GfsVariable * v;
} ComputeErrorParams;

static void init_vble (FttCell * cell, FttCell * root, GfsVariable * v)
{
  FttVector pos; 
  ftt_cell_pos (cell, &pos);
  GFS_VALUE (cell, v) = gfs_interpolate (root, pos, v);
}

static void add_norm_weighted (FttCell * cell, ComputeErrorParams * cep)
{

  gint np = pow(2., (cep->o->maxlevel-ftt_cell_level(cell)));
  FttVector posroot, pos; 
  ftt_cell_pos (cell, &posroot);
  gdouble L = ftt_cell_size (cell);
  gdouble dx = L/np;
  gdouble vol, pointerr;
  gint i,j,k;
  GFS_VALUE(cell, cep->v) = 0;

  vol = pow(dx, FTT_DIMENSION);
#if FTT_2D
  pos.z = 0.;
  for (i = 0; i < np; i++) {  
    pos.x = posroot.x - L/2. + dx/2. + dx*i;
    for (j = 0; j < np; j++) {  
      pos.y = posroot.y - L/2. + dx/2. + dx*j;
      pointerr = ABS(gfs_interpolate (cell, pos, cep->o->v) -
                 gfs_function_spatial_value (GFS_OUTPUT_ADAPTIVE_NORM (cep->o)->s, &pos));
      gfs_norm_add (cep->n, pointerr, vol*gfs_function_spatial_value (cep->w, &pos));
      GFS_VALUE(cell, cep->v) += pointerr*vol*gfs_function_spatial_value (cep->w, &pos);
    }
  }
#else  /* FTT_3D */
  for (i = 0; i < np; i++) {  
    pos.x = posroot.x - L/2. + dx/2. + dx*i;
    for (j = 0; j < np; j++) {  
      pos.y = posroot.y - L/2. + dx/2. + dx*j;
      for (k = 0; k < np; k++) {  
        pos.z = posroot.z - L/2. + dx/2. + dx*k;
        pointerr = gfs_interpolate (cell, pos, cep->o->v) -
                   gfs_function_spatial_value (GFS_OUTPUT_ADAPTIVE_NORM (cep->o)->s, &pos);
        gfs_norm_add (cep->n, pointerr, vol*gfs_function_spatial_value (cep->w, &pos));
        GFS_VALUE(cell, cep->v) += pointerr*vol*gfs_function_spatial_value (cep->w, &pos);
      }
    }
  }
#endif /* FTT_3D */

  GFS_VALUE(cell, cep->v) /= gfs_cell_volume (cell, GFS_DOMAIN(cep->sim));
}

#ifdef HAVE_MPI
static void norm_reduce (void * i, void * o, 
			 int * len,
			 MPI_Datatype * type)
{
  gdouble * in = (gdouble *) i;
  gdouble * inout = (gdouble *) o;
  g_assert (*len == 5);
  
  inout[0] += in[0];    /* bias */
  inout[1] += in[1];    /* first */
  inout[2] += in[2];    /* second */
  if (in[3] > inout[3]) /* infty */
    inout[3] = in[3];    
  inout[4] += in[4];    /* w */
}


static void domain_norm_reduce (GfsDomain * domain, GfsNorm * n)
{
  if (domain->pid >= 0) {
    double in[5];
    double out[5] = { 0., 0., 0., - G_MAXDOUBLE, 0. };
    MPI_Op op;

    MPI_Op_create (norm_reduce, TRUE, &op);
    in[0] = n->bias; in[1] = n->first; in[2] = n->second; in[3] = n->infty;
    in[4] = n->w;
    MPI_Allreduce (in, out, 5, MPI_DOUBLE, op, MPI_COMM_WORLD);
    MPI_Op_free (&op);
    n->bias = out[0]; n->first = out[1]; n->second = out[2]; n->infty = out[3];
    n->w = out[4];
  }
}
#else /* not HAVE_MPI */
static void domain_norm_reduce (GfsDomain * domain, GfsNorm * n)
{
}
#endif /* not HAVE_MPI */

static GfsNorm domain_norm_error (GfsDomain * domain,
				  ComputeErrorParams * cep,
				  GfsFunction * w,
				  FttTraverseFlags flags,
				  gint max_depth,
				  gboolean (* condition) (FttCell *, gpointer),
				  gpointer cdata)
{
  GfsNorm n;

  g_return_val_if_fail (domain != NULL, n);
  
  gfs_norm_init (&n);
  cep->n = &n;
  cep->w = w;
  if (w)
    gfs_catch_floating_point_exceptions ();
  if (condition)
    gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, flags, max_depth, 
					(FttCellTraverseFunc) add_norm_weighted, cep,
					condition, cdata);
  else
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			      (FttCellTraverseFunc) add_norm_weighted, cep);
  if (w)
    gfs_restore_fpe_for_function (w);
  domain_norm_reduce (domain, &n);
  gfs_norm_update (&n);

  return n;
}



static gboolean gfs_output_adaptive_norm_event (GfsEvent * event, 
					     GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_adaptive_norm_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsOutputScalar * output = GFS_OUTPUT_SCALAR (event);
    GfsOutputAdaptiveNorm * enorm = GFS_OUTPUT_ADAPTIVE_NORM (event);
    GfsVariable * v = enorm->v;
    GfsNorm norm, snorm = { 0., 0., 0., 0. };
    ComputeErrorParams cep = { sim, output, NULL };
    if (v == NULL)
      enorm->v = gfs_temporary_variable (domain);

    cep.v = enorm->v;
    gfs_catch_floating_point_exceptions ();
    gint maxlevel = MAX(output->maxlevel, gfs_domain_depth (domain));
    norm = domain_norm_error (domain, &cep, enorm->w,
				     FTT_TRAVERSE_LEAFS|FTT_TRAVERSE_LEVEL, 
				     maxlevel,
				     output->condition ? cell_condition : NULL,
				     output->condition);

    if (v == NULL) {
      gts_object_destroy (GTS_OBJECT (enorm->v));
      enorm->v = NULL;
    }
    gchar * format;
    if (output->format) {
      gchar * f = output->format;
      format = g_strdup_printf ("%%s time: %s first: %s second: %s infty: %s bias: %s\n",
				f, f, f, f, f);
    }
    else
      format = g_strdup ("%s time: %g first: %10.3e second: %10.3e infty: %10.3e bias: %10.3e\n");

    fprintf (GFS_OUTPUT (event)->file->fp, format,
	     output->name, sim->time.t,
	     norm.first, norm.second, norm.infty, norm.bias);

    g_free (format);
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_adaptive_norm_class_init (GfsOutputClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = output_adaptive_norm_destroy;
  GTS_OBJECT_CLASS (klass)->read = output_adaptive_norm_read;
  GTS_OBJECT_CLASS (klass)->write = output_adaptive_norm_write;
  GFS_EVENT_CLASS (klass)->event = gfs_output_adaptive_norm_event;
}

static void output_adaptive_norm_init (GfsOutputAdaptiveNorm * e)
{
  e->s = gfs_function_new (gfs_function_spatial_class (), 0.);
  e->w = gfs_function_new (gfs_function_spatial_class (), 1.);
}

GfsOutputClass * gfs_output_adaptive_norm_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_adaptive_norm_info = {
      "GfsOutputAdaptiveNorm",
      sizeof (GfsOutputAdaptiveNorm),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_adaptive_norm_class_init,
      (GtsObjectInitFunc) output_adaptive_norm_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_scalar_class ()),
				  &gfs_output_adaptive_norm_info);
  }

  return klass;
}
/* GfsAdaptConvergence: Header */

typedef struct _GfsAdaptConvergence         GfsAdaptConvergence;

struct _GfsAdaptConvergence {
  /*< private >*/
  GfsAdapt parent;
  GfsVariable * v;

  /*< public >*/
  GfsFunction * f;
};

#define GFS_ADAPT_CONVERGENCE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsAdaptConvergence,\
					         gfs_adapt_convergence_class ())
#define GFS_IS_ADAPT_CONVERGENCE(obj)         (gts_object_is_from_class (obj,\
						 gfs_adapt_convergence_class ()))

GfsEventClass * gfs_adapt_convergence_class  (void);

/* \beginobject{GfsOutputAdaptiveNorm} */
/**
 * Adapting cells depending on the value of a function.
 * \beginobject{GfsAdaptConvergence}
 */

static void gfs_adapt_convergence_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_ADAPT_CONVERGENCE (o)->f));

  (* GTS_OBJECT_CLASS (gfs_adapt_convergence_class ())->parent_class->destroy) (o);
}

void print_function (FttCell * cell, GfsAdaptConvergence * a)
{
  g_return_if_fail (cell != NULL);

  gdouble funval = ABS(gfs_function_value (a->f, cell));
  GFS_VALUE(cell, a->v) = funval;

  FttCell * parent = ftt_cell_parent (cell);
  if (parent)
    if ( funval > 1.e-10 )
      GFS_VALUE(cell, a->v) = ABS(funval - gfs_function_value (a->f, parent))/funval;
    else
      GFS_VALUE(cell, a->v) = 0.;
  else
    GFS_VALUE(cell, a->v) = 1.;

}

static gboolean gfs_adapt_convergence_event (GfsEvent * event, 
				       GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_adapt_convergence_class ())->parent_class)->event) 
      (event, sim)) {
    GfsAdaptConvergence * a = GFS_ADAPT_CONVERGENCE (event);
    GfsDomain * domain = GFS_DOMAIN(sim);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) print_function, a);
    return TRUE;
  }
  return FALSE;
}

static void none (FttCell * cell, GfsVariable * v) {}

static void gfs_adapt_convergence_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_adapt_convergence_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  gfs_function_read (GFS_ADAPT_CONVERGENCE (*o)->f, gfs_object_simulation (*o), fp);

  GFS_ADAPT_CONVERGENCE (*o)->v = GFS_ADAPT (*o)->c ? GFS_ADAPT (*o)->c :
    gfs_temporary_variable (GFS_DOMAIN (gfs_object_simulation (*o)));
  GFS_ADAPT_CONVERGENCE (*o)->v->coarse_fine = none;
  GFS_ADAPT_CONVERGENCE (*o)->v->fine_coarse = none;
}

static void gfs_adapt_convergence_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_adapt_convergence_class ())->parent_class->write) (o, fp);
  gfs_function_write (GFS_ADAPT_CONVERGENCE (o)->f, fp);
}

static void gfs_adapt_convergence_class_init (GtsObjectClass * klass)
{
  klass->destroy = gfs_adapt_convergence_destroy;  
  klass->read = gfs_adapt_convergence_read;
  klass->write = gfs_adapt_convergence_write;
  GFS_EVENT_CLASS (klass)->event = gfs_adapt_convergence_event;
}

static gdouble function_cost (FttCell * cell, GfsAdaptConvergence * a)
{
  return GFS_VALUE (cell, a->v);
}

static void gfs_adapt_convergence_init (GfsAdaptConvergence * object)
{
  object->f = gfs_function_new (gfs_function_class (), 0.);
  GFS_ADAPT (object)->cost = (GtsKeyFunc) function_cost;
}

GfsEventClass * gfs_adapt_convergence_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_adapt_convergence_info = {
      "GfsAdaptConvergence",
      sizeof (GfsAdaptConvergence),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_adapt_convergence_class_init,
      (GtsObjectInitFunc) gfs_adapt_convergence_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_adapt_class ()),
				  &gfs_adapt_convergence_info);
  }

  return klass;
}

/** \endobject{GfsAdaptConvergence} */


/* Initialize module */

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "adaptiveintegration";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  gfs_output_adaptive_norm_class ();
  gfs_adapt_convergence_class ();
  return NULL;
} 
