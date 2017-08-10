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
#include "cartesian.h"
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_wavelet2d.h>
    
/* GfsOutputWavelet: header */

typedef struct _GfsOutputWavelet                     GfsOutputWavelet;

struct _GfsOutputWavelet {
  /*< private >*/
  GfsOutput parent;

  /*< public >*/
  GfsFunction * f;
};

#define GFS_OUTPUT_WAVELET(obj)            GTS_OBJECT_CAST (obj,\
							    GfsOutputWavelet, \
							    gfs_output_wavelet_class ())
#define GFS_IS_OUTPUT_WAVELET(obj)         (gts_object_is_from_class (obj,\
								      gfs_output_wavelet_class ()))

GfsOutputClass * gfs_output_wavelet_class  (void);

/** \beginobject{GfsOutputWavelet} */

typedef struct {
  GfsVariable * v;
  GfsFunction * f;
} WData;

static void init_v (FttCell * cell, WData * wd)
{
  g_return_if_fail (cell != NULL);

  GFS_VALUE(cell, wd->v) = gfs_function_value (wd->f, cell);

}

static void gfs_wavelet_transform_forward ( GfsCartesianGrid * cgd, gsl_wavelet * w )
{
    gint i,j,k;
    gsl_wavelet_workspace * work = gsl_wavelet_workspace_alloc (cgd->n[0] * cgd->n[1] * cgd->n[2]);
    
    /*wavelets for rows (x)*/
    gint k1 = 0;
    for (j = 0; j < cgd->n[1]; j++)
      for (k = 0; k < cgd->n[2]; k++) {
        gsl_wavelet_transform_forward (w, cgd->v + k1*cgd->n[1]*cgd->n[2], 1, cgd->n[0], work);
        k1++;
      }

    /*tranpose x-y*/
    for (k = 0; k < cgd->n[2]; k++) 
      for (i = 0; i < cgd->n[0]; i++) 
        for (j = i+1; j < cgd->n[1]; j++) {
          gdouble tmp = cgd->v[i + cgd->n[0] * ( j + cgd->n[1] * k )];
          cgd->v[i + cgd->n[0] * ( j + cgd->n[1] * k )] = cgd->v[j + cgd->n[1] * ( i + cgd->n[1] * k )];
          cgd->v[j + cgd->n[1] * ( i + cgd->n[1] * k )] = tmp;
        }
 
    /*wavelets for columns (y)*/
    k1 = 0;
    for (j = 0; j < cgd->n[0]; j++)
      for (k = 0; k < cgd->n[2]; k++) {
        gsl_wavelet_transform_forward (w, cgd->v + k1*cgd->n[0]*cgd->n[2], 1, cgd->n[1], work);
        k1++;
      }

    gsl_wavelet_workspace_free (work);

}

static void initialize_cgd ( GfsCartesianGrid * cgd, gint np )
{
    gint i,j;
    cgd->N = 3;
    cgd->x = g_malloc0 (3*sizeof (gdouble *));
    cgd->n = g_malloc (3*sizeof (guint));
    cgd->name = g_malloc0 (3*sizeof (char *));
    cgd->name[0] = g_strdup ("x");
    cgd->name[1] = g_strdup ("y");
    cgd->name[2] = g_strdup ("z");

#if FTT_2D
    cgd->n[0] = cgd->n[1] = np;
    cgd->n[2] = 1;
#else
    cgd->n[0] = cgd->n[1] = cgd->n[2] = np;
#endif

    cgd->v = g_malloc0( sizeof ( gdouble ) * cgd->n[0] * cgd->n[1] * cgd->n[2] );
    gdouble dx = 1./np;

    for (i = 0; i < 2; i++) {
      cgd->x[i] = g_malloc (cgd->n[i]*sizeof (gdouble));
      for (j = 0; j < cgd->n[i]; j++) {
        cgd->x[i][j] = (gdouble) j/(cgd->n[i]-1)*(1. - dx) - (1. - dx)/2.;
      }
    }
#if FTT_2D
    cgd->x[2] = g_malloc (sizeof (gdouble));
    cgd->x[2][0] = 0.;
#else
    cgd->x[2] = g_malloc (np*sizeof (gdouble));
    for (j = 0; j < cgd->n[2]; j++)
      cgd->x[2][j] = ( (gdouble) j/(cgd->n[i]-1) - 0.5 ) ;
#endif

}

static GfsCartesianGrid * cartesian_grid_from_variable ( GfsDomain * domain, GfsVariable * v, gint level)
{
    gint i,j,k;
    GfsCartesianGrid * cgd = gfs_cartesian_grid_new (gfs_cartesian_grid_class ());

    gint np = pow(2,level);
    initialize_cgd(cgd, np);
   
    /* interpolate values*/
    FttVector pos;
    for (i = 0; i < cgd->n[0]; i++)
      for (j = 0; j < cgd->n[1]; j++)
        for (k = 0; k < cgd->n[2]; k++) {
          pos.x = cgd->x[0][i];
          pos.y = cgd->x[1][j];
          pos.z = cgd->x[2][k];
          FttCell * cell = gfs_domain_locate (domain, pos, -1, NULL);
          if (cell)
            cgd->v[i + cgd->n[0] * ( j + cgd->n[1] * k ) ] = gfs_interpolate (cell, pos, v);
          else 
            cgd->v[i + cgd->n[0] * ( j + cgd->n[1] * k ) ] = 0.;
        }

    return cgd;
}

static void write_wavelet  ( FILE * fp, GfsCartesianGrid * cgd)
{
    gint i,j,k,k1;
    fputs ("# 1:i 2:j 3:k 4:coeff\n", fp);
    k1 = 0;
    for (k = 0; k < cgd->n[2]; k++) 
      for (i = 0; i < cgd->n[0]; i++) {
        for (j = 0; j < cgd->n[1]; j++) 
          fprintf(fp, "%i %i %i %g \n", i, j, k, cgd->v[j+k1*cgd->n[0]*cgd->n[2]]);
        k1++;
      }
}

static gboolean output_wavelet_event (GfsEvent * event, 
				      GfsSimulation * sim) 
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_wavelet_class ())->parent_class)->event)
      (event, sim)) {
 
    GfsDomain * domain = GFS_DOMAIN(sim);
    GfsOutputWavelet * a = GFS_OUTPUT_WAVELET (event);
    GfsVariable * v = gfs_temporary_variable (domain);
    WData wd = { v, a->f};
    
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_v, &wd);

    gint level = gfs_domain_depth (domain);
    GfsCartesianGrid * cgd = cartesian_grid_from_variable (domain,v,level);

    gsl_wavelet * w = gsl_wavelet_alloc (gsl_wavelet_daubechies, 4);
    gfs_wavelet_transform_forward ( cgd, w );

    write_wavelet(GFS_OUTPUT (event)->file->fp, cgd);

    gsl_wavelet_free (w);
    gts_object_destroy (GTS_OBJECT (cgd));
    gts_object_destroy (GTS_OBJECT (v));

    return TRUE;
  }
  return FALSE;
}

static void output_wavelet_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_wavelet_class ())->parent_class->read) (o, fp); 
  if (fp->type == GTS_ERROR)
    return;

  gfs_function_read (GFS_OUTPUT_WAVELET (*o)->f, gfs_object_simulation (*o), fp);

}

static void output_wavelet_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_wavelet_class ())->parent_class->write) (o, fp); 
  gfs_function_write (GFS_OUTPUT_WAVELET (o)->f, fp);
}

static void output_wavelet_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = output_wavelet_event;
  klass->read =  output_wavelet_read;
  klass->write = output_wavelet_write;
}

static void output_wavelet_init (GfsOutputWavelet * object)
{
  object->f = gfs_function_new (gfs_function_class (), 0.);
}

GfsOutputClass * gfs_output_wavelet_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsOutputWavelet",
      sizeof (GfsOutputWavelet),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) output_wavelet_class_init,
      (GtsObjectInitFunc) output_wavelet_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsOutputWavelet} */

/**
 * Wavelet analysis.
 * \beginobject{GfsVariableWavelet1}
 */

/* GfsVariableWavelet1: header */

typedef struct _GfsVariableWavelet1                GfsVariableWavelet1;

struct _GfsVariableWavelet1 {
  /*< private >*/
  GfsVariable parent;

  /*< public >*/
  GfsVariable * v, * coeff;
};

#define GFS_VARIABLE_WAVELET1(obj)            GTS_OBJECT_CAST (obj,\
					           GfsVariableWavelet1,\
					           gfs_variable_wavelet1_class ())
#define GFS_IS_VARIABLE_WAVELET1(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_wavelet1_class ()))

GfsVariableClass * gfs_variable_wavelet1_class  (void);

static void none (FttCell * parent, GfsVariable * v)
{
}

static void variable_wavelet1_read (GtsObject ** o, GtsFile * fp)
{
  GfsVariableWavelet1 * v = GFS_VARIABLE_WAVELET1 (*o);
  GfsDomain * domain;

  (* GTS_OBJECT_CLASS (gfs_variable_wavelet1_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (v)");
    return;
  }
  domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(v->v = gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }

  if (GFS_VARIABLE (v)->description)
    g_free (GFS_VARIABLE (v)->description);
  GFS_VARIABLE (v)->description = g_strjoin (" ", "Wavelet coefficients for variable", v->v->name,
					     NULL);
  GFS_VARIABLE (v)->units = v->v->units;
  
  gts_file_next_token (fp);

  if (fp->type == GTS_STRING) {
      v->coeff = gfs_domain_get_or_add_variable (domain, fp->token->str, "Wavelet intensity");
      if (v->coeff) {
	v->coeff->coarse_fine = none;
	v->coeff->fine_coarse = none;
	gts_file_next_token (fp);
        GFS_VARIABLE (v->coeff)->units = v->v->units;
      }
   }
  
}

static void variable_wavelet1_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_wavelet1_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s", GFS_VARIABLE_WAVELET1 (o)->v->name);

  if (GFS_VARIABLE_WAVELET1 (o)->coeff)
    fprintf (fp, " %s", GFS_VARIABLE_WAVELET1 (o)->coeff->name);
    
}

static void wavelet_octree_coefficients (FttCell * cell, GfsVariableWavelet1 * w)
{
  /*fixme: think what to do when there are no child cells (solids)*/
  FttCellChildren child;
  ftt_cell_children(cell,&child);
  
  guint i;
  GFS_VALUE (cell, GFS_VARIABLE (w->coeff)) = 0.;
  for (i = 0; i < (FTT_CELLS-1); i++) 
      GFS_VALUE (cell, GFS_VARIABLE (w->coeff)) += ABS(GFS_VALUE (child.c[i], GFS_VARIABLE (w)));
  
}

static void wavelet_coefficients (FttCell * cell, GfsVariableWavelet1 * w)
{
  /*fixme: think what to do when there are no child cells (solids)*/
  FttCellChildren child;
  ftt_cell_children(cell,&child);
  
#if FTT_2D
  static gint cmatrix[4][4] = {
    {1.414213562, -1.41421356, 0, 0}, 
    {0, 0, 1.41421356, -1.41421356}, 
    {1, 1, -1, -1}, 
    {1, 1, 1, 1}
  };
#else
  static gint cmatrix[8][8] = {
    {2, -2, 0, 0, 0, 0, 0, 0}, 
    {0, 0, 0, 0, 2, -2, 0, 0}, 
    {0, 0, 2, -2, 0, 0, 0, 0}, 
    {0, 0, 0, 0, 0, 0, 2, -2}, 
    {1.41421356, 1.41421356, -1.41421356, -1.41421356, 0, 0, 0, 0}, 
    {0, 0, 0, 0, 1.41421356, 1.41421356, -1.41421356, -1.41421356}, 
    {1, 1, 1, 1, -1, -1, -1, -1}, 
    {1, 1, 1, 1, 1, 1, 1, 1},
  };
#endif

  guint i,j;
  gdouble values[FTT_CELLS];
  for (i = 0; i < FTT_CELLS; i++) {
    if (child.c[i])
      values[i] = GFS_VALUE (child.c[i], w->v);
    else
      values[0] = 0.;
  }

  for (i = 0; i < FTT_CELLS; i++) {
    GFS_VALUE (child.c[i], GFS_VARIABLE (w)) = 0.;
    for (j = 0; j < FTT_CELLS; j++) 
      GFS_VALUE (child.c[i], GFS_VARIABLE (w)) += cmatrix[i][j]*values[j]/pow(2.,0.5*FTT_DIMENSION*ftt_cell_level(child.c[i]));
  }
}

static void wavelet_weight (FttCell * cell, GfsVariable * v)
{
   FttCell * parent = ftt_cell_parent (cell);
   gdouble val = 0.;

   while (parent) {
      val += GFS_VALUE(parent, v);
      parent = ftt_cell_parent (parent);
   }

   GFS_VALUE(cell,v) += val;
}

typedef struct {
  GfsVariable * v1, * v2;
} CopyData;

static gboolean variable_wavelet1_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_wavelet1_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsVariableWavelet1 * w = GFS_VARIABLE_WAVELET1 (event);

    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) w->v->fine_coarse, w->v);
    gfs_domain_bc (domain, FTT_TRAVERSE_NON_LEAFS, -1, w->v);

    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) wavelet_coefficients, w);
    gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, GFS_VARIABLE(w));

    if (w->coeff) {
      gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) wavelet_octree_coefficients, w);
      gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) wavelet_weight, w->coeff);
      gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, w->coeff);     
    }

    return TRUE;
  }
  return FALSE;
}

static void variable_wavelet1_class_init (GtsObjectClass * klass)
{
  klass->read = variable_wavelet1_read;
  klass->write = variable_wavelet1_write;
  GFS_EVENT_CLASS (klass)->event = variable_wavelet1_event;
}

static void variable_wavelet1_init (GfsVariable * v)
{
  v->fine_coarse = none;
  v->coarse_fine = none;
}

GfsVariableClass * gfs_variable_wavelet1_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsVariableWavelet1",
      sizeof (GfsVariableWavelet1),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_wavelet1_class_init,
      (GtsObjectInitFunc) variable_wavelet1_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsVariableWavelet1} */

/**
 * \beginobject{GfsVariableErrorWavelet}
 */

/* GfsVariableErrorWavelet: header */

typedef struct _GfsVariableErrorWavelet                GfsVariableErrorWavelet;

struct _GfsVariableErrorWavelet {
  /*< private >*/
  GfsVariable parent;

  /*< public >*/
  GfsVariableWavelet1 * v1, * v2;
};

#define GFS_VARIABLE_ERROR_WAVELET(obj)            GTS_OBJECT_CAST (obj,\
					           GfsVariableErrorWavelet,\
					           gfs_variable_error_wavelet_class ())
#define GFS_IS_VARIABLE_ERROR_WAVELET(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_error_wavelet_class ()))

GfsVariableClass * gfs_variable_error_wavelet_class  (void);

static void variable_error_wavelet_read (GtsObject ** o, GtsFile * fp)
{
  GfsVariableErrorWavelet * v = GFS_VARIABLE_ERROR_WAVELET (*o);
  GfsDomain * domain;

  (* GTS_OBJECT_CLASS (gfs_variable_error_wavelet_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (v)");
    return;
  }
  domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(v->v1 = gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }

  if (GFS_VARIABLE (v)->description)
    g_free (GFS_VARIABLE (v)->description);
  GFS_VARIABLE (v)->description = g_strjoin (" ", "WaveletError coefficients for variable", NULL);
  
  gts_file_next_token (fp);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (v)");
    return;
  }
  if (!(v->v2 = gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  
  gts_file_next_token (fp);
  
}

static void variable_error_wavelet_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_error_wavelet_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s %s", GFS_VARIABLE(GFS_VARIABLE_ERROR_WAVELET (o)->v1)->name, 
			 GFS_VARIABLE(GFS_VARIABLE_ERROR_WAVELET (o)->v2)->name);
    
}

typedef struct {
  GfsVariable * v1, * v2, * v3;
} SubstractData;

static void substract_data (FttCell * cell, SubstractData * sd)
{
  GFS_VALUE (cell, sd->v1) = ABS(GFS_VALUE (cell, sd->v2)-GFS_VALUE (cell, sd->v3));
}

static gboolean variable_error_wavelet_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_error_wavelet_class ())->parent_class)->event)
      (event, sim)) {

    GfsVariableErrorWavelet * vew = GFS_VARIABLE_ERROR_WAVELET(event);
    SubstractData sd = { GFS_VARIABLE(vew), GFS_VARIABLE(vew->v1), GFS_VARIABLE(vew->v2) };
    GfsDomain * domain = GFS_DOMAIN (sim);

    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) substract_data, &sd);
    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) wavelet_weight, sd.v1);
    gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, sd.v1);     

    return TRUE;
  }
  return FALSE;
}

static void variable_error_wavelet_class_init (GtsObjectClass * klass)
{
  klass->read = variable_error_wavelet_read;
  klass->write = variable_error_wavelet_write;
  GFS_EVENT_CLASS (klass)->event = variable_error_wavelet_event;
}

static void variable_error_wavelet_init (GfsVariable * v)
{
  v->fine_coarse = none;
  v->coarse_fine = none;
}


GfsVariableClass * gfs_variable_error_wavelet_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsVariableErrorWavelet",
      sizeof (GfsVariableErrorWavelet),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_error_wavelet_class_init,
      (GtsObjectInitFunc) variable_error_wavelet_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsVariableErrorWavelet} */

/* GfsVariableDegraded: header */

typedef struct _GfsVariableDegraded                GfsVariableDegraded;

struct _GfsVariableDegraded {
  /*< private >*/
  GfsVariable parent;

  /*< public >*/
  GfsVariable * v;
  guint niter;
};

#define GFS_VARIABLE_DEGRADED(obj)            GTS_OBJECT_CAST (obj,\
					           GfsVariableDegraded,\
					           gfs_variable_degraded_class ())
#define GFS_IS_VARIABLE_DEGRADED(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_degraded_class ()))

GfsVariableClass * gfs_variable_degraded_class  (void);

/**
 * Spatial filtering.
 * \beginobject{GfsVariableDegraded}
 */

static void variable_degraded_read (GtsObject ** o, GtsFile * fp)
{
  GfsVariableDegraded * v = GFS_VARIABLE_DEGRADED (*o);
  GfsDomain * domain;

  (* GTS_OBJECT_CLASS (gfs_variable_degraded_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (v)");
    return;
  }
  domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(v->v = gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting a number (niter)");
    return;
  }
  v->niter = atoi (fp->token->str);
/*  if (v->niter == 0) {
    gts_file_error (fp, "niter must be strictly positive");
    return;
  }*/
  gts_file_next_token (fp);  

  if (GFS_VARIABLE (v)->description)
    g_free (GFS_VARIABLE (v)->description);
  GFS_VARIABLE (v)->description = g_strjoin (" ", "Variable", v->v->name, "degraded", NULL);

  GFS_VARIABLE (v)->units = v->v->units;
}

static void variable_degraded_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_degraded_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s %d", GFS_VARIABLE_DEGRADED (o)->v->name, GFS_VARIABLE_DEGRADED (o)->niter);
}

static void degrade_solution (FttCell * cell, GfsVariableDegraded * vd)
{
   FttCell * parent = ftt_cell_parent (cell);
   gint i;
   if (vd->niter > 0) {
   for (i = 0; i < vd->niter-1; i++)
    parent = ftt_cell_parent (parent);

   if (parent) {
    FttVector p;
    ftt_cell_pos (cell, &p);
    gdouble f[4*(FTT_DIMENSION - 1) + 1];
    gfs_cell_corner_values (parent, vd->v, ftt_cell_level (parent), f);
    GFS_VALUE (cell, GFS_VARIABLE(vd)) = gfs_interpolate_from_corners (parent, p, f);
   } 
   } else
     GFS_VALUE (cell, GFS_VARIABLE(vd)) = GFS_VALUE (parent, vd->v);
}

static void copy_data (FttCell * cell, CopyData * cd)
{
  GFS_VALUE (cell, cd->v1) = GFS_VALUE (cell, cd->v2);
}

static void variable_degraded_event_half (GfsEvent * event, GfsSimulation * sim)
{
    GfsVariable * v = GFS_VARIABLE (event);
    GfsVariableDegraded * vd = GFS_VARIABLE_DEGRADED (event);
    GfsDomain * domain = GFS_DOMAIN (sim);
    CopyData cd = { v, vd->v };

    gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, vd->v);
    gfs_domain_cell_traverse (domain,
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) v->fine_coarse, vd->v);
    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) copy_data, &cd);
    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) degrade_solution, vd);
    gfs_domain_cell_traverse (domain,
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) v->fine_coarse, v);

}

static gboolean variable_degraded_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_degraded_class ())->parent_class)->event)
      (event, sim)) {
    variable_degraded_event_half (event, sim);
    return TRUE;
  }
  return FALSE;
}

static void variable_degraded_class_init (GtsObjectClass * klass)
{
  klass->read = variable_degraded_read;
  klass->write = variable_degraded_write;
  GFS_EVENT_CLASS (klass)->event = variable_degraded_event;
  GFS_EVENT_CLASS (klass)->event_half = variable_degraded_event_half;
}

static void variable_degraded_init (GfsEvent * v)
{
  /* the variable/event may need to be initialised at the start */
  v->start = -1;
}

GfsVariableClass * gfs_variable_degraded_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_degraded_info = {
      "GfsVariableDegraded",
      sizeof (GfsVariableDegraded),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_degraded_class_init,
      (GtsObjectInitFunc) variable_degraded_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_degraded_info);
  }

  return klass;
}

/** \endobject{GfsVariableDegraded}
 *

/* GfsVariableLevel: header */

typedef struct _GfsVariableLevel                GfsVariableLevel;

struct _GfsVariableLevel {
  /*< private >*/
  GfsVariable parent;

  /*< public >*/
  GfsVariable * v;
  guint niter;
};

#define GFS_VARIABLE_LEVEL(obj)            GTS_OBJECT_CAST (obj,\
					           GfsVariableLevel,\
					           gfs_variable_level_class ())
#define GFS_IS_VARIABLE_LEVEL(obj)         (gts_object_is_from_class (obj,\
					     gfs_variable_level_class ()))

GfsVariableClass * gfs_variable_level_class  (void);

/**
 * Spatial filtering.
 * \beginobject{GfsVariableLevel}
 */

static void variable_level_read (GtsObject ** o, GtsFile * fp)
{
  GfsVariableLevel * v = GFS_VARIABLE_LEVEL (*o);
  GfsDomain * domain;

  (* GTS_OBJECT_CLASS (gfs_variable_level_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (v)");
    return;
  }
  domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(v->v = gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting a number (niter)");
    return;
  }
  v->niter = atoi (fp->token->str);
/*  if (v->niter == 0) {
    gts_file_error (fp, "niter must be strictly positive");
    return;
  }*/
  gts_file_next_token (fp);  

  if (GFS_VARIABLE (v)->description)
    g_free (GFS_VARIABLE (v)->description);
  GFS_VARIABLE (v)->description = g_strjoin (" ", "Variable", v->v->name, "level", NULL);

  GFS_VARIABLE (v)->units = v->v->units;
}

static void variable_level_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_level_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s %d", GFS_VARIABLE_LEVEL (o)->v->name, GFS_VARIABLE_LEVEL (o)->niter);
}

static void solution_at_level (FttCell * cell, GfsVariableLevel * vd)
{
  if (vd->niter < ftt_cell_level (cell)) {
   FttCell * parent = ftt_cell_parent (cell);
   while ( ftt_cell_level (parent) > vd->niter)
    parent = ftt_cell_parent (parent);
   GFS_VALUE (cell, GFS_VARIABLE(vd)) = GFS_VALUE (parent,vd->v);
  } 
  else
    GFS_VALUE (cell, GFS_VARIABLE(vd)) = GFS_VALUE (cell,vd->v);
}

static void variable_level_event_half (GfsEvent * event, GfsSimulation * sim)
{
    GfsVariable * v = GFS_VARIABLE (event);
    GfsVariableLevel * vd = GFS_VARIABLE_LEVEL (event);
    GfsDomain * domain = GFS_DOMAIN (sim);
    CopyData cd = { v, vd->v };

    gfs_domain_bc (domain, FTT_TRAVERSE_ALL, -1, vd->v);
    gfs_domain_cell_traverse (domain,
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) v->fine_coarse, vd->v);
    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) copy_data, &cd);
    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) solution_at_level, vd);
    gfs_domain_cell_traverse (domain,
			    FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			    (FttCellTraverseFunc) v->fine_coarse, v);

}

static gboolean variable_level_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_level_class ())->parent_class)->event)
      (event, sim)) {
    variable_level_event_half (event, sim);
    return TRUE;
  }
  return FALSE;
}

static void variable_level_class_init (GtsObjectClass * klass)
{
  klass->read = variable_level_read;
  klass->write = variable_level_write;
  GFS_EVENT_CLASS (klass)->event = variable_level_event;
  GFS_EVENT_CLASS (klass)->event_half = variable_level_event_half;
}

static void variable_level_init (GfsEvent * v)
{
  /* the variable/event may need to be initialised at the start */
  v->start = -1;
}

GfsVariableClass * gfs_variable_level_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_level_info = {
      "GfsVariableLevel",
      sizeof (GfsVariableLevel),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_level_class_init,
      (GtsObjectInitFunc) variable_level_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
				  &gfs_variable_level_info);
  }

  return klass;
}

/** \endobject{GfsVariableLevel}


/* Initialize module */

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "wavelets";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  gfs_output_wavelet_class ();
  gfs_variable_wavelet1_class ();
  gfs_variable_degraded_class ();
  gfs_variable_level_class ();
  gfs_variable_error_wavelet_class ();
  return NULL;
} 
