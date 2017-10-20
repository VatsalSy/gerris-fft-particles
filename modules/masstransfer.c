/* Gerris - The GNU Flow Solver
 * Copyright (C) 2012 Gaël Guédon, Riccardo Mereu
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.  
 */

#include <gfs.h>
#include "simulation.h"
#include "source.h"
#include "adaptive.h"
#include "output.h"
#include "solid.h"
#include "mpi_boundary.h"
/**
 * Notes:
 * GfsSourceVOFVelocity
 * GfsMassTransfer
 * GfsMassTransferAxi
 */


/* GfsSourceVOFVelocity: Header */

typedef struct _GfsSourceVOFVelocity        GfsSourceVOFVelocity;

struct _GfsSourceVOFVelocity {
  /*< private >*/
  GfsSourceGeneric parent;

  /*< public >*/
  GfsVariable * c;
  GfsFunction * intensity;
  GfsVariable * nc[FTT_DIMENSION];
};

#define GFS_SOURCE_VOF_VELOCITY(obj)        GTS_OBJECT_CAST (obj,\
                                            GfsSourceVOFVelocity,\
                                            gfs_source_vof_velocity_class ())
#define GFS_IS_SOURCE_VOF_VELOCITY(obj)    (gts_object_is_from_class (obj,\
                                            gfs_source_vof_velocity_class ()))

GfsSourceGenericClass * gfs_source_vof_velocity_class  (void);

/**
 * Add a normal velocity to the advection of a TracerVOF.
 * \beginobject{GfsSourceVOFVelocity}
 */

static void source_vof_velocity_destroy (GtsObject * o)
{
  GfsSourceVOFVelocity * s = GFS_SOURCE_VOF_VELOCITY (o);
  FttComponent c;

  if (s->intensity)
    gts_object_destroy (GTS_OBJECT (s->intensity));

  for (c = 0; c < FTT_DIMENSION; c++) {
    if (s->nc[c])
      gts_object_destroy (GTS_OBJECT (s->nc[c]));
  }

  (* GTS_OBJECT_CLASS (gfs_source_vof_velocity_class ())->parent_class->destroy) (o);
}

static void allocate_normal_vof (GfsSourceVOFVelocity * s)
{
  GfsVariable * t = GFS_VARIABLE (s->c);
  FttComponent c;
  static gchar index[][2] = {"x", "y", "z"};
  for (c = 0; c < FTT_DIMENSION; c++) {
    gchar * name = g_strdup_printf ("%s_n%s", t->name, index[c]);
    gchar * description =
      g_strdup_printf ("%s-component of the normal to the interface computed from grad(%s)",
                  index[c], t->name);
    s->nc[c] = gfs_domain_get_or_add_variable (t->domain, name, description);
    s->nc[c]->units = 0.;
    g_free (name);
    g_free (description);
  }
  gfs_variable_set_vector (s->nc, FTT_DIMENSION);
}

static void source_vof_velocity_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceVOFVelocity * s = GFS_SOURCE_VOF_VELOCITY (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));

  (* GTS_OBJECT_CLASS (gfs_source_vof_velocity_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a variable (TracerVOF)");
    return;
  }
  if ((s->c = gfs_variable_from_name (domain->variables, fp->token->str)) == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  if (!GFS_IS_VARIABLE_TRACER_VOF(s->c)) {
    gts_file_error (fp, "variable `%s' is not a VOF tracer", fp->token->str);
    return;
  }
  if (s->c->sources == NULL)
    s->c->sources = gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ()));
  gts_container_add (s->c->sources, GTS_CONTAINEE (s));
  gts_file_next_token (fp);

  gfs_function_read (s->intensity, domain, fp);
  if (fp->type == GTS_ERROR)
    return;
  gfs_function_set_units (s->intensity, 1.);

  allocate_normal_vof (s);
}

static void source_vof_velocity_write (GtsObject * o, FILE * fp)
{
  GfsSourceVOFVelocity * s = GFS_SOURCE_VOF_VELOCITY (o);
  (* GTS_OBJECT_CLASS (gfs_source_vof_velocity_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s", s->c->name);
  gfs_function_write (s->intensity, fp);
}

static void normal_vof (FttCell * cell, GfsSourceVOFVelocity * s)
{
  GfsVariable ** nc = s->nc;
  GtsVector n = { 0., 0., 0. };
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++)
    n[c] = gfs_center_gradient (cell, c, GFS_VARIABLE (s->c)->i);
  gts_vector_normalize (n);
  for (c = 0; c < FTT_DIMENSION; c++)
    GFS_VALUE (cell, nc[c]) = n[c];
}

static gboolean source_vof_velocity_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_source_vof_velocity_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsSourceVOFVelocity * s = GFS_SOURCE_VOF_VELOCITY (event);
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                (FttCellTraverseFunc) normal_vof, s);
    return TRUE;
  }
  return FALSE;
}

static void source_vof_velocity_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = source_vof_velocity_destroy;
  GTS_OBJECT_CLASS (klass)->read =    source_vof_velocity_read;
  GTS_OBJECT_CLASS (klass)->write =   source_vof_velocity_write;
  GFS_EVENT_CLASS (klass)->event =    source_vof_velocity_event;
}

static void source_vof_velocity_init (GfsSourceVOFVelocity * s)
{
  s->intensity = gfs_function_new (gfs_function_class (), 0.);
}

GfsSourceGenericClass * gfs_source_vof_velocity_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo source_vof_velocity_info = {
      "GfsSourceVOFVelocity",
      sizeof (GfsSourceVOFVelocity),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_vof_velocity_class_init,
      (GtsObjectInitFunc) source_vof_velocity_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = 
      gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_generic_class ()),
                  &source_vof_velocity_info);
  }

  return klass;
}

/** \endobject{GfsSourceVOFVelocity} */


/* GfsMassTransfer: Header */

#define GFS_IS_MASS_TRANSFER(obj)      (gts_object_is_from_class (obj, gfs_mass_transfer_class ()))

GfsSimulationClass * gfs_mass_transfer_class  (void);

/**
 * Solver for simulating mass transfer.
 * \beginobject{GfsMassTransfer}
 */

static void mass_transfer_run (GfsSimulation * sim);

static void mass_transfer_class_init (GfsSimulationClass * klass) 
{
  klass->run = mass_transfer_run;
}

GfsSimulationClass * gfs_mass_transfer_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo mass_transfer_info = {
      "GfsMassTransfer",
      sizeof (GfsSimulation),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) mass_transfer_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()),
                &mass_transfer_info);
  }

  return klass;
}

/* GfsMassTransfer: Object */

typedef struct {
  GfsSourceVOFVelocity * s;
  gdouble coeff;
} CorrectPar;

static void correct_normal_velocity (FttCellFace * face,
            CorrectPar * par)
{
  gdouble du;

  if (GFS_FACE_FRACTION_RIGHT (face) == 0.)
    return;

  du = gfs_function_face_value (par->s->intensity, face)*
    gfs_face_interpolated_value (face, par->s->nc[face->d/2]->i);

  GFS_FACE_NORMAL_VELOCITY_LEFT (face) -= du*par->coeff;

  if (ftt_face_type (face) == FTT_FINE_COARSE)
    du *= 2/FTT_CELLS;
  GFS_FACE_NORMAL_VELOCITY_RIGHT (face) -= du*par->coeff;
}

static void vof_correct_normal_velocities (GfsDomain * domain,
            guint dimension,
            GfsSourceVOFVelocity * s,
            gdouble coeff)
{
  CorrectPar par;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (s != NULL);

  par.s = s;
  par.coeff = coeff;
  gfs_domain_face_traverse (domain, dimension == 2 ? FTT_XY : FTT_XYZ,
              FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
              (FttFaceTraverseFunc) correct_normal_velocity, &par);
}

static void vof_velocity_face_sources (GfsDomain * domain, GfsAdvectionParams * par,
            gdouble coeff)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (par != NULL);

  if (par->v->sources) {
    GSList * i = GTS_SLIST_CONTAINER (par->v->sources)->items;

    while (i) {
      if (GFS_IS_SOURCE_VOF_VELOCITY (i->data)) {
        GfsSourceVOFVelocity * s = i->data;

        vof_correct_normal_velocities (domain, FTT_DIMENSION, s, coeff);
      }
      i = i->next;
    }
  }
}

static void mass_transfer_advance_tracers (GfsSimulation * sim, gdouble dt)
{
  g_return_if_fail (sim != NULL);

  GfsDomain * domain = GFS_DOMAIN (sim);
  GSList * i = domain->variables;
  while (i) {
    if (GFS_IS_VARIABLE_TRACER_VOF (i->data)) {
      GfsVariableTracer * t = i->data;

      t->advection.dt = dt;
      vof_velocity_face_sources (domain, &t->advection, -1.);
      gfs_tracer_vof_advection (domain, &t->advection);
      vof_velocity_face_sources (domain, &t->advection, 1.);
      gfs_domain_variable_centered_sources (domain, i->data, i->data, t->advection.dt);
    }
    i = i->next;
  }

  i = domain->variables;
  while (i) {
    if (GFS_IS_VARIABLE_TRACER (i->data) && !GFS_IS_VARIABLE_TRACER_VOF (i->data)) {
      GfsVariableTracer * t = i->data;

      t->advection.dt = dt;
      gfs_tracer_advection_diffusion (domain, &t->advection, sim->physical_params.alpha);
      gfs_domain_cell_traverse (domain,
                  FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
                  (FttCellTraverseFunc) GFS_VARIABLE (t)->fine_coarse, t);
    }
    i = i->next;
  }
}

static void mass_transfer_run (GfsSimulation * sim)
{
  GfsVariable * p, * pmac, * res = NULL, * g[FTT_DIMENSION], * gmac[FTT_DIMENSION];
  GfsVariable ** gc = sim->advection_params.gc ? g : NULL;
  GfsDomain * domain;
  GSList * i;

  domain = GFS_DOMAIN (sim);

  p = gfs_variable_from_name (domain->variables, "P");
  g_assert (p);
  pmac = gfs_variable_from_name (domain->variables, "Pmac");
  g_assert (pmac);
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++) {
    gmac[c] = gfs_temporary_variable (domain);
    if (sim->advection_params.gc)
      g[c] = gfs_temporary_variable (domain);
    else
      g[c] = gmac[c];
  }
  gfs_variable_set_vector (gmac, FTT_DIMENSION);
  gfs_variable_set_vector (g, FTT_DIMENSION);

  gfs_simulation_refine (sim);
  gfs_simulation_init (sim);

  i = domain->variables;
  while (i) {
    if (GFS_IS_VARIABLE_RESIDUAL (i->data))
      res = i->data;
    i = i->next;
  }

  gfs_simulation_set_timestep (sim);
  if (sim->time.i == 0) {
	//vatsal change hint from masstransfer.c
    gfs_approximate_projection (domain,
                &sim->approx_projection_params,
                sim->advection_params.dt,
                p, sim->physical_params.alpha, res, g, NULL);
    gfs_simulation_set_timestep (sim);
    mass_transfer_advance_tracers (sim, sim->advection_params.dt/2.);
  }
  else if (sim->advection_params.gc)
    gfs_update_gradients (domain, p, sim->physical_params.alpha, g);

  while (sim->time.t < sim->time.end &&
    sim->time.i < sim->time.iend) {
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

    if (sim->advection_params.linear) {
      /* linearised advection */
      gfs_domain_face_traverse (domain, FTT_XYZ,
                  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                  (FttFaceTraverseFunc) gfs_face_reset_normal_velocity, NULL);
      gfs_domain_face_traverse (domain, FTT_XYZ,
                  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                  (FttFaceTraverseFunc) gfs_face_interpolated_normal_velocity,
                  sim->u0);
    }
    else
      gfs_predicted_face_velocities (domain, FTT_DIMENSION, &sim->advection_params);

    gfs_variables_swap (p, pmac);
	//vatsal: check from electrohydro.c
    gfs_mac_projection (domain,
    			&sim->projection_params, 
    			sim->advection_params.dt/2.,
			p, sim->physical_params.alpha, gmac, NULL);
    gfs_variables_swap (p, pmac);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_half_do, sim);

    gfs_centered_velocity_advection_diffusion (domain,
                FTT_DIMENSION,
                &sim->advection_params,
                gmac,
                sim->time.i > 0 || !gc ? gc : gmac,
                sim->physical_params.alpha);
    if (gc) {
      gfs_source_coriolis_implicit (domain, sim->advection_params.dt);
      gfs_correct_centered_velocities (domain, FTT_DIMENSION, sim->time.i > 0 ? gc : gmac, 
                  -sim->advection_params.dt);
    }
    else if (gfs_has_source_coriolis (domain)) {
      gfs_correct_centered_velocities (domain, FTT_DIMENSION, gmac, sim->advection_params.dt);
      gfs_source_coriolis_implicit (domain, sim->advection_params.dt);
      gfs_correct_centered_velocities (domain, FTT_DIMENSION, gmac, -sim->advection_params.dt);
    }

    gfs_domain_cell_traverse (domain,
                FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
                (FttCellTraverseFunc) gfs_cell_coarse_init, domain);
    gfs_simulation_adapt (sim);
	// vatsal: changed using hint from electrohydro.c
    gfs_approximate_projection (domain,
                &sim->approx_projection_params, 
                sim->advection_params.dt, 
                p, sim->physical_params.alpha, res, g, NULL);

    sim->time.t = sim->tnext;
    sim->time.i++;

    gfs_simulation_set_timestep (sim);
    mass_transfer_advance_tracers (sim, sim->advection_params.dt);

    gts_range_add_value (&domain->timestep, gfs_clock_elapsed (domain->timer) - tstart);
    gts_range_update (&domain->timestep);
    gts_range_add_value (&domain->size, gfs_domain_size (domain, FTT_TRAVERSE_LEAFS, -1));
    gts_range_update (&domain->size);
  }
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);  
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gts_object_destroy, NULL);

  for (c = 0; c < FTT_DIMENSION; c++) {
    gts_object_destroy (GTS_OBJECT (gmac[c]));
    if (sim->advection_params.gc)
      gts_object_destroy (GTS_OBJECT (g[c]));
  }
}

/** \endobject{GfsMassTransfer} */


/* GfsMassTransferAxi: Header */

#define GFS_IS_MASS_TRANSFER_AXI(obj)  (gts_object_is_from_class (obj, gfs_mass_transfer_axi_class ()))

GfsSimulationClass * gfs_mass_transfer_axi_class  (void);

/**
 * The axisymmetric solver using the MassTransfer model.
 * \beginobject{GfsMassTransferAxi}
 */

static void mass_transfer_axi_class_init (GfsSimulationClass * klass) 
{
  klass->run = mass_transfer_run;
}

GfsSimulationClass * gfs_mass_transfer_axi_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo mass_transfer_axi_info = {
      "GfsMassTransferAxi",
      sizeof (GfsSimulation),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) mass_transfer_axi_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_axi_class ()), &mass_transfer_axi_info);
  }

  return klass;
}

/** \endobject{GfsMassTransferAxi} */


/* Initialize module */

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "masstransfer";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  gfs_source_vof_velocity_class ();
  gfs_mass_transfer_class ();
  gfs_mass_transfer_axi_class ();
  return NULL;
}

