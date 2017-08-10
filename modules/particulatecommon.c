/* Gerris - The GNU Flow Solver
 * Copyright (C) 2009-2012 National Institute of Water and Atmospheric Research
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

#include <stdlib.h>

#include "particulatecommon.h"
#include "source.h"
#include "refine.h"
#include "adaptive.h"
#include "domain.h"
#include "fluid.h"
#include "solid.h"
#include "vof.h"
#include "mpi_boundary.h"

#ifdef HAVE_MPI
typedef struct {
  GfsParticleList *plist;
  guint *idadd;
} ParticleMpiId;

static void modify_particle_id(GfsEvent * event, ParticleMpiId *part_id)
{

  GfsSimulation *sim = gfs_object_simulation(event);
  GfsDomain *domain = GFS_DOMAIN(sim);

  GfsParticle *p = GFS_PARTICLE(event);

  if(p->id > part_id->plist->maxid)
    p->id = part_id->plist->maxid + part_id->idadd[domain->pid]--;

  return;
}
static void mpi_particle_numbering(GfsDomain *domain, GfsParticleList *plist)
{

  if(domain->pid < 0)
    return;

  guint comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  guint i;

  guint *idadd = g_malloc0(comm_size*sizeof(guint));

  idadd[domain->pid] = plist->idlast - plist->maxid;

  gint tot_adds = 0;
  for(i = 0; i < comm_size; i++){
    gfs_all_reduce(domain, idadd[i], MPI_INT, MPI_MAX);
    tot_adds += idadd[i];
  }

  for(i = 1; i <= domain->pid; i++)
      idadd[i] += idadd[i-1];

  ParticleMpiId part_id;
  part_id.plist = plist;
  part_id.idadd = idadd;

  gts_container_foreach (GTS_CONTAINER (GFS_EVENT_LIST(plist)->list), (GtsFunc) modify_particle_id, &part_id);

  plist->maxid += tot_adds;

  plist->idlast = plist->maxid;

  g_free(idadd);

  return;
}
#endif /*HAVE_MPI*/

/* Forces on the Particle */
static void copy_cell_gfsvariable(FttCell *c, gpointer * data)
{
  GfsVariable * v1 = (GfsVariable *) data[0];
  GfsVariable * v2 = (GfsVariable *) data[1];

  GFS_VALUE(c, v1) = GFS_VALUE(c, v2);
}

static void store_domain_previous_vel(GfsDomain *d, GfsVariable **un)
{
  GfsVariable ** u = gfs_domain_velocity (d);
  FttComponent c;
  for(c = 0; c < FTT_DIMENSION; c++){
    gpointer data[2];

    data[0] = un[c];
    data[1] = u[c];
    gfs_domain_cell_traverse (d, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc)copy_cell_gfsvariable,
			      data);
    gfs_domain_bc (d, FTT_TRAVERSE_LEAFS, -1, un[c]);
  }
}

static FttVector subs_fttvectors (FttVector *a, FttVector *b)
{
  FttVector result;
  FttComponent c;
  for(c = 0; c< FTT_DIMENSION; c++)
    (&result.x)[c]  = (&a->x)[c] - (&b->x)[c];
  return result;
}

/* Same as in source.c used here to obtained viscosity */
static GfsSourceDiffusion * source_diffusion_viscosity (GfsVariable * v)
{
  if (v->sources) {
    GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;

    while (i) {
      GtsObject * o = i->data;

      if (GFS_IS_SOURCE_DIFFUSION (o))
	return GFS_SOURCE_DIFFUSION (o);
      i = i->next;
    }
  }
  return NULL;
}

/* Similar to gfs_vorticity which returns norm of the vorticity */
static void vorticity_vector (FttCell *cell, GfsVariable **v,
			      FttVector *vort)
{
  gdouble size;

  if (cell == NULL) return;
  if (v == NULL) return;

  size = ftt_cell_size (cell);
#if FTT_2D
  vort->x = 0.;
  vort->y = 0.;
  vort->z = (gfs_center_gradient (cell, FTT_X, v[1]->i) -
	     gfs_center_gradient (cell, FTT_Y, v[0]->i))/size;
#else  /* FTT_3D */
  vort->x = (gfs_center_gradient (cell, FTT_Y, v[2]->i) -
	     gfs_center_gradient (cell, FTT_Z, v[1]->i))/size;
  vort->y = (gfs_center_gradient (cell, FTT_Z, v[0]->i) -
	     gfs_center_gradient (cell, FTT_X, v[2]->i))/size;
  vort->z = (gfs_center_gradient (cell, FTT_X, v[1]->i) -
	     gfs_center_gradient (cell, FTT_Y, v[0]->i))/size;
#endif
}

/* GfsForceCoeff: object */

static void gfs_force_coeff_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_force_coeff_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_force_coeff_class ())->parent_class->read)
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsForceCoeff * force = FORCE_COEFF (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));

  force->Uold[0] = gfs_domain_get_or_add_variable(domain, "Un", "x-component of Velocity at previous time step");
  force->Uold[1] = gfs_domain_get_or_add_variable(domain, "Vn", "y-component of Velocity at previous time step");
#if !FTT_2D
  force->Uold[2] = gfs_domain_get_or_add_variable(domain, "Wn", "z-component of Velocity at previous time step");
#endif

  store_domain_previous_vel(domain,force->Uold);

  if (fp->type != '\n') {
    /* fixme: "Rep", "Urelp" etc... should be derived variables not
       straight variables (i.e. there is no need to allocate memory
       for these as they are only used temporarily to compute the
       coefficient) */
    force->re_p = gfs_domain_get_or_add_variable (domain, "Rep",
						  "Particle Reynolds number");
    force->u_rel = gfs_domain_get_or_add_variable (domain, "Urelp",
						   "Particle x - relative velocity");
    force->v_rel = gfs_domain_get_or_add_variable (domain, "Vrelp",
						   "Particle y - relative velocity");
#if !FTT_2D
    force->w_rel = gfs_domain_get_or_add_variable (domain, "Wrelp",
						   "Particle z - relative velocity");
#endif
    force->pdia = gfs_domain_get_or_add_variable (domain, "Pdia",
						  "Particle radii");

    force->coefficient = gfs_function_new (gfs_function_class (), 0.);
    gfs_function_read (force->coefficient, gfs_object_simulation (*o), fp);

  }
}

static void gfs_force_coeff_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_force_coeff_class ())->parent_class->write) (o, fp);
  GfsForceCoeff * force = FORCE_COEFF (o);
  if (force->coefficient)
    gfs_function_write (force->coefficient, fp);
}

static void gfs_force_coeff_destroy (GtsObject * o)
{
  if (FORCE_COEFF (o)->coefficient)
    gts_object_destroy (GTS_OBJECT (FORCE_COEFF (o)->coefficient));

  (* GTS_OBJECT_CLASS (gfs_force_coeff_class ())->parent_class->destroy) (o);
}

static void gfs_force_coeff_class_init (GtsObjectClass * klass)
{
  klass->read = gfs_force_coeff_read;
  klass->write = gfs_force_coeff_write;
  klass->destroy = gfs_force_coeff_destroy;
}

GtsSListContaineeClass * gfs_force_coeff_class (void)
{
  static GtsSListContaineeClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_force_coeff_info = {
      "GfsForceCoeff",
      sizeof (GfsForceCoeff),
      sizeof (GtsSListContaineeClass),
      (GtsObjectClassInitFunc) gfs_force_coeff_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particle_force_class ()),
				  &gfs_force_coeff_info);
  }
  return klass;
}

/* GfsForceInertial: object */
static FttVector compute_inertial_force (GfsParticle * p, GfsParticleForce * inertialforce)
{
  GfsParticulate * particulate = GFS_PARTICULATE (p);
  GfsForceCoeff * coeff = FORCE_COEFF (inertialforce);

  GfsSimulation * sim = gfs_object_simulation (particulate);
  GfsDomain * domain = GFS_DOMAIN (sim);

  FttVector force;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] = 0;

  FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
  if (cell == NULL) return force;

  gdouble size = ftt_cell_size(cell);

  gdouble fluid_rho = sim->physical_params.alpha ? 1./
    gfs_function_value (sim->physical_params.alpha, cell) : 1.;
  GfsVariable ** u = gfs_domain_velocity (domain);

  gdouble viscosity = 0.;
  GfsSourceDiffusion * d = source_diffusion_viscosity (u[0]);
  if (d) viscosity = gfs_diffusion_cell (d->D, cell);

  FttVector fluid_vel;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&fluid_vel.x)[c] = gfs_interpolate (cell, p->pos, u[c]);

  FttVector fluid_veln;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&fluid_veln.x)[c] = gfs_interpolate (cell, p->pos, coeff->Uold[c]);


  if(sim->advection_params.dt > 0.)
    for (c = 0; c < FTT_DIMENSION; c++)
      (&force.x)[c] = fluid_rho*((&fluid_vel.x)[c]-(&fluid_veln.x)[c])/sim->advection_params.dt;
  else
    return force;

  FttComponent c2;
  for(c = 0; c < FTT_DIMENSION; c++)
    for(c2 = 0; c2 < FTT_DIMENSION; c2++)
      (&force.x)[c] += fluid_rho*gfs_center_gradient(cell, c2, u[c]->i)*
	GFS_VALUE(cell, u[c2])/size;

  return force;
}

static void gfs_force_inertial_init (GfsParticleForce * force)
{
  force->force = compute_inertial_force;
}

GtsSListContaineeClass * gfs_force_inertial_class (void)
{
  static GtsSListContaineeClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_force_inertial_info = {
      "GfsForceInertial",
      sizeof (GfsForceCoeff),
      sizeof (GtsSListContaineeClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_force_inertial_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_force_coeff_class ()),
				  &gfs_force_inertial_info);
  }
  return klass;
}

/* GfsForceAddedMass: object */
static FttVector compute_addedmass_force (GfsParticle * p, GfsParticleForce * addedmassforce)
{
  GfsParticulate * particulate = GFS_PARTICULATE (p);
  GfsForceCoeff * coeff = FORCE_COEFF (addedmassforce);

  GfsSimulation * sim = gfs_object_simulation (particulate);
  GfsDomain * domain = GFS_DOMAIN (sim);

  FttVector force;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] = 0;

  FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
  if (cell == NULL) return force;

  force = compute_inertial_force (p, addedmassforce);

  gdouble fluid_rho = sim->physical_params.alpha ? 1./
    gfs_function_value (sim->physical_params.alpha, cell) : 1.;

  gdouble cm = 0.5;
  if (coeff->coefficient) {

    GfsVariable ** u = gfs_domain_velocity (domain);

    gdouble viscosity = 0.;
    GfsSourceDiffusion * d = source_diffusion_viscosity (u[0]);
    if (d) viscosity = gfs_diffusion_cell (d->D, cell);

    FttVector fluid_vel;
    for (c = 0; c < FTT_DIMENSION; c++)
      (&fluid_vel.x)[c] = gfs_interpolate (cell, p->pos, u[c]);

    FttVector relative_vel = subs_fttvectors (&fluid_vel, &particulate->vel);

    gdouble norm_relative_vel = sqrt (relative_vel.x*relative_vel.x +
				      relative_vel.y*relative_vel.y +
				      relative_vel.z*relative_vel.z);
    gdouble dia =  2.*pow(3.0*(particulate->volume)/4.0/M_PI, 1./3.);

    if (viscosity == 0) {
      g_warning ("Viscosity is 0. cannot compute Reynolds number based addedmass force coefficient on particulate: Setting viscosity to 0.001 for Re_p calculation\n");
      viscosity = 0.001;
    }
    gdouble Re = norm_relative_vel*dia*fluid_rho/viscosity;

    GFS_VALUE (cell, coeff->re_p) = Re;
    GFS_VALUE (cell, coeff->pdia) = dia;
    GFS_VALUE (cell, coeff->u_rel) = relative_vel.x;
    GFS_VALUE (cell, coeff->v_rel) = relative_vel.y;
#if !FTT_2D
    GFS_VALUE (cell, coeff->w_rel) = relative_vel.z;
#endif
    cm = gfs_function_value (coeff->coefficient, cell);
  }

  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] *= cm;

  particulate->mass += fluid_rho*particulate->volume*cm;

  return force;
}

static void gfs_force_addedmass_init (GfsParticleForce * force)
{
  force->force = compute_addedmass_force;
}

GtsSListContaineeClass * gfs_force_addedmass_class (void)
{
  static GtsSListContaineeClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_force_addedmass_info = {
      "GfsForceAddedMass",
      sizeof (GfsForceCoeff),
      sizeof (GtsSListContaineeClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_force_addedmass_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_force_coeff_class ()),
				  &gfs_force_addedmass_info);
  }
  return klass;
}

/* GfsForceLift: object */

static FttVector compute_lift_force (GfsParticle * p, GfsParticleForce * liftforce)
{
  GfsParticulate * particulate = GFS_PARTICULATE (p);
  GfsForceCoeff * coeff = FORCE_COEFF (liftforce);

  GfsSimulation * sim = gfs_object_simulation (particulate);
  GfsDomain * domain = GFS_DOMAIN (sim);

  FttVector force;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] = 0;

  FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
  if (cell == NULL) return force;

  gdouble fluid_rho = sim->physical_params.alpha ? 1./
    gfs_function_value (sim->physical_params.alpha, cell) : 1.;
  GfsVariable ** u = gfs_domain_velocity (domain);

  gdouble viscosity = 0.;
  GfsSourceDiffusion * d = source_diffusion_viscosity (u[0]);
  if (d) viscosity = gfs_diffusion_cell (d->D, cell);

  FttVector fluid_vel;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&fluid_vel.x)[c] = gfs_interpolate (cell, p->pos, u[c]);

  FttVector relative_vel = subs_fttvectors (&fluid_vel, &particulate->vel);
  FttVector vorticity;
  vorticity_vector (cell, u, &vorticity);

  gdouble cl = 0.5;
  if (coeff->coefficient) {
    gdouble norm_relative_vel = sqrt (relative_vel.x*relative_vel.x +
				      relative_vel.y*relative_vel.y +
				      relative_vel.z*relative_vel.z);
    gdouble dia =  2.*pow(3.0*(particulate->volume)/4.0/M_PI, 1./3.);
    if (viscosity == 0) {
      g_warning ("Viscosity is 0. cannot compute Reynolds number based lift force coefficient on particulate: Setting viscosity to 0.001 for Re_p calculation\n");
      viscosity = 0.001;
    }
    gdouble Re = norm_relative_vel*dia*fluid_rho/viscosity;

    GFS_VALUE (cell, coeff->re_p) = Re;
    GFS_VALUE (cell, coeff->pdia) = dia;
    GFS_VALUE (cell, coeff->u_rel) = relative_vel.x;
    GFS_VALUE (cell, coeff->v_rel) = relative_vel.y;
#if !FTT_2D
    GFS_VALUE (cell, coeff->w_rel) = relative_vel.z;
#endif
    cl = gfs_function_value (coeff->coefficient, cell);
  }

#if FTT_2D
  force.x = fluid_rho*cl*relative_vel.y*vorticity.z;
  force.y = -fluid_rho*cl*relative_vel.x*vorticity.z;
#else
  force.x = fluid_rho*cl*(relative_vel.y*vorticity.z
			  -relative_vel.z*vorticity.y);
  force.y = fluid_rho*cl*(relative_vel.z*vorticity.x
			  -relative_vel.x*vorticity.z);
  force.z = fluid_rho*cl*(relative_vel.x*vorticity.y
			  -relative_vel.y*vorticity.x);
#endif

  return force;
}

static void gfs_force_lift_init (GfsParticleForce * force)
{
  force->force = compute_lift_force;
}

GtsSListContaineeClass * gfs_force_lift_class (void)
{
  static GtsSListContaineeClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_force_lift_info = {
      "GfsForceLift",
      sizeof (GfsForceCoeff),
      sizeof (GtsSListContaineeClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_force_lift_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_force_coeff_class ()),
				  &gfs_force_lift_info);
  }
  return klass;
}

/* GfsForceDrag: object */

static FttVector compute_drag_force (GfsParticle * p, GfsParticleForce * dragforce)
{
  GfsParticulate * particulate = GFS_PARTICULATE (p);
  GfsForceCoeff * coeff = FORCE_COEFF (dragforce);
  GfsSimulation * sim = gfs_object_simulation (particulate);
  GfsDomain * domain = GFS_DOMAIN (sim);

  FttVector force;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] = 0;

  FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
  if (cell == NULL) return force;

  gdouble fluid_rho = sim->physical_params.alpha ? 1./
    gfs_function_value (sim->physical_params.alpha,cell) : 1.;
  GfsVariable ** u = gfs_domain_velocity (domain);

  gdouble viscosity = 0.;

  GfsSourceDiffusion * d = source_diffusion_viscosity (u[0]);
  if (d) viscosity = gfs_diffusion_cell (d->D, cell);

  FttVector fluid_vel;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&fluid_vel.x)[c] = gfs_interpolate (cell, p->pos, u[c]);

  FttVector relative_vel = subs_fttvectors (&fluid_vel, &particulate->vel);

  gdouble dia = 2.*pow(3.0*(particulate->volume)/4.0/M_PI, 1./3.);
#if !FTT_2D
  gdouble norm_relative_vel = sqrt (relative_vel.x*relative_vel.x +
				    relative_vel.y*relative_vel.y +
				    relative_vel.z*relative_vel.z);
#else
  gdouble norm_relative_vel = sqrt (relative_vel.x*relative_vel.x +
				    relative_vel.y*relative_vel.y);
#endif

  gdouble cd = 0.;
  gdouble Re;
  if (viscosity == 0)
    return force;
  else
    Re = norm_relative_vel*dia*fluid_rho/viscosity;

  if (coeff->coefficient) {
    GFS_VALUE (cell, coeff->re_p) = Re;
    GFS_VALUE (cell, coeff->u_rel) = relative_vel.x;
    GFS_VALUE (cell, coeff->v_rel) = relative_vel.y;
#if !FTT_2D
    GFS_VALUE (cell, coeff->w_rel) = relative_vel.z;
#endif
    GFS_VALUE (cell, coeff->pdia) = dia;
    cd = gfs_function_value (coeff->coefficient, cell);
  }
  else {
    if (Re < 1e-8)
      return force;
    else if (Re < 50.0)
      cd = 16.*(1. + 0.15*pow(Re,0.5))/Re;
    else
      cd = 48.*(1. - 2.21/pow(Re,0.5))/Re;
  }
  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] += 3./(4.*dia)*cd*norm_relative_vel*(&relative_vel.x)[c]*fluid_rho;

  return force;
}

static void gfs_force_drag_init (GfsParticleForce * force)
{
  force->force = compute_drag_force;
}

GtsSListContaineeClass * gfs_force_drag_class (void)
{
  static GtsSListContaineeClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_force_drag_info = {
      "GfsForceDrag",
      sizeof (GfsForceCoeff),
      sizeof (GtsSListContaineeClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_force_drag_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_force_coeff_class ()),
				  &gfs_force_drag_info);
  }
  return klass;
}

/* GfsForceBuoy: object */

static FttVector compute_buoyancy_force (GfsParticle * p, GfsParticleForce * buoyforce)
{
  GfsParticulate * particulate = GFS_PARTICULATE (p);
  GfsSimulation * sim = gfs_object_simulation (particulate);
  GfsDomain * domain = GFS_DOMAIN (sim);

  FttVector force;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] = 0;

  FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
  if (cell == NULL) return force;

  gdouble fluid_rho = sim->physical_params.alpha ? 1./
    gfs_function_value (sim->physical_params.alpha, cell) : 1.;
  GfsVariable ** u = gfs_domain_velocity (domain);

  gdouble g[3];
  for (c = 0; c < FTT_DIMENSION; c++) {
    g[c] = 0.;
    if (u[c]->sources) {
      GSList * i = GTS_SLIST_CONTAINER (u[c]->sources)->items;

      while (i) {
	if (GFS_IS_SOURCE (i->data)) {
	  g[c] += gfs_function_value (GFS_SOURCE ((GfsSourceGeneric *) i->data)->intensity,
				      cell);
	}
	i = i->next;
      }
    }
  }

  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] += (particulate->mass/particulate->volume-fluid_rho)*g[c];

  return force;
}

static void gfs_force_buoy_init (GfsParticleForce * force)
{
  force->force = compute_buoyancy_force;
}

GtsSListContaineeClass * gfs_force_buoy_class (void)
{
  static GtsSListContaineeClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_force_buoy_info = {
      "GfsForceBuoy",
      sizeof (GfsParticleForce),
      sizeof (GtsSListContaineeClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_force_buoy_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particle_force_class ()),
				  &gfs_force_buoy_info);
  }
  return klass;
}

/* GfsParticleForce: object */

static void gfs_particle_force_read (GtsObject ** o, GtsFile * fp)
{
  GtsObjectClass *klass;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsParticleClass)");
    return;
  }

  klass = gfs_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_particle_force_class ())) {
    gts_file_error (fp, "`%s' is not a GfsParticleForce", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
}

static void gfs_particle_force_write (GtsObject * o, FILE * fp)
{
  fprintf (fp, "%s", o->klass->info.name);
}

static void gfs_particle_force_class_init (GtsObjectClass * klass)
{
  GTS_OBJECT_CLASS(klass)->read = gfs_particle_force_read;
  GTS_OBJECT_CLASS(klass)->write = gfs_particle_force_write;
}

GtsSListContaineeClass * gfs_particle_force_class (void)
{
  static GtsSListContaineeClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_particle_force_info = {
      "GfsParticleForce",
      sizeof (GfsParticleForce),
      sizeof (GtsSListContaineeClass),
      (GtsObjectClassInitFunc) gfs_particle_force_class_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_slist_containee_class ()),
				  &gfs_particle_force_info);
  }
  return klass;
}

/* GfsParticulate: Object */

void compute_forces (GfsParticleForce * event, GfsParticulate * p)
{
  FttComponent c;
  FttVector new_force = (event->force) (GFS_PARTICLE (p), event);
  FttVector total_force;

  for ( c = 0 ; c < FTT_DIMENSION; c++)
    (&total_force.x)[c] = (&new_force.x)[c]*p->volume + (&p->force.x)[c];

#if FTT_2D
    (&total_force.x)[2] = 0.;
#endif

  p->force = total_force;
}

static void compute_forces_onfluid (GfsParticleForce * event, GfsParticulate * p)
{
  FttComponent c;
  if (!GFS_IS_FORCE_BUOY(event)) {
    FttVector new_force = (event->force) (GFS_PARTICLE (p), event);
    FttVector total_force;

    for ( c = 0 ; c < FTT_DIMENSION; c++)
      (&total_force.x)[c] = (&new_force.x)[c]*p->volume + (&p->force.x)[c];

    p->force = total_force;
  }
}


static gboolean gfs_particulate_event (GfsEvent * event,
				       GfsSimulation * sim)
{
  GfsParticle * p = GFS_PARTICLE (event);
  GfsParticulate * particulate = GFS_PARTICULATE (event);
  GfsDomain * domain = GFS_DOMAIN (sim);
  FttComponent c;

  if (particulate->forces == NULL)
    (* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_particulate_class ())->parent_class)->event)
      (event, sim);
  else {

  /*if (sim->time.i == 0) {
    GfsVariable ** u = gfs_domain_velocity (domain);
    FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
    for (c = 0; c < FTT_DIMENSION; c++)
      (&particulate->vel.x)[c] = gfs_interpolate (cell, p->pos, u[c]);
  }*/
    /* I get an approximation for the Stokes number */
/*    gdouble St =0. ;
    gdouble viscosity = 0.;
    GfsVariable ** u = gfs_domain_velocity (domain);
    GfsSourceDiffusion * d = source_diffusion_viscosity (u[0]);
    FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
    if (d) {
      viscosity = gfs_diffusion_cell (d->D, cell);
      gdouble dia = 2.*pow(3.0*(particulate->volume)/4.0/M_PI, 1./3.);
      gdouble vnorm = 0.;
      for (c = 0; c < FTT_DIMENSION; c++)
        vnorm += pow(GFS_VALUE(cell,u[c]),2);
      St = particulate->mass*dia*sqrt(vnorm)/(18.*viscosity*particulate->volume);
    }
    else
      St = 1.e10;*/

    FttVector pos = p->pos;
    p->pos_old = pos;
    /* Velocity Verlet Algorithm */
/*    for (c = 0; c < FTT_DIMENSION; c++) {
       (&pos.x)[c] += (&particulate->force.x)[c]*sim->advection_params.dt*sim->advection_params.dt
      	/particulate->mass/2.+ (&particulate->vel.x)[c]*sim->advection_params.dt;
      (&particulate->vel.x)[c] += (&particulate->force.x)[c]*sim->advection_params.dt
      	/(2.*particulate->mass);
    }*/

    /* Compute forces */
    for (c = 0; c < 3; c++)
      (&particulate->force.x)[c] = 0.;

    gts_container_foreach (GTS_CONTAINER (particulate->forces),
    			   (GtsFunc) compute_forces, particulate);

    GfsVariable ** u = gfs_domain_velocity (domain);
    FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
    gdouble vliq = 0.;
    for (c = 0; c < FTT_DIMENSION; c++)
      vliq += pow(gfs_interpolate (cell, p->pos, u[c]),2.);
    vliq = sqrt(vliq);

    for (c = 0; c < FTT_DIMENSION; c++) {
     /* (&particulate->vel.x)[c] +=
    	(&particulate->force.x)[c]*sim->advection_params.dt/(2.*particulate->mass);*/
     (&pos.x)[c] +=
    	(&particulate->vel.x)[c]*sim->advection_params.dt/2.;
     (&particulate->vel.x)[c] +=
    	(&particulate->force.x)[c]*sim->advection_params.dt/particulate->mass;
     (&pos.x)[c] +=
    	(&particulate->vel.x)[c]*sim->advection_params.dt/2.;
    }

    p->pos = pos;
  }
  return TRUE;
}

static void gfs_particulate_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_particulate_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_particulate_class ())->parent_class->read)
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  GfsParticulate * p = GFS_PARTICULATE (*o);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (mass)");
    return;
  }
  p->mass = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (volume)");
    return;
  }
  gdouble L = gfs_object_simulation (*o)->physical_params.L;
  p->volume = atof (fp->token->str);
  p->volume /= pow(L, FTT_DIMENSION);

  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (v.x)");
    return;
  }
  p->vel.x = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (v.y)");
    return;
  }
  p->vel.y = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (v.z)");
    return;
  }
  p->vel.z = atof (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type == GTS_INT || fp->type == GTS_FLOAT) {
    p->force.x = atof (fp->token->str);
    gts_file_next_token (fp);
  }

  if (fp->type == GTS_INT || fp->type == GTS_FLOAT) {
    p->force.y = atof (fp->token->str);
    gts_file_next_token (fp);
  }

  if (fp->type == GTS_INT || fp->type == GTS_FLOAT) {
    p->force.z = atof (fp->token->str);
    gts_file_next_token (fp);
  }

  gfs_simulation_map(gfs_object_simulation(*o), &p->vel);
  gfs_simulation_map(gfs_object_simulation(*o), &p->force);
}

static void gfs_particulate_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_particulate_class ())->parent_class->write) (o, fp);

 GfsParticulate * p = GFS_PARTICULATE (o);
 gdouble L = gfs_object_simulation (o)->physical_params.L;

 gfs_simulation_map_inverse(gfs_object_simulation(o), &p->vel);
 gfs_simulation_map_inverse(gfs_object_simulation(o), &p->force);

 fprintf (fp, " %g %g %g %g %g", p->mass, p->volume*pow(L, FTT_DIMENSION),
	  p->vel.x, p->vel.y, p->vel.z);
 fprintf (fp, " %g %g %g", p->force.x, p->force.y, p->force.z);

 gfs_simulation_map(gfs_object_simulation(o), &p->vel);
 gfs_simulation_map(gfs_object_simulation(o), &p->force);
}

static void gfs_particulate_class_init (GfsEventClass * klass)
{
  klass->event = gfs_particulate_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_particulate_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_particulate_write;
}

GfsEventClass * gfs_particulate_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_particulate_info = {
      "GfsParticulate",
      sizeof (GfsParticulate),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_particulate_class_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particle_class ()),
				  &gfs_particulate_info);
  }
  return klass;
}

/* GfsParticleList: Object */
static void remove_particles_not_in_domain(GfsEvent * event, GfsParticleList *plist)
{

  GfsSimulation *sim = gfs_object_simulation(plist);
  GfsDomain *domain = GFS_DOMAIN(sim);

  GfsParticle *p = GFS_PARTICLE(event);
  FttCell * cell = gfs_domain_locate(domain, p->pos, -1, NULL);

  if(cell==NULL){
       gts_container_remove(GTS_CONTAINER(GFS_EVENT_LIST(plist)->list),GTS_CONTAINEE(event));
       gts_object_destroy(GTS_OBJECT(event));
  }

}

struct _GfsRequest {
  FILE * fp;
  long length;
#ifdef HAVE_MPI
  MPI_Request request[2];
  void * buf;
#endif
};

static gboolean gfs_particle_list_event (GfsEvent * event,
    GfsSimulation * sim)
{

  GfsParticleList *p = GFS_PARTICLE_LIST(event);
  GfsEventList *l = GFS_EVENT_LIST(event);

  gts_container_foreach (GTS_CONTAINER (l->list), (GtsFunc)remove_particles_not_in_domain, p);

  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_particle_list_class ())->parent_class)->event)
      (event, sim)) {


    gfs_particle_bc (GFS_PARTICLE_LIST(event));

#ifdef HAVE_MPI
    mpi_particle_numbering(GFS_DOMAIN(sim), p);
#endif
  }
  else {
    return FALSE;
  }

  GfsDomain * domain = GFS_DOMAIN(sim);
  GSList * i = p->forces->items;
  while (i) {
    GfsParticleForce * pf = GFS_PARTICLE_FORCE(i->data);
    if (GFS_IS_FORCE_INERTIAL(pf)) {
      GfsForceCoeff * coeff = FORCE_COEFF (pf);
      store_domain_previous_vel(domain,coeff->Uold);
    }
    i = i->next;
  }
  return TRUE;

}

void assign_forces (GfsParticulate *particulate, GtsSListContainer *forces)
{
  particulate->forces = forces;
}

static void gfs_particle_list_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_particle_list_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_particle_list_class ())->parent_class->read)
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsParticleList * p = GFS_PARTICLE_LIST (*o);
  GfsEventList * l = GFS_EVENT_LIST (p);
  if (fp->type == '{') {
    fp->scope_max++;
    gts_file_next_token (fp);

    while (fp->type == '\n')
      gts_file_next_token (fp);

    GfsSimulation * sim = gfs_object_simulation (*o);
    GtsObjectClass * klass;
    while (fp->type != '}') {

      if (fp->type != GTS_STRING) {
	gts_file_error (fp, "expecting a keyword (GfsParticleForce)");
	break;
      }
      klass = gfs_object_class_from_name (fp->token->str);

      if (klass == NULL) {
	gts_file_error (fp, "unknown class `%s'", fp->token->str);
	break;
      }
      if (!gts_object_class_is_from_class (klass, gfs_particle_force_class ())) {
	gts_file_error (fp, "'%s' is not a GfsParticleForce", fp->token->str);
	break;
      }

      GtsObject * object = gts_object_new (klass);
      gfs_object_simulation_set (object, sim);

      (* klass->read) (&object, fp);

      if (fp->type == GTS_ERROR) {
	gts_object_destroy (object);
	break;
      }

      while (fp->type == '\n')
	gts_file_next_token (fp);

      gts_container_add (GTS_CONTAINER (p->forces), GTS_CONTAINEE (object));
    }

    if (fp->type != '}') {
      gts_file_error (fp, "expecting a closing brace");
      return;
    }
    fp->scope_max--;
    gts_file_next_token (fp);
  }

  if (p->forces->items != NULL) {
    p->forces->items = g_slist_reverse (p->forces->items);
    gts_container_foreach (GTS_CONTAINER (l->list), (GtsFunc) assign_forces, p->forces);
  }

  if(fp->type == GTS_INT){
    p->idlast = atoi (fp->token->str);
    gts_file_next_token (fp);
  }

  p->first_call = TRUE;

}

static void gfs_particle_list_write (GtsObject * o, FILE * fp)
{

  GfsParticleList *p = GFS_PARTICLE_LIST(o);

  (* GTS_OBJECT_CLASS (gfs_particle_list_class ())->parent_class->write) (o, fp);
  fputs (" {\n", fp);
  GSList * i = p->forces->items;
  while (i) {
    fputs ("    ", fp);
    (* GTS_OBJECT (i->data)->klass->write) (i->data, fp);
    fputc ('\n', fp);
    i = i->next;
  }
  fputc ('}', fp);

  fprintf (fp, " %d", p->idlast);
}

static void gfs_particle_list_init (GtsObject *o){

  GfsParticleList * plist = GFS_PARTICLE_LIST(o);

  plist->forces =
    GTS_SLIST_CONTAINER (gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ())));

  plist->maxid = 0;
}

static void gfs_particle_list_destroy (GtsObject * o)
{
  GfsParticleList * plist = GFS_PARTICLE_LIST(o);

  gts_container_foreach (GTS_CONTAINER (plist->forces), (GtsFunc) gts_object_destroy, NULL);
  gts_object_destroy (GTS_OBJECT (plist->forces));

  (* GTS_OBJECT_CLASS (gfs_particle_list_class ())->parent_class->destroy) (o);
}

static void gfs_particle_list_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_particle_list_event;
  klass->read = gfs_particle_list_read;
  klass->write = gfs_particle_list_write;
  klass->destroy = gfs_particle_list_destroy;
}

GfsEventClass * gfs_particle_list_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_particle_list_info = {
      "GfsParticleList",
      sizeof (GfsParticleList),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_particle_list_class_init,
      (GtsObjectInitFunc) gfs_particle_list_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_list_class ()),
				  &gfs_particle_list_info);
  }
  return klass;
}

/** \beginobject{GfsDropletToParticle} */

typedef struct {
  FttVector pos, vel;
  gdouble volume;
  guint boundary;
} Droplets;

typedef struct {
  GfsVariable * tag, * c, *t;
  Droplets * drops;
  GfsVariable **u;
  guint * sizes;
  guint n, min;
  gdouble resetval;
  gdouble density;
  GfsFunction *fc;
} DropletsPar;

static int greater (const void * a, const void * b)
{
  return *((guint *)a) > *((guint *)b) ? -1 : 1;
}

static void reset_small_fraction (FttCell * cell, DropletsPar * p)
{
  gint i = GFS_VALUE (cell, p->tag);
  if (i > 0 && p->sizes[i - 1] < p->min)
    GFS_VALUE (cell, p->c) = p->resetval;
}

static void compute_droplet_properties (FttCell * cell, DropletsPar * p)
{
  g_assert(cell!=NULL);
  gint i = GFS_VALUE (cell, p->tag);

  if(i > 0) {
    FttDirection d ;
    for(d = 0; d < FTT_NEIGHBORS; d++){
      FttCell * neighbor = ftt_cell_neighbor (cell, d);
      if(neighbor && GFS_CELL_IS_BOUNDARY(neighbor)) {
	p->drops[i - 1].boundary = 1;
	return;
      }
      else
	p->drops[i - 1].boundary = 0;
    }
  }
  else
    return;

  gdouble h = ftt_cell_size (cell), vol;
  FttVector pos;
  ftt_cell_pos (cell, &pos);
  GfsVariable ** u = p->u;

  if (i > 0) {
    p->sizes[i - 1]++;
    vol = pow (h, FTT_DIMENSION);
    p->drops[i-1].volume += vol*GFS_VALUE (cell, p->c);
    FttComponent c;
    for(c = 0; c < FTT_DIMENSION; c++){
      (&(p->drops[i-1].pos.x))[c] +=  (&pos.x)[c];
      (&(p->drops[i-1].vel.x))[c] += GFS_VALUE (cell,u[c]);
    }
  }
}

gint particle_id(GfsParticleList * plist)
{
  return ++plist->idlast;
}


static void add_particulate (GfsParticulate data_part,
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
  if ( data_part.mass < 1.e-20 ) return;
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

  part->vel = data_part.vel;
  p->pos = data_part.parent.pos;
  part->volume = data_part.volume;

  p->id = particle_id(plist);
  part->mass = data_part.mass;
  for (c = 0; c < 3; c++)
    (&part->force.x)[c] = 0.;
  assign_forces ( part , plist->forces);

  part->force.x = 0.;
  part->force.y = 0.;
  part->force.z = 0.;

}

static void convert_droplets (GfsDomain * domain,
			      DropletsPar * pars, GfsDropletToParticle * d)
{
  GfsSimulation * sim = gfs_object_simulation (d);
  guint i;

  pars->sizes = g_malloc0 (pars->n*sizeof (guint));
  pars->drops = g_malloc0 (pars->n*sizeof (Droplets));

  FttComponent c;
  /* Initialize drops */
  for (i = 0; i < pars->n; i++){
    pars->drops[i].volume = 0.;
    pars->sizes[i] = 0;
    pars->drops[i].boundary = 0;
    for(c = 0; c < FTT_DIMENSION; c++) {
      (&(pars->drops[i].pos.x))[c] = 0.;
      (&(pars->drops[i].vel.x))[c] = 0.;
    }
  }

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) compute_droplet_properties, pars);

#ifdef HAVE_MPI
  if (domain->pid >= 0) {
    guint * sizes = g_malloc0 (pars->n*sizeof (guint));
    MPI_Allreduce (pars->sizes, sizes, pars->n, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    g_free (pars->sizes);
    pars->sizes = sizes;
  }
#endif

  if (d->min >= 0)
    pars->min = d->min;
  else {
    guint * tmp = g_malloc (pars->n*sizeof (guint));
    memcpy (tmp, pars->sizes, pars->n*sizeof (guint));
    qsort (tmp, pars->n, sizeof (guint), greater);
    g_assert (-1 - d->min < pars->n);
    /* fixme: this won't work for parallel jobs */
    pars->min = tmp[-1 - d->min];
    g_free (tmp);
  }

  for (i = 0; i < pars->n; i++) {
    if (pars->drops[i].boundary) continue;
    if (pars->sizes[i] < pars->min && pars->sizes[i] > 1){
      GfsParticulate newpart;
      for (c = 0; c < FTT_DIMENSION; c++) {
      	(&newpart.parent.pos.x)[c] = (&pars->drops[i].pos.x)[c]/pars->sizes[i];
	(&newpart.vel.x)[c] = (&pars->drops[i].vel.x)[c]/pars->sizes[i];
      }
      FttCell * cell = gfs_domain_locate (domain, newpart.parent.pos, -1, NULL);
      if (cell) {
	newpart.volume = pars->drops[i].volume;
	newpart.mass = sim->physical_params.alpha ? 1./
	  gfs_function_value (sim->physical_params.alpha, cell) : 1.;
	newpart.mass *= newpart.volume;
        add_particulate (newpart,d->plist);
      }
    }
  }
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) reset_small_fraction, pars);


  g_free (pars->drops);
  g_free (pars->sizes);

}

/* GfsDropletToParticle: object */

typedef struct {
  GfsVariable * v;
  GfsFunction * fc;
} compute_v_data;

static void compute_v (FttCell * cell, compute_v_data * d)
{
  GFS_VALUE (cell, d->v) = gfs_function_value (d->fc, cell);
}

static gboolean gfs_droplet_to_particle_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsDropletToParticle *d = DROPLET_TO_PARTICLE (event);
    GfsVariable * v = d->fc ? gfs_function_get_variable (d->fc) : d->c;
    DropletsPar p ;

    p.resetval = d->resetwith;
    p.tag = gfs_temporary_variable (domain);
    p.u = gfs_domain_velocity (domain);
    p.density = d->density;
    p.t = d->c;

    if (v){
      p.c = v;
      p.n = gfs_domain_tag_droplets (domain, p.c, p.tag);

      if (p.n > 0 && -d->min < (gint) p.n){
	p.c = d->c;
	convert_droplets (domain, &p, d);
      }
    }
    else {
      v = gfs_temporary_variable (domain);
      compute_v_data cvd = { v , d->fc };
      gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
				(FttCellTraverseFunc) compute_v, &cvd);
      p.c = v;

      p.n = gfs_domain_tag_droplets (domain, p.c, p.tag);
      if (p.n > 0 && -d->min < (gint) p.n){
	p.c = d->c;
	convert_droplets (domain, &p, d);
      }
      gts_object_destroy (GTS_OBJECT (v));
    }

    gts_object_destroy (GTS_OBJECT (p.tag));
    return TRUE;
  }
  return FALSE;
}

static void gfs_droplet_to_particle_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (object name)");
    return;
  }

  GfsDropletToParticle * r = DROPLET_TO_PARTICLE (*o);

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
  gts_file_next_token (fp);

  r->plist = GFS_PARTICLE_LIST (object);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (variable)");
    return;
  }

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (r));

  if ((r->c = gfs_variable_from_name (domain->variables, fp->token->str)) == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_INT, "min",  TRUE},
      {GTS_DOUBLE, "reset",    TRUE},
      {GTS_DOUBLE, "density",   TRUE},
      {GTS_NONE}
    };

    var[0].data = &r->min;
    var[1].data = &r->resetwith;
    var[2].data = &r->density;

    gts_file_assign_variables (fp, var);
  }

  if (fp->type != '\n') {
    r->fc = gfs_function_new (gfs_function_class (), 0.);
    gfs_function_read (r->fc, gfs_object_simulation (r), fp);
  }
}

static void gfs_droplet_to_particle_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class->write) (o, fp);

  GfsDropletToParticle * r = DROPLET_TO_PARTICLE(o);

  fprintf (fp, " %s", GFS_EVENT (r->plist)->name);

  fprintf (fp, " %s { min = %d reset = %g density = %g } ",
	   r->c->name, r->min, r->resetwith, r->density);
  if (r->fc)
    gfs_function_write (r->fc, fp);
}

static void gfs_droplet_to_particle_destroy (GtsObject * o)
{
  GfsDropletToParticle * drops = DROPLET_TO_PARTICLE (o);
  if (drops->fc)
    gts_object_destroy (GTS_OBJECT (drops->fc));

  (* GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class->destroy) (o);
}

static void gfs_droplet_to_particle_init (GfsDropletToParticle * r)
{
  r->resetwith = 0.;
  r->min = 20;
  r->density = 1.;
}

static void gfs_droplet_to_particle_class_init (GfsEventClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_droplet_to_particle_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_droplet_to_particle_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_droplet_to_particle_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_droplet_to_particle_destroy;
}

GfsEventClass * gfs_droplet_to_particle_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_droplet_to_particle_info = {
      "GfsDropletToParticle",
      sizeof (GfsDropletToParticle),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_droplet_to_particle_class_init,
      (GtsObjectInitFunc) gfs_droplet_to_particle_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_droplet_to_particle_info);
  }
  return klass;
}

/** \endobject{GfsDropletToParticle} */

/** \beginobject{GfsParticleToDroplet} */

typedef struct {
  guint maxlevel;
  GfsDomain * domain;
  GfsGenericSurface * surface;
  gboolean check;
} RefineCut;

static void save_p_solid (FttCell * cell, gpointer * data)
{
  GfsVariable * c = data[0];
  GfsVariable * save_prev = data[3];

  GFS_VALUEI (cell, save_prev->i) = GFS_VALUEI (cell, c->i);
  GFS_DOUBLE_TO_POINTER (GFS_VALUEI (cell, c->i)) = GFS_STATE (cell)->solid;
  GFS_STATE (cell)->solid = NULL;
}

static void restore_p_solid (FttCell * cell, gpointer * data)
{
  GfsVariable * c = data[0];
  gboolean * not_cut = data[1];
  GfsVariable * status = data[2];

  GfsSolidVector * solid = GFS_STATE (cell)->solid;

  GFS_STATE (cell)->solid = GFS_DOUBLE_TO_POINTER (GFS_VALUEI (cell, c->i));

  if (solid) {
    GFS_VALUEI (cell, c->i) =  solid->a;
    g_free (solid);
    *not_cut = FALSE;
  }
  else if (GFS_VALUEI (cell, status->i) == 0.) {
    /* fixme: this can fail for non-contiguous domains (e.g. non-connected GfsBoxes) */
    g_assert (*not_cut);
    GFS_VALUEI (cell, c->i) = 0.;
  }
  else {
    g_assert (GFS_VALUEI (cell, status->i) == 1. || GFS_VALUEI (cell, status->i) == 2.);
    GFS_VALUEI (cell, c->i) =  GFS_VALUEI (cell, status->i) - 1.;
  }
}

static void add_to_prev_void (FttCell * cell, gpointer * data)
{
  GfsVariable * c = data[0];
  GfsVariable * save_prev = data[3];
  GfsParticulate *p = GFS_PARTICULATE((GfsParticle *)data[4]);
  GfsDomain *domain = data[5];
  GfsVariable **u = gfs_domain_velocity(domain);

  if(GFS_VALUE (cell, c) > 0.){
    guint j;
    for(j = 0; j < FTT_DIMENSION; j++)
      GFS_VALUE(cell, u[j]) = (&p->vel.x)[j];
  }

  GFS_VALUEI (cell, c->i) += GFS_VALUEI (cell, save_prev->i);

  if(GFS_VALUE (cell, c) > 1.0)
    GFS_VALUEI (cell, c->i) = 1.0;
  else if (GFS_VALUE (cell, c) < 0.0)
      GFS_VALUEI (cell, c->i) = 0.0;
}

static void gfs_domain_assign_fraction (GfsDomain * domain,
			       GfsGenericSurface * s,
			       GfsVariable * c,  GfsParticle *p)
{
  gboolean not_cut = TRUE;
  gpointer data[6];
  GfsVariable * status, *save_prev;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (s != NULL);
  g_return_if_fail (c != NULL);

  status = gfs_temporary_variable (domain);

  save_prev = gfs_temporary_variable (domain);
  data[0] = c;
  data[1] = &not_cut;
  data[2] = status;
  data[3] = save_prev;
  data[4] = p;
  data[5] = domain;

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) save_p_solid, data);

  GfsSolid tmp;
  tmp.s = s;
  GSList * l = g_slist_prepend (NULL, &tmp);

  gfs_domain_init_solid_fractions (domain, l, FALSE, NULL, NULL, status);

  if(l)  g_slist_free (l);

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) restore_p_solid, data);

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
  			    (FttCellTraverseFunc) add_to_prev_void, data);

  gts_object_destroy (GTS_OBJECT (status));
  gts_object_destroy (GTS_OBJECT (save_prev));

}


static void generate_surface(GfsSurface *surface, GfsSurface *shape, gdouble rad, GfsParticle *p)
{

  memcpy(surface, shape, sizeof(GfsSurface));

  surface->translate[0] = p->pos.x;
  surface->translate[1] = p->pos.y;
  surface->translate[2] = p->pos.z;

  GtsMatrix * m = gts_matrix_translate (NULL, surface->translate);

  surface->scale[0] = rad;
  surface->scale[1] = rad;
  surface->scale[2] = rad;

  GtsMatrix * ms = gts_matrix_scale (NULL, surface->scale);
  if (surface->m)
    gts_matrix_destroy (surface->m);
  surface->m = gts_matrix_product (m, ms);
  gts_matrix_destroy (m);
  gts_matrix_destroy (ms);

  if (surface->s) {
    gts_surface_foreach_vertex (surface->s, (GtsFunc) gts_point_transform, surface->m);
    gts_matrix_destroy (surface->m);
    surface->m = NULL;
    if (surface->flip)
      gts_surface_foreach_face (surface->s, (GtsFunc) gts_triangle_revert, NULL);
  }
  else {
    GtsMatrix * i = gts_matrix_inverse (surface->m);
    gts_matrix_destroy (surface->m);
    surface->m = i;
  }
}

static double compute_xyz (FttCell * cell, GfsParticleToDroplet * d)
{
  if(d->fc)
    return gfs_function_value (d->fc, cell);
  else
    return -1.0;
}

#define PROXIMITY_FRAC 0.25

typedef struct {
  GfsVariable *v;
  gboolean check;
} ProximityCheck;


static void check_cut_cells(FttCell *cell, GfsGenericSurface *s, ProximityCheck *proximity)
{
  GfsVariable *v = proximity->v;

  if(!proximity->check){
    if(GFS_VALUE(cell, v) > PROXIMITY_FRAC){
      proximity->check = TRUE;
    }
  }
}

static gboolean check_proximity(GfsVariable *v, GfsGenericSurface *surface)
{

  GfsDomain *domain = GFS_DOMAIN(gfs_object_simulation(v));
  ProximityCheck proximity;
  proximity.v = v;
  proximity.check = FALSE;

  gfs_domain_traverse_cut (domain,
			   surface,
			   FTT_PRE_ORDER,
			   FTT_TRAVERSE_LEAFS,
			   (FttCellTraverseCutFunc) check_cut_cells,
			   &proximity);

  return proximity.check;
}

static void refine_implicit_p_cell (FttCell * cell, RefineCut * p)
{
  guint maxlevel = p->maxlevel;

  if (ftt_cell_level (cell) < maxlevel && (gfs_cell_is_cut (cell, p->surface, FALSE, maxlevel))){
    ftt_cell_refine_single (cell, p->domain->cell_init, p->domain->cell_init_data);
    p->check = TRUE;
  }
}

static void convert_particles(GfsEvent * event, GfsParticleToDroplet *d)
{
  GfsSimulation *sim = gfs_object_simulation(d);
  GfsDomain *domain = GFS_DOMAIN(sim);

  GfsParticulate * p = GFS_PARTICULATE(event);
  GfsParticle *particle = GFS_PARTICLE(event);

  gdouble rad = pow(p->volume/M_PI,1./2.);
#if !FTT_2D
  rad = pow(3.0*(p->volume)/4.0/M_PI, 1./3.);
#endif

  /*Fix me: Create surface only when you are sure of converting the particle
            One can use simple radius based check to find the proximity criterion (probably Kernel based)*/

  GfsGenericSurface * surface = GFS_GENERIC_SURFACE (gts_object_new
					      (GTS_OBJECT_CLASS (gfs_surface_class ())));
  generate_surface (GFS_SURFACE(surface), d->shape, rad, particle);

  FttCell * cellpart = gfs_domain_locate (domain, particle->pos, -1, NULL);

  g_assert(cellpart!=NULL);

  gboolean check = FALSE;
  if(rad > ftt_cell_size(cellpart))
    check = check_proximity(d->c, surface);
  else
    if(GFS_VALUE(cellpart, d->c) > PROXIMITY_FRAC)
      check = TRUE;

  if(cellpart && (compute_xyz(cellpart, d) < 0 || check) ){

    RefineCut prefine;
    prefine.domain = domain;
    prefine.maxlevel = d->maxlevel;
    prefine.surface = surface;
    prefine.check = TRUE;
    while(prefine.check){
      prefine.check = FALSE;
      gfs_domain_cell_traverse (domain,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) refine_implicit_p_cell, &prefine);
    }
      gfs_domain_assign_fraction (domain,
				  surface,
				  d->c, particle);

      gts_container_remove(GTS_CONTAINER(GFS_EVENT_LIST(d->plist)->list),GTS_CONTAINEE(event));
      gts_object_destroy (GTS_OBJECT (event));

  }

  if (GFS_SURFACE(surface)->s)
    gts_object_destroy (GTS_OBJECT (GFS_SURFACE(surface)->s));

  if (GFS_SURFACE(surface)->m)
    gts_matrix_destroy (GFS_SURFACE(surface)->m);
}

static gboolean gfs_particle_to_droplet_event (GfsEvent * event,
					 GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_particle_to_droplet_class ())->parent_class)->event)
      (event, sim)) {

    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsParticleToDroplet *d = GFS_PARTICLE_TO_DROPLET(event);
    GfsParticleList *plist = d->plist;

    GfsEventList *l = GFS_EVENT_LIST(plist);

    gts_container_foreach (GTS_CONTAINER (l->list), (GtsFunc)remove_particles_not_in_domain, plist);
    gts_container_foreach (GTS_CONTAINER (l->list), (GtsFunc)convert_particles, d);

    gfs_domain_reshape (domain, gfs_domain_depth(domain));

    return TRUE;
  }
  return FALSE;
}

static void gfs_particle_to_droplet_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_particle_to_droplet_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (object name)");
    return;
  }

  GfsParticleToDroplet * p = GFS_PARTICLE_TO_DROPLET (*o);
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
  gts_file_next_token (fp);

  p->plist = GFS_PARTICLE_LIST (object);

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (p));

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (variable)");
    return;
  }

  if ((p->c = gfs_variable_from_name (domain->variables, fp->token->str)) == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_INT, "maxlevel",  TRUE},
      {GTS_NONE}
    };

    var[0].data = &p->maxlevel;
    gts_file_assign_variables (fp, var);
  }

  while(fp->type == '\n')
    gts_file_next_token (fp);

  p->shape = GFS_SURFACE (gts_object_new (GTS_OBJECT_CLASS (gfs_surface_class ())));
  gfs_generic_surface_read (GFS_GENERIC_SURFACE(p->shape), gfs_object_simulation (*o), fp);

  if (fp->type != '\n') {
    p->fc = gfs_function_new (gfs_function_class (), 0.);
    gfs_function_read (p->fc, gfs_object_simulation (p), fp);
  }

}

static void gfs_particle_to_droplet_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_particle_to_droplet_class ())->parent_class->write) (o, fp);

  GfsParticleToDroplet * r = GFS_PARTICLE_TO_DROPLET(o);

  fprintf (fp, " %s", GFS_EVENT (r->plist)->name);

  fprintf (fp, " %s { maxlevel = %d } ",
	   r->c->name, r->maxlevel);
  if (r->shape) {
    fputs (" ( ", fp);
    gfs_function_write (r->shape->f, fp);
    fputs (" ) ", fp);
  }

  if (r->fc)
    gfs_function_write (r->fc, fp);

}

static void particle_to_droplet_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_particle_to_droplet_event;
  klass->read  = gfs_particle_to_droplet_read;
  klass->write = gfs_particle_to_droplet_write;
}

GfsEventClass * gfs_particle_to_droplet_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_particle_to_droplet_info = {
      "GfsParticleToDroplet",
      sizeof (GfsParticleToDroplet),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) particle_to_droplet_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
                                  &gfs_particle_to_droplet_info);
  }

  return klass;
}

/** \endobject{GfsParticleToDroplet} */

/** \beginobject{GfsParticulateField} */

static void voidfraction_from_particles (FttCell * cell, GfsVariable * v, GfsParticulate * part)
{
  GFS_VALUE (cell, v) += part->volume/ftt_cell_volume (cell);
}

static gboolean particulate_field_event (GfsEvent * event,
					 GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_particulate_field_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsVariable * v = GFS_VARIABLE (event);
    GfsParticulateField * pfield = GFS_PARTICULATE_FIELD (v);

    /* Reset variable */
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_reset, v);
    /* Loop over the list of particles in the selected object */
    GSList * i = GFS_EVENT_LIST (pfield->plist)->list->items;
    while (i) {
      FttCell * cellpart = gfs_domain_locate (domain, GFS_PARTICLE (i->data)->pos, -1, NULL);
      if (cellpart)
	pfield->voidfraction_func (cellpart, v, i->data);
      i = i->next;
    }
    return TRUE;
  }
  return FALSE;
}

static void particulate_field_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_particulate_field_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (object name)");
    return;
  }

  GfsParticulateField * pfield = GFS_PARTICULATE_FIELD (*o);
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
  gts_file_next_token (fp);

  pfield->plist = GFS_PARTICLE_LIST (object);
}

static void particulate_field_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_particulate_field_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s", GFS_EVENT (GFS_PARTICULATE_FIELD (o)->plist)->name);
}

static void particulate_field_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = particulate_field_event;
  klass->read =  particulate_field_read;
  klass->write = particulate_field_write;
}

static void particulate_field_init (GfsVariable * v)
{
  v->units = -FTT_DIMENSION;
  GFS_PARTICULATE_FIELD (v)->voidfraction_func = voidfraction_from_particles;
}

GfsVariableClass * gfs_particulate_field_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_particulate_field_info = {
      "GfsParticulateField",
      sizeof (GfsParticulateField),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) particulate_field_class_init,
      (GtsObjectInitFunc) particulate_field_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()),
                                  &gfs_particulate_field_info);
  }

  return klass;
}
/** \endobject{GfsParticulateField} */

/** \beginobject{GfsSourceParticulate} */

static gdouble source_particulate_value(GfsSourceGeneric * s,
                                           FttCell * cell,
                                           GfsVariable * v)
{

  GfsSourceParticulate *sp = GFS_SOURCE_PARTICULATE(s);

  FttCellFace f;
  f.cell = cell;

  FttCellNeighbors n;
  ftt_cell_neighbors (cell, &n);

  switch(v->component){
  case FTT_X:
    f.d = FTT_RIGHT;
    f.neighbor = n.c[f.d];
    break;
  case FTT_Y:
    f.d = FTT_TOP;
    f.neighbor = n.c[f.d];
    break;
#if !FTT_2D
  case FTT_Z:
    f.d = FTT_FRONT;
    f.neighbor = n.c[f.d];
    break;
#endif
  default: g_assert_not_reached ();
  }

  if(sp->u[v->component]!=NULL)
    return  gfs_face_interpolated_value_generic(&f, (sp->u[v->component]));
  else
    g_assert_not_reached ();

}

static gdouble source_particulate_centered_value(GfsSourceGeneric * s,
                                                FttCell * cell,
                                                GfsVariable * v)
{

  GfsSourceParticulate *sp = GFS_SOURCE_PARTICULATE(s);

  if(sp->u[v->component]!=NULL)
    return GFS_VALUE(cell,sp->u[v->component]);
  else
    g_assert_not_reached ();

}


static void source_particulate_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_SOURCE_PARTICULATE (o)->kernel_function));

  (* GTS_OBJECT_CLASS (gfs_source_particulate_class ())->parent_class->destroy) (o);
}

void distance_normalization (FttVector * pos1, GfsParticulate * p)
{
  gdouble rb = pow (3.*p->volume/(4.*M_PI), 1./3.);
  FttVector * pos2 = &(GFS_PARTICLE (p)->pos);
  pos1->x = (pos1->x - pos2->x)/rb;
  pos1->y = (pos1->y - pos2->y)/rb;
  pos1->z = 0.;
#if !FTT_2D
  pos1->z = (pos1->z - pos2->z)/rb;
#endif
}

typedef struct {
  gdouble correction,volume;
  GfsParticulate * p;
  GfsFunction * kernel_function;
  GfsVariable **v;
} KernelData;

static void kernel_volume (FttCell * cell, KernelData * kd)
{
  gdouble cellvol = gfs_cell_volume (cell,GFS_DOMAIN(gfs_object_simulation(kd->p)));

  kd->volume += cellvol;

  /* correction term to make a discretely conservative kernel */
  FttVector pos;
  ftt_cell_pos (cell, &pos);
  distance_normalization (&pos, kd->p);
  kd->correction += gfs_function_spatial_value (kd->kernel_function, &pos)*cellvol;
}

typedef struct {
  FttVector * pos;
  gdouble distance;
} CondData;

static gboolean cond_kernel (FttCell * cell, gpointer data)
{

  g_assert(cell!=NULL);

  CondData * p = data;
  FttVector pos;
  ftt_cell_pos (cell, &pos);
  gdouble radeq;
  gdouble size = ftt_cell_size(cell)/2.;

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

static void diffuse_force (FttCell * cell, KernelData * kd)
{
  FttVector pos;
  FttComponent c;
  ftt_cell_pos (cell, &pos);
  distance_normalization (&pos, kd->p);
  gdouble cellvol = gfs_cell_volume (cell,GFS_DOMAIN(gfs_object_simulation(kd->p)));

  GfsSimulation * sim = gfs_object_simulation(kd->p);
  gdouble liq_rho = sim->physical_params.alpha ? 1./
    gfs_function_value (sim->physical_params.alpha, cell) : 1.;

  if (kd->correction > 1.e-10 )
    for (c = 0; c <  FTT_DIMENSION; c++)
      GFS_VALUE (cell, kd->v[c]) -= (&(kd->p->force.x))[c]/liq_rho/cellvol*
        gfs_function_spatial_value (kd->kernel_function, &pos)/kd->correction;

}

static gboolean source_particulate_event (GfsEvent * event,
    GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_source_particulate_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsSourceParticulate * sp = GFS_SOURCE_PARTICULATE (event);
    gdouble influencerad;

    FttComponent c;

    /* Reset variable */
    for (c = 0; c <  FTT_DIMENSION; c++)
      gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
          (FttCellTraverseFunc) gfs_cell_reset, sp->u[c]);

    /* Loop over the list of particles in the selected object */
    /* Compute forces */

    GSList * i = GFS_EVENT_LIST (sp->plist)->list->items;
    GfsParticulate * particulate;

    while (i) {
      particulate = GFS_PARTICULATE(i->data);
      for (c = 0; c < 3; c++)
        (&particulate->force.x)[c] = 0.;
      gts_container_foreach (GTS_CONTAINER (particulate->forces),
    			   (GtsFunc) compute_forces_onfluid, particulate);
      i = i->next;
    }

    i = GFS_EVENT_LIST (sp->plist)->list->items;
    while (i) {
      influencerad = pow (GFS_PARTICULATE(i->data)->volume*3./(4.*M_PI), 1./3.)*sp->rkernel;
      CondData cd = { &GFS_PARTICLE (i->data)->pos, sp->rkernel };
      KernelData kd = { 0., 0., GFS_PARTICULATE(i->data), sp->kernel_function,sp->u};
      gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
          (FttCellTraverseFunc) kernel_volume, &kd,
          cond_kernel, &cd);
      kd.correction /= kd.volume;
      gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
          (FttCellTraverseFunc) diffuse_force, &kd,
          cond_kernel, &cd);

      i = i->next;
    }


    return TRUE;
  }
  return FALSE;
}

static void source_particulate_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_particulate_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (object name)");
    return;
  }
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  GfsSourceParticulate * sp = GFS_SOURCE_PARTICULATE (*o);
  GtsObject * object = gfs_object_from_name (domain, fp->token->str);
  if (object == NULL) {
    gts_file_error (fp, "unknown object '%s'", fp->token->str);
    return;
  }
  if (!GFS_IS_PARTICLE_LIST (object)) {
    gts_file_error (fp, "object '%s' is not a GfsParticleList", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  sp->plist = GFS_PARTICLE_LIST (object);

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
      sp->rkernel = atof (fp->token->str);
      gts_file_next_token (fp);
    }
    else if (!strcmp (fp->token->str, "kernel")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting '='");
	return;
      }
      gts_file_next_token (fp);
      gfs_function_read (sp->kernel_function, gfs_object_simulation (*o), fp);
    }
    else {
      gts_file_error (fp, "unknown keyword `%s'", fp->token->str);
      return;
    }
  }

  /*FixME:
           A. User specifies the Forces which should be evaluated for the sourceterm
           B. Forces are computed for each particle.
           C. Also modify write method
  */
  if (fp->type == GTS_ERROR)
    return;
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return;
  }
  fp->scope_max--;
  gts_file_next_token (fp);

  /* force variable declaration*/
  sp->u[0] = gfs_domain_get_or_add_variable (domain, g_strconcat(GFS_EVENT(sp->plist)->name,"_Fx", NULL), "Fx-component of the source particulate");
  sp->u[1] = gfs_domain_get_or_add_variable (domain, g_strconcat(GFS_EVENT(sp->plist)->name,"_Fy", NULL), "Fy-component of the source particulate");

#if !FTT_2D
  sp->u[2] = gfs_domain_get_or_add_variable (domain, g_strconcat(GFS_EVENT(sp->plist)->name,"_Fz", NULL), "Fz-component of the source particulate");
#endif

}

static void source_particulate_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_particulate_class ())->parent_class->write) (o, fp);

  GfsSourceParticulate * sp = GFS_SOURCE_PARTICULATE (o);

  fprintf (fp, " %s", GFS_EVENT (GFS_SOURCE_PARTICULATE (o)->plist)->name);
  fprintf (fp, " { rkernel = %g ", sp->rkernel);
  fputs (" kernel =", fp);
  gfs_function_write (GFS_SOURCE_PARTICULATE (o)->kernel_function, fp);
  fputc ('}', fp);

}

static void source_particulate_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = source_particulate_event;
  klass->destroy =  source_particulate_destroy;
  klass->read    =  source_particulate_read;
  klass->write   =  source_particulate_write;
}

static void source_particulate_init (GfsSourceParticulate * sp )
{
  sp->kernel_function = gfs_function_new (gfs_function_spatial_class (), 0.);
  GfsSourceGeneric * s = GFS_SOURCE_GENERIC (sp);
  s->mac_value = source_particulate_value;
  s->centered_value = source_particulate_centered_value;
  s->flux = NULL;
  GFS_EVENT (sp)->istep = 1;
}

GfsSourceGenericClass * gfs_source_particulate_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_particulate_info = {
      "GfsSourceParticulate",
      sizeof (GfsSourceParticulate),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_particulate_class_init,
      (GtsObjectInitFunc) source_particulate_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_velocity_class ()),
                                  &gfs_source_particulate_info);
  }

  return klass;
}
/** \endobject{GfsSourceParticulate} */


/** \beginobject{GfsFeedParticle} */

static void feed_particulate (GfsDomain * domain,
			     GfsFeedParticle * feedp)
{
  GfsParticulate newpart;
  newpart.parent.pos.x = gfs_function_value (feedp->posx, NULL);
  newpart.parent.pos.y = gfs_function_value (feedp->posy, NULL);
  newpart.parent.pos.z = gfs_function_value (feedp->posz, NULL);

  GfsVariable ** u = gfs_domain_velocity (domain);
  FttCell * cell = gfs_domain_locate (domain, newpart.parent.pos, -1, NULL);
  if (cell) {
    //newpart.vel.x  = gfs_function_value (feedp->velx, cell);
    newpart.vel.x  = gfs_interpolate (cell, newpart.parent.pos, u[0]);
//    newpart.vel.y  = gfs_function_value (feedp->vely, cell);
    newpart.vel.y  = gfs_interpolate (cell, newpart.parent.pos, u[1]);
//    newpart.vel.z  = gfs_function_value (feedp->velz, cell);
#if FTT_2D
    newpart.vel.z  = 0.;
#else
    newpart.vel.z  = gfs_interpolate (cell, newpart.parent.pos, u[2]);
#endif
    newpart.volume = gfs_function_value (feedp->vol, cell);
    newpart.mass = gfs_function_value (feedp->mass, cell);
    add_particulate (newpart,feedp->plist);
  }
}

static gboolean gfs_feed_particle_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_feed_particle_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsFeedParticle * feedp = GFS_FEED_PARTICLE (event);
    gint i;
    guint np = gfs_function_value (feedp->np, NULL);

    for (i = 0; i < np; i++)
      feed_particulate (domain, feedp);
    return TRUE;
  }
  return FALSE;
}

static void gfs_feed_particle_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->np));
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->posx));
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->posy));
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->posz));
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->velx));
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->vely));
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->velz));
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->mass));
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->vol));

  (* GTS_OBJECT_CLASS (gfs_feed_particle_class ())->parent_class->destroy) (o);
}

static void gfs_feed_particle_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_feed_particle_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_feed_particle_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsFeedParticle * feedp = GFS_FEED_PARTICLE(*o);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (object name)");
    return;
  }

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
  gts_file_next_token (fp);

  feedp->plist = GFS_PARTICLE_LIST (object);


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
  else if (!strcmp (fp->token->str, "nparts")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      gfs_function_read (feedp->np, gfs_object_simulation (*o), fp);
    }
    else if (!strcmp (fp->token->str, "xfeed")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      gfs_function_read (feedp->posx, gfs_object_simulation (*o), fp);
    }
    else if (!strcmp (fp->token->str, "yfeed")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      gfs_function_read (feedp->posy, gfs_object_simulation (*o), fp);
    }
    else if (!strcmp (fp->token->str, "zfeed")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      gfs_function_read (feedp->posz, gfs_object_simulation (*o), fp);
    }
    else if (!strcmp (fp->token->str, "velx")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      gfs_function_read (feedp->velx, gfs_object_simulation (*o), fp);
    }
    else if (!strcmp (fp->token->str, "vely")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      gfs_function_read (feedp->vely, gfs_object_simulation (*o), fp);
    }
    else if (!strcmp (fp->token->str, "velz")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      gfs_function_read (feedp->velz, gfs_object_simulation (*o), fp);
    }
    else if (!strcmp (fp->token->str, "mass")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      gfs_function_read (feedp->mass, gfs_object_simulation (*o), fp);
    }
    else if (!strcmp (fp->token->str, "volume")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      gfs_function_read (feedp->vol, gfs_object_simulation (*o), fp);
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

static void gfs_feed_particle_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_feed_particle_class ())->parent_class->write) (o, fp);

  GfsFeedParticle * feedp = GFS_FEED_PARTICLE(o);

  fprintf (fp, " %s", GFS_EVENT (feedp->plist)->name);

  fputs (" {\n  nparts = ", fp);
  gfs_function_write (feedp->np, fp);
  fputs ("  xfeed =", fp);
  gfs_function_write (feedp->posx, fp);
  fputs (" yfeed =", fp);
  gfs_function_write (feedp->posy, fp);
  fputs (" zfeed =", fp);
  gfs_function_write (feedp->posz, fp);
  fputs ("\n  velx =", fp);
  gfs_function_write (feedp->velx, fp);
  fputs (" vely =", fp);
  gfs_function_write (feedp->vely, fp);
  fputs (" velz =", fp);
  gfs_function_write (feedp->velz, fp);
  fputs ("\n  mass =", fp);
  gfs_function_write (feedp->mass, fp);
  fputs ("\n  volume =", fp);
  gfs_function_write (feedp->vol, fp);
  fputs ("\n}", fp);
}

static void gfs_feed_particle_class_init (GfsEventClass * klass)
{
  GFS_EVENT_CLASS (klass)->event    = gfs_feed_particle_event;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_feed_particle_destroy;
  GTS_OBJECT_CLASS (klass)->read    = gfs_feed_particle_read;
  GTS_OBJECT_CLASS (klass)->write   = gfs_feed_particle_write;
}

static void gfs_feed_particle_init ( GfsFeedParticle * feedp)
{
  feedp->np   = gfs_function_new (gfs_function_class (), 1);
  feedp->posx = gfs_function_new (gfs_function_class (), 0.);
  feedp->posy = gfs_function_new (gfs_function_class (), 0.);
  feedp->posz = gfs_function_new (gfs_function_class (), 0.);
  feedp->velx = gfs_function_new (gfs_function_class (), 0.);
  feedp->vely = gfs_function_new (gfs_function_class (), 0.);
  feedp->velz = gfs_function_new (gfs_function_class (), 0.);
  feedp->mass = gfs_function_new (gfs_function_class (), 0.);
  feedp->vol  = gfs_function_new (gfs_function_class (), 0.);
}

GfsEventClass * gfs_feed_particle_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_feed_particle_info = {
      "GfsFeedParticle",
      sizeof (GfsFeedParticle),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_feed_particle_class_init,
      (GtsObjectInitFunc) gfs_feed_particle_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_class ()),
				  &gfs_feed_particle_info);
  }
  return klass;
}

/** \endobject{GfsFeedParticle} */

/** \beginobject{GfsOutputParticleList} */

static void fill_output_location(GfsEvent * event, GfsOutputLocation * l)
{
    GfsParticle *p = GFS_PARTICLE(event);
    FttVector pos;
    pos = p->pos;
    pos.r = pow(3.*GFS_PARTICULATE(event)->volume/(4.*M_PI),1./3.);
    g_array_append_val (l->p, pos);
}

static void gfs_output_particle_list_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

    GfsOutputParticleList * p = GFS_OUTPUT_PARTICLE_LIST (*o);
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
    gts_file_next_token (fp);

    p->plist = GFS_PARTICLE_LIST (object);

    GfsOutputLocation * l = GFS_OUTPUT_LOCATION (p);
    GfsEventList *el = GFS_EVENT_LIST(p->plist);
    gts_container_foreach (GTS_CONTAINER (el->list), (GtsFunc)fill_output_location, l);

}

static void gfs_output_particle_list_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s",  GFS_EVENT (GFS_OUTPUT_PARTICLE_LIST (o)->plist)->name);
}


static gboolean gfs_output_particle_list_event (GfsEvent * event,
					 GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_class ())->parent_class)->event)
      (event, sim)) {
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_particle_list_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event  = gfs_output_particle_list_event;
  GTS_OBJECT_CLASS (klass)->read  = gfs_output_particle_list_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_output_particle_list_write;
}

GfsOutputClass * gfs_output_particle_list_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_particle_list_info = {
      "GfsOutputParticleList",
      sizeof (GfsOutputParticleList),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_particle_list_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_particle_class ()),
				  &gfs_output_particle_list_info);
  }

  return klass;
}

/** \endobject{GfsOutputParticleList} */

/**
 * Source terms for particle variables.
 * \beginobject{GfsSourceParticulateVol}
 */

static void source_particulatevol_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_particulatevol_class ())->parent_class->write) (o, fp);

  GfsSourceParticulateVol * pv = GFS_SOURCE_PARTICULATEVOL (o);

  fprintf (fp, " %s", GFS_EVENT (pv->plist)->name);
  gfs_function_write (pv->source, fp);

  if (pv->volsource)
    fprintf (fp, " %s", pv->volsource->name);
}

static void source_particulatevol_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceParticulateVol * pv;
  GfsDomain * domain;

  if (GTS_OBJECT_CLASS (gfs_source_particulatevol_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_particulatevol_class ())->parent_class->read)
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  pv = GFS_SOURCE_PARTICULATEVOL (*o);
  domain =  GFS_DOMAIN (gfs_object_simulation (pv));

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (object name)");
    return;
  }

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

  pv->plist = GFS_PARTICLE_LIST (object);
  gts_file_next_token (fp);

  pv->rad = gfs_domain_get_or_add_variable (domain, "Rad",
						  "Particle Radius");
  pv->u_rel = gfs_domain_get_or_add_variable (domain, "Urelp",
      "Particle x - relative velocity");
  pv->v_rel = gfs_domain_get_or_add_variable (domain, "Vrelp", "Particle y - relative velocity");
#if !FTT_2D
  pv->w_rel = gfs_domain_get_or_add_variable (domain, "Wrelp",
      "Particle z - relative velocity");
#endif

  gfs_function_read (pv->source, gfs_object_simulation (*o), fp);

  if (fp->type == GTS_STRING) {
    if (!(pv->volsource = gfs_variable_from_name (domain->variables, fp->token->str))) {
      gts_file_error (fp, "unknown variable `%s'", fp->token->str);
      return;
    }
    gts_file_next_token (fp);
  }
}

static void update_vol (GfsEvent * event, GfsSourceParticulateVol * spv)
{
  GfsParticulate *p = GFS_PARTICULATE(event);
  GfsSimulation *sim = gfs_object_simulation(p);
  GfsDomain * domain = GFS_DOMAIN(sim);
  FttCell * cell = gfs_domain_locate (domain, GFS_PARTICLE(p)->pos, -1, NULL);
  GFS_VALUE (cell, spv->rad) = pow(3.0*p->volume/4.0/M_PI,1./3.);

  GfsVariable ** u = gfs_domain_velocity (domain);
  FttVector fluid_vel;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&fluid_vel.x)[c] = gfs_interpolate (cell, GFS_PARTICLE(p)->pos, u[c]);
  FttVector relative_vel = subs_fttvectors (&fluid_vel, &p->vel);
  GFS_VALUE (cell, spv->u_rel) = relative_vel.x;
  GFS_VALUE (cell, spv->v_rel) = relative_vel.y;
#if !FTT_2D
  GFS_VALUE (cell, spv->w_rel) = relative_vel.z;
#endif

  gdouble source = gfs_function_value (spv->source, cell);

  p->volume += source*sim->advection_params.dt;
  if (spv->volsource)
    GFS_VALUE(cell,spv->volsource) += source;

}

static gboolean source_particulatevol_event (GfsEvent * event,
					 GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_source_particulatevol_class ())->parent_class)->event)
      (event, sim)) {

    GfsSourceParticulateVol * spv = GFS_SOURCE_PARTICULATEVOL(event);
    GfsEventList *l = GFS_EVENT_LIST(spv->plist);

    if (spv->volsource)
      gfs_domain_cell_traverse (GFS_DOMAIN(sim), FTT_PRE_ORDER, FTT_TRAVERSE_ALL, 1,
                               (FttCellTraverseFunc) gfs_cell_reset, spv->volsource);

    gts_container_foreach (GTS_CONTAINER (l->list), (GtsFunc)update_vol, spv);

    return TRUE;
  }
  return FALSE;
}

static void source_particulatevol_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read =  source_particulatevol_read;
  GTS_OBJECT_CLASS (klass)->write = source_particulatevol_write;
  GFS_EVENT_CLASS (klass)->event = source_particulatevol_event;
}

static void source_particulatevol_init (GfsSourceGeneric * s)
{
    GFS_SOURCE_PARTICULATEVOL(s)->source = gfs_function_new (gfs_function_class (), 0.);
}

GfsSourceGenericClass * gfs_source_particulatevol_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_particulatevol_info = {
      "GfsSourceParticulateVol",
      sizeof (GfsSourceParticulateVol),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_particulatevol_class_init,
      (GtsObjectInitFunc) source_particulatevol_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_generic_class ()),
				  &gfs_source_particulatevol_info);
  }

  return klass;
}

/** \endobject{GfsSourceParticulateVol} */

/**
 * Source terms for particle variables.
 * \beginobject{GfsSourceParticulateMass}
 */

static void source_particulatemass_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_source_particulatemass_class ())->parent_class->write) (o, fp);

  GfsSourceParticulateMass * pv = GFS_SOURCE_PARTICULATEMASS (o);

  fprintf (fp, " %s", GFS_EVENT (pv->plist)->name);
  gfs_function_write (pv->source, fp);

  if (pv->dmdtvariable)
    fprintf (fp, " %s", pv->dmdtvariable->name);

}

static void source_particulatemass_read (GtsObject ** o, GtsFile * fp)
{
  GfsSourceParticulateMass * pv;
  GfsDomain * domain;

  if (GTS_OBJECT_CLASS (gfs_source_particulatemass_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_source_particulatemass_class ())->parent_class->read)
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  pv = GFS_SOURCE_PARTICULATEMASS (*o);
  domain =  GFS_DOMAIN (gfs_object_simulation (pv));

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (object name)");
    return;
  }

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

  pv->plist = GFS_PARTICLE_LIST (object);
  gts_file_next_token (fp);

  pv->rad = gfs_domain_get_or_add_variable (domain, "Rad",
      "Particle Radius");
  pv->u_rel = gfs_domain_get_or_add_variable (domain, "Urelp",
      "Particle x - relative velocity");
  pv->v_rel = gfs_domain_get_or_add_variable (domain, "Vrelp",
      "Particle y - relative velocity");
#if !FTT_2D
  pv->w_rel = gfs_domain_get_or_add_variable (domain, "Wrelp",
      "Particle z - relative velocity");
#endif

  gfs_function_read (pv->source, gfs_object_simulation (*o), fp);

  if (fp->type == GTS_STRING) {
    if (!(pv->dmdtvariable = gfs_domain_get_or_add_variable (domain, fp->token->str, "mass source"))) {
      gts_file_error (fp, "unknown variable `%s'", fp->token->str);
      return;
    }
    gts_file_next_token (fp);
  }

}

static void update_mass (GfsEvent * event, GfsSourceParticulateMass * spv)
{
  GfsParticulate *p = GFS_PARTICULATE(event);
  GfsSimulation *sim = gfs_object_simulation(p);
  GfsDomain * domain = GFS_DOMAIN(sim);
  FttCell * cell = gfs_domain_locate (domain, GFS_PARTICLE(p)->pos, -1, NULL);
  GFS_VALUE (cell, spv->rad) = pow(3.0*p->volume/4.0/M_PI,1./3.);

  GfsVariable ** u = gfs_domain_velocity (domain);
  FttVector fluid_vel;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&fluid_vel.x)[c] = gfs_interpolate (cell, GFS_PARTICLE(p)->pos, u[c]);
  FttVector relative_vel = subs_fttvectors (&fluid_vel, &p->vel);
  GFS_VALUE (cell, spv->u_rel) = relative_vel.x;
  GFS_VALUE (cell, spv->v_rel) = relative_vel.y;
#if !FTT_2D
  GFS_VALUE (cell, spv->w_rel) = relative_vel.z;
#endif

  gdouble source = gfs_function_value (spv->source, cell);

  p->mass += source*sim->advection_params.dt;

  if (spv->dmdtvariable)
    GFS_VALUE(cell, spv->dmdtvariable) += source;

}

static gboolean source_particulatemass_event (GfsEvent * event,
					 GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_source_particulatemass_class ())->parent_class)->event)
      (event, sim)) {

    GfsDomain * domain = GFS_DOMAIN(sim);
    GfsSourceParticulateMass * spv = GFS_SOURCE_PARTICULATEMASS(event);
    GfsEventList *l = GFS_EVENT_LIST(spv->plist);

    if (spv->dmdtvariable)
        gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, 1,
                           (FttCellTraverseFunc) gfs_cell_reset, spv->dmdtvariable);

    gts_container_foreach (GTS_CONTAINER (l->list), (GtsFunc)update_mass, spv);

    return TRUE;
  }
  return FALSE;
}

static void source_particulatemass_class_init (GfsSourceGenericClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read =  source_particulatemass_read;
  GTS_OBJECT_CLASS (klass)->write = source_particulatemass_write;
  GFS_EVENT_CLASS (klass)->event = source_particulatemass_event;
}

static void source_particulatemass_init (GfsSourceGeneric * s)
{
    GFS_SOURCE_PARTICULATEMASS(s)->source = gfs_function_new (gfs_function_class (), 0.);
}

GfsSourceGenericClass * gfs_source_particulatemass_class (void)
{
  static GfsSourceGenericClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_source_particulatemass_info = {
      "GfsSourceParticulateMass",
      sizeof (GfsSourceParticulateMass),
      sizeof (GfsSourceGenericClass),
      (GtsObjectClassInitFunc) source_particulatemass_class_init,
      (GtsObjectInitFunc) source_particulatemass_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_source_generic_class ()),
				  &gfs_source_particulatemass_info);
  }

  return klass;
}

/** \endobject{GfsSourceParticulateMass} */

/** /beginobject{GfsBcParticle} */

typedef struct {
  gint boxid;
  FttDirection d;
  GfsParticle *p;
} Particle_send;

/*Checks the intersection of a particle path ray(a line segment) with the boundaries of the containee cell*/
static gboolean check_intersetion(FttVector cellpos, FttVector p0, FttVector p1,
				  FttDirection *dstore, gdouble size)
{

  gdouble t;
  FttDirection d;


  for(d = 0; d < FTT_NEIGHBORS; d++){

    gdouble normal = ((gdouble) FTT_OPPOSITE_DIRECTION(d) - (gdouble)d);

#if FTT_2D
    switch(d/2){
    case 0:
      if((p1.x - p0.x)!=0 && normal*(p1.x-p0.x) > 0){
	t = (cellpos.x + normal*size*0.5 - p0.x)/(p1.x - p0.x);
	gdouble py = p0.y + t*(p1.y - p0.y);
	if((py - cellpos.y + size*0.5)*(py - cellpos.y - size*0.5) <= 0  &&  t*(t-1)<= 0 ){
	  *dstore = d;
	  return TRUE;
	}
      }
      break;
    case 1:
      if((p1.y - p0.y)!=0 && normal*(p1.y-p0.y) > 0){
	t = (cellpos.y + normal*size*0.5- p0.y)/(p1.y - p0.y);

	gdouble px = p0.x + t*(p1.x - p0.x);
	if((px - cellpos.x + size*0.5)*(px - cellpos.x - size*0.5) <= 0  &&  t*(t-1) <= 0){
	  *dstore = d;
	  return TRUE;
	}
      }
      break;
    default: g_assert_not_reached();
    }
#else
    switch(d/2){
    case 0:
      if((p1.x - p0.x)!=0 && normal*(p1.x-p0.x) > 0){
	t = (cellpos.x + normal*size*0.5 - p0.x)/(p1.x - p0.x);
	gdouble py = p0.y + t*(p1.y - p0.y);
	gdouble pz = p0.z + t*(p1.z - p0.z);
	if((py - cellpos.y + size*0.5)*(py - cellpos.y - size*0.5) <= 0  &&
	   (pz - cellpos.z + size*0.5)*(pz - cellpos.z - size*0.5) <= 0
	   &&  t*(t-1)<= 0 )
	  {
	  *dstore = d;
	  return TRUE;
	}
      }
      break;
    case 1:
      if((p1.y - p0.y)!=0 && normal*(p1.y-p0.y) > 0){
	t = (cellpos.y + normal*size*0.5- p0.y)/(p1.y - p0.y);

	gdouble px = p0.x + t*(p1.x - p0.x);
	gdouble pz = p0.z + t*(p1.z - p0.z);
	if((px - cellpos.x + size*0.5)*(px - cellpos.x - size*0.5) <= 0  &&
	   (pz - cellpos.z + size*0.5)*(pz - cellpos.z - size*0.5) <= 0
	   &&  t*(t-1)<= 0 )
	  {
	  *dstore = d;
	  return TRUE;
	}
      }
      break;
    case 2:
      if((p1.z - p0.z)!=0 && normal*(p1.z-p0.z) > 0){
	t = (cellpos.z + normal*size*0.5- p0.z)/(p1.z - p0.z);

	gdouble px = p0.x + t*(p1.x - p0.x);
	gdouble py = p0.y + t*(p1.y - p0.y);
	if((px - cellpos.x + size*0.5)*(px - cellpos.x - size*0.5) <= 0  &&
	   (py - cellpos.y + size*0.5)*(py - cellpos.y - size*0.5) <= 0
	   &&  t*(t-1)<= 0 )
	  {
	  *dstore = d;
	  return TRUE;
	}
      }
      break;
    default: g_assert_not_reached();
    }
#endif
  }

  g_warning("Intersection finding algorithm failed\n");
  return FALSE;
}

/*Tracks the particle path ray to identify boundary cell for the application of the Boundary Conditions*/
static FttCell * boundarycell (GfsDomain *domain, GfsParticle *p, FttDirection *dstore)
{

  FttCell * cell = gfs_domain_locate(domain, p->pos_old, -1, NULL);
  g_assert(cell!=NULL);

  FttVector cellpos;
  gdouble size;
  ftt_cell_pos(cell, &cellpos);
  size = ftt_cell_size(cell);

  check_intersetion(cellpos, p->pos_old, p->pos, dstore, size);
  FttCellFace face = ftt_cell_face(cell, *dstore);

  if(!face.neighbor)
    return cell;

  while(!GFS_CELL_IS_BOUNDARY(face.neighbor)) {
    cell = face.neighbor;
    g_assert(cell!=NULL);
    ftt_cell_pos(cell, &cellpos);
    size = ftt_cell_size(cell);

    check_intersetion(cellpos, p->pos_old, p->pos, dstore, size);
    face = ftt_cell_face(cell, *dstore);

    if(!face.neighbor)
      return cell;
  };

  if(cell)
    return cell;
  else
    g_assert_not_reached();

}

/*Periodic boundary conditions*/
static void periodic_bc_particle(FttDirection d, GfsBox * box, GSList *particles, GfsBcParticle *bc)
{
  FttVector box_face, box_face_nbr;
  GfsBox *box_nbr = GFS_BOX(GFS_BOUNDARY_PERIODIC(box->neighbor[d])->matching);
  g_assert(box->root!=NULL);
  g_assert(box_nbr->root!=NULL);
  ftt_cell_pos(box->root, &box_face);
  ftt_cell_pos(box_nbr->root, &box_face_nbr);
  gdouble size = ftt_cell_size(box->root);

  gdouble normal = (gdouble)FTT_OPPOSITE_DIRECTION(d) - (gdouble) d;
  (&box_face.x)[d/2] += (gdouble)normal * size/2.;
  (&box_face_nbr.x)[d/2] -= (gdouble)normal * size/2.;

  GSList *i = particles;
  gdouble tolerance = size/1.e8;
  while(i){
      GfsParticle *p = (GfsParticle *) (i->data);
      i = i->next;
      gdouble distance = ((&p->pos.x)[d/2] - (&box_face.x)[d/2])*normal;
      (&p->pos.x)[d/2] = (&box_face_nbr.x)[d/2] + distance + normal*tolerance;
      (&p->pos_old.x)[d/2] = (&p->pos.x)[d/2];
      particles = g_slist_remove(particles,p);
      gts_container_add(GTS_CONTAINER(GFS_EVENT_LIST(bc->plist)->list),GTS_CONTAINEE(p));
  }
}

#ifdef HAVE_MPI
/*MPI Boundary conditions*/
static void mpi_send_particle(gint dest, GSList *particles)
{
  GfsRequest * request = gfs_send_objects(particles, dest);
  gfs_wait(request);
}

static void mpi_rcv_particle(gint src, GfsParticleList *plist)
{
  GfsSimulation *sim = gfs_object_simulation(plist);
  GfsDomain *domain = GFS_DOMAIN(sim);
  GfsEventList *l = GFS_EVENT_LIST(plist);

  GSList *list = gfs_receive_objects (domain, src);

  while(list){
    GtsObject * object = (GtsObject *)(list->data);
    GfsEvent * event = GFS_EVENT (l);
    gfs_event_set (GFS_EVENT (object),
		   event->start, event->end, event->step, event->istart, event->iend, event->istep);
    gts_container_add(GTS_CONTAINER(l->list),GTS_CONTAINEE(object));
    if(GFS_IS_PARTICULATE(object))
      GFS_PARTICULATE(object)->forces = plist->forces;
    list = list->next;
  }

  g_slist_free(list);
}

#endif /*HAVE_MPI*/

static void send_particles(FttDirection d, GfsBoundary *b, GfsBox *box, gint nsends,
			   GSList *particles, GfsBcParticle *bc)
{
#ifdef HAVE_MPI
  GfsBoundaryMpi *mpi = GFS_BOUNDARY_MPI (b);
  if(GFS_IS_BOUNDARY_MPI(mpi)){
    /*MPI BC*/
    mpi_send_particle(mpi->process, particles);
    return;
  }
#endif /*HAVE_MPI*/
  if(nsends > 0){
    if(box->neighbor[d] && GFS_IS_BOUNDARY_PERIODIC(box->neighbor[d])){
      /*Periodic BC*/
      g_assert(GFS_IS_BOX(GFS_BOUNDARY(box->neighbor[d])->box));
      periodic_bc_particle(d, box, particles, bc);
    }
  }
}

static void rcv_particles(GfsBoundary *b, GfsParticleList *plist )
{
#ifdef HAVE_MPI
  GfsBoundaryMpi *mpi = GFS_BOUNDARY_MPI (b);
  if(GFS_IS_BOUNDARY_MPI(mpi)){
    /*MPI BC*/
    mpi_rcv_particle(mpi->process, plist);
  }
#endif /*HAVE_MPI*/
}

static void box_send_bc(GfsBox *box, GfsBcParticle *bc)
{

  gint nsends[FTT_NEIGHBORS];
  GSList *packet_send[FTT_NEIGHBORS];
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++){
    nsends[d] = 0;
    packet_send[d] = NULL;
  }


  GSList * i = bc->boundary;
  while(i){
    Particle_send *psend = (Particle_send *) (i->data);
    i = i->next;
    if(psend->boxid == box->id){
      nsends[psend->d]++;
      packet_send[psend->d] = g_slist_append(packet_send[psend->d], psend->p);
      bc->boundary = g_slist_remove(bc->boundary,psend);
      g_free(psend);
    }
  }

  for (d = 0; d < FTT_NEIGHBORS; d++){
    /*Fixme: If neighboring box does not exist -> slip-wall (Reflect)*/
    if(GFS_IS_BOUNDARY(box->neighbor[d])){
      GfsBoundary *b = GFS_BOUNDARY(box->neighbor[d]);
      send_particles(d, b, box, nsends[d], packet_send[d], bc);
    }
  }

}

static void box_rcv_bc(GfsBox *box, GfsBcParticle *bc)
{
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++){
    if(box->neighbor[d] && GFS_IS_BOUNDARY(box->neighbor[d])){
      GfsBoundary *b = GFS_BOUNDARY(box->neighbor[d]);
      rcv_particles(b, bc->plist);
    }
  }
}

static void list_boundary_particles(GfsEvent * event, GfsBcParticle *bc)
{

  GfsSimulation *sim = gfs_object_simulation(bc->plist);
  GfsDomain *domain = GFS_DOMAIN(sim);

  GfsParticle *p = GFS_PARTICLE(event);
  FttCell * cellp = gfs_domain_locate(domain, p->pos, -1, NULL);
  if(cellp==NULL){
    gts_container_remove(GTS_CONTAINER(GFS_EVENT_LIST(bc->plist)->list),GTS_CONTAINEE(event));
    FttDirection d;
    FttCell *cell = boundarycell (domain, p, &d);
    if(cell){
      while(cell && !FTT_CELL_IS_ROOT(cell))
	cell = ftt_cell_parent(cell);

      GfsBox *box = GFS_BOX(FTT_ROOT_CELL(cell)->parent);

      if(d < 0 || d > FTT_NEIGHBORS)
	g_assert_not_reached();

      if(!box->neighbor[d] || GFS_IS_BOUNDARY(box->neighbor[d])){
      	Particle_send *psend = g_malloc(sizeof(Particle_send));
	psend->boxid = box->id;
	psend->d = d;
	psend->p = p;
	bc->boundary = g_slist_append(bc->boundary, psend);
      }
    }
    else
      gts_object_destroy(GTS_OBJECT(event));
  }

}

static void free_boundary_particles(GfsBcParticle *bc)
{
  GSList *i = bc->boundary;

  while(i){
    Particle_send *psend = (Particle_send *) (i->data);
    i = i->next;
    bc->boundary = g_slist_remove(bc->boundary,psend);
    g_free(psend);
  }
  g_free(bc->boundary);
}


void gfs_particle_bc (GfsParticleList *plist)
{
    GfsEventList *l = GFS_EVENT_LIST(plist);
    GfsBcParticle *bc = g_malloc0 (sizeof(GfsBcParticle));

    GfsSimulation *sim = gfs_object_simulation(plist);
    GfsDomain *domain = GFS_DOMAIN(sim);

    bc->plist = plist;
    bc->boundary = NULL;
    gts_container_foreach (GTS_CONTAINER (l->list), (GtsFunc)list_boundary_particles, bc);

    gts_container_foreach (GTS_CONTAINER (domain),
    			   (GtsFunc) box_send_bc, bc);

    gts_container_foreach (GTS_CONTAINER (domain),
    			   (GtsFunc) box_rcv_bc, bc);

    free_boundary_particles(bc);
    g_free(bc);
}

/** \end gfs_particle_bc */
