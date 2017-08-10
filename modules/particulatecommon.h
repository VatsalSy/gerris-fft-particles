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

#include "particle.h"
#include "refine.h"
#include "adaptive.h"
#include "domain.h"
#include "fluid.h"
#include "solid.h"
#include "vof.h"
#include "source.h"
#include "output.h"
#include "mpi_boundary.h"

/* GfsParticulate: Header */

typedef struct _GfsParticulate GfsParticulate;

struct _GfsParticulate {
  GfsParticle parent;
  FttVector vel;
  gdouble mass, volume;
  FttVector force;
  GtsSListContainer * forces;
};

#define GFS_PARTICULATE(obj)            GTS_OBJECT_CAST (obj,		\
							 GfsParticulate, gfs_particulate_class ())
#define GFS_IS_PARTICULATE(obj)         (gts_object_is_from_class (obj, gfs_particulate_class ()))

GfsEventClass * gfs_particulate_class  (void);

/* GfsParticleList: Header */

typedef struct _GfsParticleList GfsParticleList;

struct _GfsParticleList {
  /*< private >*/
  GfsEventList parent;

  /*< public >*/
  gint idlast, maxid;
  GtsSListContainer * forces;
  gboolean first_call;
};

#define GFS_PARTICLE_LIST(obj)            GTS_OBJECT_CAST (obj,		\
							   GfsParticleList, \
							   gfs_particle_list_class ())

#define GFS_IS_PARTICLE_LIST(obj)         (gts_object_is_from_class (obj, \
								     gfs_particle_list_class ()))

GfsEventClass * gfs_particle_list_class  (void);

/* GfsParticleForce: header */

typedef struct _GfsParticleForce GfsParticleForce;

struct _GfsParticleForce{
  GtsSListContainee parent;
  FttVector (* force) (GfsParticle *p, GfsParticleForce *force);
};

#define GFS_PARTICLE_FORCE(obj)            GTS_OBJECT_CAST (obj,		\
							GfsParticleForce, \
							gfs_particle_force_class ())
#define GFS_IS_PARTICLE_FORCE(obj)         (gts_object_is_from_class (obj, \
								      gfs_particle_force_class ()))

GtsSListContaineeClass * gfs_particle_force_class  (void);

/* GfsForceCoeff: header */

typedef struct _GfsForceCoeff GfsForceCoeff;

struct _GfsForceCoeff{
  /*< private >*/
  GfsParticleForce parent;

  /*< public >*/
  GfsFunction * coefficient;
  GfsVariable *re_p, *u_rel, *v_rel, *w_rel, *pdia;
  GfsVariable *Uold[3];
  GfsParticulate *p;
};

#define FORCE_COEFF(obj)            GTS_OBJECT_CAST (obj,		\
						    GfsForceCoeff,		\
						    gfs_force_coeff_class ())
#define GFS_IS_FORCE_COEFF(obj)         (gts_object_is_from_class (obj,	\
								  gfs_force_coeff_class ()))
GtsSListContaineeClass * gfs_force_coeff_class  (void);

/* GfsForceAddedMass: header */

#define GFS_IS_FORCE_ADDEDMASS(obj)         (gts_object_is_from_class (obj,	\
								  gfs_force_addedmass_class ()))
GtsSListContaineeClass * gfs_force_addedmass_class  (void);

/* GfsForceInertial: header */

#define GFS_IS_FORCE_INERTIAL(obj)         (gts_object_is_from_class (obj,	\
								  gfs_force_inertial_class ()))
GtsSListContaineeClass * gfs_force_inertial_class  (void);

/* GfsForceLift: header */

#define GFS_IS_FORCE_LIFT(obj)         (gts_object_is_from_class (obj,	\
								  gfs_force_lift_class ()))
GtsSListContaineeClass * gfs_force_lift_class  (void);

/* GfsForceDrag: header */

#define GFS_IS_FORCE_DRAG(obj)         (gts_object_is_from_class (obj,	\
								  gfs_force_drag_class ()))
GtsSListContaineeClass * gfs_force_drag_class  (void);

/* GfsForceBuoy: header */

#define GFS_IS_FORCE_BUOY(obj)         (gts_object_is_from_class (obj,	\
								  gfs_force_buoy_class ()))
GtsSListContaineeClass * gfs_force_buoy_class  (void);

/* GfsDropletToParticle: header */

typedef struct _GfsDropletToParticle GfsDropletToParticle;

struct _GfsDropletToParticle{
  /*< private >*/
  GfsEvent parent;

  /*< public >*/
  GfsParticleList * plist;
  GfsFunction * fc;
  GfsVariable * c;
  gint min;
  gdouble resetwith;
  gdouble density;
};

#define DROPLET_TO_PARTICLE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsDropletToParticle,\
					         gfs_droplet_to_particle_class ())

#define IS_DROPLET_TO_PARTICLE(obj)         (gts_object_is_from_class (obj,\
						 gfs_droplet_to_particle_class ()))

GfsEventClass * gfs_droplet_to_particle_class  (void);

/* GfsParticleToDroplet: header */

typedef struct _GfsParticleToDroplet                GfsParticleToDroplet;

struct _GfsParticleToDroplet {
  /*< private >*/
  GfsEvent parent;

  /*< public >*/
  GfsParticleList * plist;
  GfsFunction * fc;
  GfsSurface *shape;
  gint maxlevel;
  gdouble resetwith;
  GfsVariable * c;
};


#define GFS_PARTICLE_TO_DROPLET(obj)            GTS_OBJECT_CAST (obj,\
                                                   GfsParticleToDroplet,\
                                                   gfs_particle_to_droplet_class ())
#define GFS_IS_PARTICLE_TO_DROPLET(obj)         (gts_object_is_from_class (obj,\
                                                   gfs_particle_to_droplet_class ()))

GfsEventClass * gfs_particle_to_droplet_class  (void);


/* GfsParticulateField: header */

typedef struct _GfsParticulateField                GfsParticulateField;

struct _GfsParticulateField {
  /*< private >*/
  GfsVariable parent;

  /*< public >*/
  GfsParticleList * plist;
  void (* voidfraction_func) (FttCell *,
			      GfsVariable *,
			      GfsParticulate *);
};


#define GFS_PARTICULATE_FIELD(obj)            GTS_OBJECT_CAST (obj,\
                                                   GfsParticulateField,\
                                                   gfs_particulate_field_class ())
#define GFS_IS_PARTICULATE_FIELD(obj)         (gts_object_is_from_class (obj,\
                                                   gfs_particulate_field_class ()))

GfsVariableClass * gfs_particulate_field_class  (void);

/* GfsSourceLangrangian: header */

typedef struct _GfsSourceParticulate                GfsSourceParticulate;

struct _GfsSourceParticulate {
  /*< private >*/
  GfsSourceVelocity parent;
  GfsVariable * u[FTT_DIMENSION];

  /*< public >*/
  GfsParticleList * plist;
  gdouble rkernel;
  GfsFunction * kernel_function;
};


#define GFS_SOURCE_PARTICULATE(obj)            GTS_OBJECT_CAST (obj,\
                                                   GfsSourceParticulate,\
                                                   gfs_source_Particulate_class ())
#define GFS_IS_SOURCE_PARTICULATE(obj)         (gts_object_is_from_class (obj,\
                                                   gfs_source_Particulate_class ()))

GfsSourceGenericClass * gfs_source_particulate_class  (void);

/* GfsFeedParticle: header */

typedef struct _GfsFeedParticle GfsFeedParticle;

struct _GfsFeedParticle{
  /*< private >*/
  GfsEvent parent;

  /*< public >*/
  GfsParticleList * plist;
  GfsFunction * posx, * posy, * posz;
  GfsFunction * velx, * vely, * velz;
  GfsFunction * np,   * mass, * vol;
};

#define GFS_FEED_PARTICLE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsFeedParticle,\
					         gfs_feed_particle_class ())

#define GFS_IS_FEED_PARTICLE(obj)         (gts_object_is_from_class (obj,\
					         gfs_feed_particle_class ()))

GfsEventClass * gfs_feed_particle_class  (void);

/* GfsOutputParticleList: Header */

typedef struct _GfsOutputParticleList GfsOutputParticleList;

struct _GfsOutputParticleList {
  /*< private >*/
  GfsOutputLocation parent;

  /*< public >*/
  GfsParticleList * plist;
};

#define GFS_OUTPUT_PARTICLE_LIST(obj)            GTS_OBJECT_CAST (obj,\
					         GfsOutputParticleList,\
					         gfs_output_particle_list_class ())

#define GFS_IS_OUTPUT_PARTICLE_LIST(obj)     (gts_object_is_from_class (obj,\
								   gfs_output_particle_list_class ()))

GfsOutputClass * gfs_output_particle_list_class  (void);

/* GfsSourceParticulateVol: Header */

typedef struct _GfsSourceParticulateVol         GfsSourceParticulateVol;

struct _GfsSourceParticulateVol {
  /*< private >*/
  GfsSourceGeneric parent;
  GfsVariable * rad, *volsource, * u_rel, * v_rel, * w_rel;

  /*< public >*/
  GfsParticleList * plist;
  GfsFunction * source;
};

#define GFS_SOURCE_PARTICULATEVOL(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceParticulateVol,\
					         gfs_source_particulatevol_class ())
#define GFS_IS_SOURCE_PARTICULATEVOL(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_particulatevol_class ()))

GfsSourceGenericClass * gfs_source_particulatevol_class  (void);

/* GfsSourceParticulateMass: Header */

typedef struct _GfsSourceParticulateMass         GfsSourceParticulateMass;

struct _GfsSourceParticulateMass {
  /*< private >*/
  GfsSourceGeneric parent;
  GfsVariable * rad, * u_rel, * v_rel, * w_rel;

  /*< public >*/
  GfsParticleList * plist;
  GfsFunction * source;
  GfsVariable * dmdtvariable;
};

#define GFS_SOURCE_PARTICULATEMASS(obj)            GTS_OBJECT_CAST (obj,\
					         GfsSourceParticulateMass,\
					         gfs_source_particulatemass_class ())
#define GFS_IS_SOURCE_PARTICULATEMASS(obj)         (gts_object_is_from_class (obj,\
						 gfs_source_particulatemass_class ()))

GfsSourceGenericClass * gfs_source_particulatemass_class  (void);


/* gfs_particle_bc: header */
typedef struct _GfsBcParticle GfsBcParticle;
struct _GfsBcParticle{
  GfsParticleList * plist;
  GSList * boundary;
};

void gfs_particle_bc (GfsParticleList *plist);
void distance_normalization (FttVector * pos1, GfsParticulate * p);
gint particle_id (GfsParticleList * plist);
void assign_forces (GfsParticulate *particulate, GtsSListContainer *forces);
void compute_forces (GfsParticleForce * event, GfsParticulate * p);
