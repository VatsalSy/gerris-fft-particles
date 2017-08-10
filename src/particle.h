/* Gerris - The GNU Flow Solver
 * Copyright (C) 2009 National Institute of Water and Atmospheric Research
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
 * Edited by Vatsal Sanjay 9/08 @ 22:50 hours
 */

#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#ifdef __cplusplus // if this is C++ compiler
extern "C" { // force C compiler
#endif /* __cplusplus */

#include "event.h"

/* Particle: header -- definitions */

typedef struct _GfsParticle GfsParticle;

struct _GfsParticle {
  GfsEvent parent;
  FttVector pos;
  FttVector pos_old;
  guint id;
};

#define GFS_PARTICLE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsParticle,\
					         gfs_particle_class ())
#define GFS_IS_PARTICLE(obj)         (gts_object_is_from_class (obj,\
						 gfs_particle_class ()))

GfsEventClass * gfs_particle_class  (void);

#ifdef __cplusplus /* if this is C++ compiler */
}
#endif /* __cplusplus */

#endif /* __PARTICLE_H__ */
