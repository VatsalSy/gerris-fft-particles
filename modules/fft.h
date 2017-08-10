/* Gerris - The GNU Flow Solver			(-*-C-*-)
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
 */

#include "output.h"
#include "cartesian.h"
#include "vof.h"
#include "tension.h"
#include "mpi_boundary.h"
//#include <fftw3.h>
#include <fftw3-mpi.h>

/* GfsCommData: header */

typedef struct _GfsCommData                     GfsCommData;

typedef struct {
  gint node, src, np;
} CommListData;

struct _GfsCommData {
  /*< private >*/
  GtsObject parent;

  /*< private >*/
  GSList * data;
};

#define GFS_COMM_DATA(obj)                 GTS_OBJECT_CAST (obj,\
                                                            GfsCommData, \
                                                            gfs_comm_data_class ())
#define GFS_IS_COMM_DATA(obj)             (gts_object_is_from_class (obj,\
                                                                     gfs_comm_data_class ()))

GtsObjectClass * gfs_comm_data_class  (void);

GfsCommData * comm_data_new (GtsObjectClass * klass);

/* GfsOutputSpectra: header */

typedef struct _GfsOutputSpectra                     GfsOutputSpectra;

struct _GfsOutputSpectra {
  /*< private >*/
  GfsOutput parent;

  /*< public >*/
  FttVector pos_min, pos_max;
  gint level;
  GfsVariable * v;
};

typedef struct {
  gint rank, ini;
} DomainData;

#define GFS_OUTPUT_SPECTRA(obj)            GTS_OBJECT_CAST (obj,\
							    GfsOutputSpectra, \
							    gfs_output_spectra_class ())
#define GFS_IS_OUTPUT_SPECTRA(obj)         (gts_object_is_from_class (obj,\
								      gfs_output_spectra_class ()))

GfsOutputClass * gfs_output_spectra_class  (void);

/* GfsOutputSpectraInterface: header */

typedef struct _GfsOutputSpectraInterface                     GfsOutputSpectraInterface;

struct _GfsOutputSpectraInterface {
  /*< private >*/
  GfsOutputSpectra parent;

};

#define GFS_OUTPUT_SPECTRA_INTERFACE(obj)            GTS_OBJECT_CAST (obj,\
							    GfsOutputSpectraInterface, \
							    gfs_output_spectra_interface_class ())
#define GFS_IS_OUTPUT_SPECTRA_INTERFACE(obj)         (gts_object_is_from_class (obj,\
								      gfs_output_spectra_interface_class ()))

GfsOutputClass * gfs_output_spectra_interface_class  (void);


/* GfsOutputEnergySpectra: header */

typedef struct _GfsOutputEnergySpectra                     GfsOutputEnergySpectra;

struct _GfsOutputEnergySpectra {
  /*< private >*/
  GfsOutput parent;

  /*< public >*/
  FttVector pos_min, pos_max;
  gint level;
  GfsVariable * v;
};


#define GFS_OUTPUT_ENERGY_SPECTRA(obj)            GTS_OBJECT_CAST (obj,\
							    GfsOutputEnergySpectra, \
							    gfs_output_energy_spectra_class ())

#define GFS_IS_OUTPUT_ENERGY_SPECTRA(obj)         (gts_object_is_from_class (obj,\
                                                   gfs_output_energy_spectra_class ()))

GfsOutputClass * gfs_output_energy_spectra_class  (void);


typedef struct {
  GfsCommData * CommData;
  GArray * ParallelList;
  gint ndest;
  gint rank, comm_size;
} TransferData;

typedef struct {
  gint coord, np, npaux;
} DirOrderData;

typedef struct {
  GfsVariable * v, * u;
  GfsCartesianGrid * cgd;
  gint levelmax, Ndim;
  gdouble dx;
  FttVector pos_min_global, pos_max_global;
  FttVector pos_min, pos_max;
  DirOrderData dirdata[3];
  /* only for parallel */
  DomainData * domain_ranges;
  TransferData td;
  ptrdiff_t local_n0, local_0_start;
//  gint local_n0, local_0_start;
  fftw_complex *out;
  gboolean plane_vof;
} SpectraData;

gint get_index ( gint i, gint j, gint k, gint np);
