/*
 * =====================================================================================
 *
 *       Filename:  turbulence.c
 *
 *    Description:  Module containing modules for LES simulations
 *
 *        Version:  1.0
 *        Created:  10/28/2011 09:59:57 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dr. Daniel Fuster , dfuster@gmail.com
 *        Company:  Institute Jean Le Rond D'Alembert (Paris)
 *
 * =====================================================================================
 */

#include <stdlib.h>

#include "output.h"
#include "variable.h"
#include "cartesian.h"
#include "mpi_boundary.h"
#include "fft.h"
#include "fluid.h"
#include "utils.h"

/* GfsParallelSpectraData: header */

typedef struct _GfsParallelSpectraData                     GfsParallelSpectraData;

struct _GfsParallelSpectraData {
  /*< private >*/
  GtsObject parent;

  /*< public >*/
  GArray * data; 
  gint node,ndata;
};

typedef struct {
  gint index[3];
  gdouble img[3];
} SendData;


#define GFS_PARALLEL_SPECTRA_DATA(obj)                 GTS_OBJECT_CAST (obj,\
                                                                        GfsParallelSpectraData, \
                                                                        gfs_parallel_spectra_data_class ())
#define GFS_IS_PARALLEL_SPECTRA_DATA(obj)             (gts_object_is_from_class (obj,\
                                                                                 gfs_parallel_spectra_data_class ()))

GtsObjectClass * gfs_parallel_spectra_data_class  (void);

/** \beginobject{GfsParallelSpectraData} */

static void gfs_parallel_spectra_data_read (GtsObject ** o, GtsFile * fp)
{
  GtsObjectClass * klass;
  if (GTS_OBJECT_CLASS (gfs_parallel_spectra_data_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_parallel_spectra_data_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsParallelSpectraData * pd = GFS_PARALLEL_SPECTRA_DATA (*o);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsParallelSpectraData)");
    return;
  }
  klass = gfs_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_parallel_spectra_data_class ())) {
    gts_file_error (fp, "`%s' is not a GfsParallelSpectraData", fp->token->str);
    return;
  }

  gts_file_next_token (fp);
  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting an integer (N)");
    return;
  }
  pd->ndata = atoi (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting an integer (N)");
    return;
  }
  pd->node = atoi (fp->token->str);
  gts_file_next_token (fp);

  gint i,j;
  SendData vals;
  for (i = 0; i < pd->ndata; i++) {
    if (fp->type == '\n')
      gts_file_next_token (fp);
    for (j = 0; j < 3; j++) {
      if (fp->type != GTS_FLOAT && fp->type != GTS_INT) {
        gts_file_error (fp, "expecting a number");
        return;
      }
      vals.index[j] = atoi (fp->token->str);
      gts_file_next_token (fp);
    }
    for (j = 0; j < 3; j++) {
      if (fp->type != GTS_FLOAT && fp->type != GTS_INT) {
        gts_file_error (fp, "expecting a number");
        return;
      }
      vals.img[j] = atof (fp->token->str);
      gts_file_next_token (fp);
    }
    g_array_append_val(pd->data, vals);
  }
}

static void gfs_parallel_spectra_data_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_parallel_spectra_data_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_parallel_spectra_data_class ())->parent_class->write) 
      (o, fp);

  GfsParallelSpectraData * pd = GFS_PARALLEL_SPECTRA_DATA (o);
  gint i;

  fprintf(fp, "%s %i %i \n", o->klass->info.name, pd->ndata, pd->node );
  for (i = 0; i < pd->ndata; i++) {
    SendData vals = g_array_index(pd->data,SendData,i);
    fprintf(fp, "%i %i %i %g %g %g \n", vals.index[0], vals.index[1], vals.index[2], 
                                        vals.img[0],   vals.img[1],   vals.img[2]);
  } 
}

static void gfs_parallel_spectra_data_destroy ( GtsObject * o )
{
   g_array_free(GFS_PARALLEL_SPECTRA_DATA (o)->data, TRUE);

  (* GTS_OBJECT_CLASS (gfs_parallel_spectra_data_class ())->parent_class->destroy) (o);
}

static void gfs_parallel_data_init ( GtsObject * o )
{
    GFS_PARALLEL_SPECTRA_DATA (o)->data = g_array_new(FALSE, TRUE, sizeof (SendData) );
    GFS_PARALLEL_SPECTRA_DATA (o)->ndata = 0;
}


static void gfs_parallel_spectra_data_class_init (GtsObjectClass * klass)
{
  klass->read  = gfs_parallel_spectra_data_read;
  klass->write = gfs_parallel_spectra_data_write;
  klass->destroy = gfs_parallel_spectra_data_destroy; 
}

GtsObjectClass * gfs_parallel_spectra_data_class (void)
{
  static GtsObjectClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_parallel_spectra_data_info = {
      "GfsParallelSpectraData",
      sizeof (GfsParallelSpectraData),
      sizeof (GtsObjectClass),
      (GtsObjectClassInitFunc) gfs_parallel_spectra_data_class_init,
      (GtsObjectInitFunc) gfs_parallel_data_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
        &gfs_parallel_spectra_data_info);
  }

  return klass;
}

static GfsParallelSpectraData * parallel_spectra_data_new (GtsObjectClass * klass)
{
  GfsParallelSpectraData * object;

  object = GFS_PARALLEL_SPECTRA_DATA (gts_object_new (GTS_OBJECT_CLASS (klass)));

  return object;
}

static void free_ParallelSpectraData_List ( GArray * ParallelList, gint ndest )
{
    gint i;
    for (i = 0; i < ndest; i++) {
        GfsParallelSpectraData * pd = ( GfsParallelSpectraData * ) g_array_index(ParallelList,gpointer,i);
        gts_object_destroy (GTS_OBJECT (pd));
    } 
    g_array_free( ParallelList, TRUE);
}

/** \endobject{GfsParallelSpectraData} */

/* GfsInitSpectra: Header */

typedef struct _GfsInitSpectra         GfsInitSpectra;

struct _GfsInitSpectra {
  /*< private >*/
  GfsGenericInit parent;
  GfsVariable * v[FTT_DIMENSION];

  /*< public >*/
  FttVector pos;
  gint level;
  gdouble totE, L;
  gdouble alpha, epsilon, c1, c2, c3, ReL, kmax, seed;
};

#define GFS_INIT_SPECTRA(obj)            GTS_OBJECT_CAST (obj,\
                                                          GfsInitSpectra,\
                                                          gfs_init_spectra_class ())
#define GFS_IS_INIT_SPECTRA(obj)         (gts_object_is_from_class (obj,\
                                                                    gfs_init_spectra_class ()))

GfsGenericInitClass * gfs_init_spectra_class  (void);

/**
 * Initialising variables.
 * \beginobject{GfsInitSpectra}
 */

static void get_domain_ranges ( gint comm_size, gint local_0_start,  DomainData * domain_ranges, gint rank )
{
  gint i;

  GPtrArray * request = g_ptr_array_new ();
  MPI_Request * r;
  for (i = 0; i < comm_size; i++) 
    if (i != rank) {
        r = g_malloc0 (sizeof (MPI_Request));
  	MPI_Isend (&(local_0_start), 1, MPI_INT, i, 1, MPI_COMM_WORLD, r); 
	g_ptr_array_add (request, r);
    }
  

  MPI_Status status;
  for (i = 0; i < comm_size; i++) 
    if (i != rank) {
    	MPI_Recv ( &(domain_ranges[i].ini), 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
    	domain_ranges[i].rank = i; 
    }

  /* Synchronize */
  MPI_Request * r2;
  for (i = 0; i < request->len; i++) {
      r2 = g_ptr_array_index (request, i);
      MPI_Wait (r2, &status);
  }
  g_ptr_array_free (request, TRUE);

}

static gint get_rank_dest ( DomainData  * domain_ranges, gint i, gint comm_size )
{
  /*fixme: improve search algorithm? */
  gint rank_dest = 0;
  while (rank_dest < comm_size) {
    if (domain_ranges[rank_dest].ini > i) {
        rank_dest--;
        return rank_dest;
    } 
    rank_dest++;
  }
  rank_dest--;
  return rank_dest;
}


static void gfs_init_spectra_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_init_spectra_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_init_spectra_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsInitSpectra * v = GFS_INIT_SPECTRA (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));

  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x0", TRUE, &v->pos.x},
    {GTS_DOUBLE, "y0", TRUE, &v->pos.y},
    {GTS_DOUBLE, "z0", TRUE, &v->pos.z},
    {GTS_DOUBLE, "L",  TRUE, &v->L},
    {GTS_DOUBLE, "E" , TRUE, &v->totE},
    {GTS_NONE}
  };
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;


  v->kmax = G_MAXDOUBLE;
  v->ReL = 0.;
  v->seed = 0.;
  GtsFileVariable var2[] = {
    {GTS_DOUBLE, "alpha"  , TRUE, &v->alpha},
    {GTS_DOUBLE, "epsilon", TRUE, &v->epsilon},
    {GTS_DOUBLE, "c1",      TRUE, &v->c1},
    {GTS_DOUBLE, "c2",      TRUE, &v->c2},
    {GTS_DOUBLE, "c3" ,     TRUE, &v->c3},
    {GTS_DOUBLE, "ReL" ,    FALSE, &v->ReL},
    {GTS_DOUBLE, "kmax" ,   FALSE, &v->kmax},
    {GTS_DOUBLE, "seed" ,   FALSE, &v->seed},
    {GTS_NONE}
  };
  gts_file_assign_variables (fp, var2);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type == GTS_INT) {
    v->level = atoi (fp->token->str);
    gts_file_next_token (fp);
  }
  else
    v->level = gfs_domain_depth (domain);

  gint i;
  for (i=0; i < FTT_DIMENSION; i++) {
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a string (v)");
      return;
    }
    domain = GFS_DOMAIN (gfs_object_simulation (*o));
    if (!(v->v[i] = gfs_variable_from_name (domain->variables, fp->token->str))) {
      gts_file_error (fp, "unknown variable `%s'", fp->token->str);
      return;
    }
    gts_file_next_token (fp);
  }


  GfsEvent * event = GFS_EVENT (*o);
  if (event->start < 0. && (event->istep < G_MAXINT || event->step < G_MAXDOUBLE))
    event->start = 0.;
}

static void gfs_init_spectra_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_init_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_init_class ())->parent_class->write) 
      (o, fp);

  GfsInitSpectra * v = GFS_INIT_SPECTRA (o);
  fprintf (fp, " { x0 = %g y0 = %g z0 = %g L = %g E = % g } ",
                v->pos.x, v->pos.y, v->pos.z, v->L, v->totE);
  
  fprintf (fp, " { alpha = %g epsilon = %g c1 = %g c2 = %g c3 = % g ReL = % g kmax = %g seed = %g}",
                v->alpha, v->epsilon, v->c1, v->c2, v->c3, v->ReL, v->kmax, v->seed);

  fprintf (fp, " %d", v->level);

  gint i;
  for (i=0; i < FTT_DIMENSION; i++)
    fprintf (fp, " %s", v->v[i]->name);
  
}

typedef struct {
  GfsCartesianGrid * cgd;
  GfsVariable * u;
} SpectraInitData;

static void get_cell_values (FttCell * cell, SpectraInitData * sid)
{
  gdouble val;
  gdouble vector[3];
  FttVector pos; 
  ftt_cell_pos (cell, &pos);
  vector[0] = pos.x;
  vector[1] = pos.y;
  vector[2] = pos.z;
  if (!gfs_cartesian_grid_interpolate (sid->cgd, vector, &val)) {
    g_warning("Interpolated value not defined (set to zero) \n");
    GFS_VALUE (cell, sid->u) = 0.;
  }

  GFS_VALUE (cell, sid->u) = val;
}

static void data_dest ( TransferData * td, gint rank_dest, SendData sd  )
{
  GfsParallelSpectraData * pd;

  gboolean found = FALSE;
  GSList * list  = td->CommData->data;
  CommListData * nodes;
  while (list && !found) {
    nodes = list->data;
    if ( nodes->node == rank_dest ) found = TRUE;
    else list = list->next;
  }    

  gint i;
  if (found) {
    i=0;
    pd = (GfsParallelSpectraData *) g_array_index(td->ParallelList,gpointer,i);
    while(pd->node != rank_dest){
      i++;
      pd = (GfsParallelSpectraData *) g_array_index(td->ParallelList,gpointer,i);
    }    
    nodes->np++;
  }    
  else {
    pd = parallel_spectra_data_new ( gfs_parallel_spectra_data_class ()); 
    CommListData * newnodes = g_malloc0( sizeof ( CommListData ) ); 
    newnodes->node = pd->node = rank_dest;
    newnodes->src  = td->rank;
    newnodes->np = 1;
    td->CommData->data = g_slist_prepend (td->CommData->data, newnodes);
    g_array_append_val(td->ParallelList,pd);
    td->ndest++;
  }    
  g_array_append_val(pd->data, sd);
  pd->ndata++;
}

static void send_recv_data ( GfsDomain * domain, TransferData * td )
{

  /*Preparing nodes that have to send and receive*/
  /* sending info to all nodes*/
  gint i;
  GPtrArray * request = g_ptr_array_new ();
  GSList * pl = NULL;
  pl = g_slist_prepend (pl, td->CommData);
  for (i = 0; i < td->comm_size; i++) 
    if ( i != td->rank)
	g_ptr_array_add (request, gfs_send_objects (pl, i));
  
  g_slist_free(pl);

  /* receiver code */
  GSList * rcv_list;
  GSList * tmp;
  GfsCommData * cd;

  gint nreceived = 0;
  GArray * rcv_nodes = g_array_new(FALSE, TRUE, sizeof (gint) );
  for (i = 0; i < td->comm_size; i++) {
    if ( i != td->rank) {
    rcv_list = gfs_receive_objects (domain, i);
    cd = GFS_COMM_DATA(rcv_list->data);
    tmp = cd->data;
    while(tmp){
      CommListData * cld = tmp->data;
      if ( td->rank == cld->node) {
        g_array_append_val(rcv_nodes, cld->src);
        nreceived++;
      }
      tmp = tmp->next;
    }
    gts_object_destroy (GTS_OBJECT (cd));
    g_slist_free(rcv_list);
   }
  }
  gts_object_destroy (GTS_OBJECT (td->CommData));

  /* Synchronize */
  for (i = 0; i < request->len; i++)
    gfs_wait (g_ptr_array_index (request, i));
  g_ptr_array_free (request, TRUE);

  /*Second step: Sending and receiving full data*/
  GfsParallelSpectraData * pd;
  request = g_ptr_array_new ();
  for (i = 0; i < td->ndest; i++) {
    pd = (GfsParallelSpectraData *) g_array_index(td->ParallelList,gpointer,i);;
    if ( pd->node != td->rank) {
      tmp = NULL;
      tmp = g_slist_prepend (tmp, pd);
      g_ptr_array_add (request, gfs_send_objects (tmp, pd->node));
      g_slist_free(tmp);
    }
  }
  /*cleaning sent data*/
  free_ParallelSpectraData_List (td->ParallelList, td->ndest);
  td->ParallelList = g_array_new(FALSE, TRUE, sizeof (gpointer) );
  td->ndest = 0;

  for (i = 0; i < nreceived; i++) {
    gint node = g_array_index(rcv_nodes,gint,i);
    if ( node != td->rank) {
    rcv_list = gfs_receive_objects (domain, node);
    while (rcv_list){
      pd = GFS_PARALLEL_SPECTRA_DATA(rcv_list->data);
      if (td->rank == pd->node) {
        g_array_append_val(td->ParallelList,pd);
        td->ndest++;
      }
      else
        gts_object_destroy (GTS_OBJECT (pd));
      rcv_list = rcv_list->next;
    }
   }
  }

  g_array_free(rcv_nodes, TRUE);
  
  /* Synchronize */
  for (i = 0; i < request->len; i++)
    gfs_wait (g_ptr_array_index (request, i));
  g_ptr_array_free (request, TRUE);
  
}

typedef struct{
  TransferData * td;
  fftw_complex ** u;
  gint local_0_start, local_n0;
  gint np;
  DomainData * domain_ranges;
  gdouble deltak, seed;
} InitSpectraData;

static void assign_data_parallel ( InitSpectraData * d )
{
  gint i, j, k, index1,rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  for (j=0; j < d->td->ndest; j++) {
    GfsParallelSpectraData * pd = (GfsParallelSpectraData *) g_array_index(d->td->ParallelList,gpointer,j);
    for (k = 0; k < pd->ndata; k++) {
      SendData vals = g_array_index(pd->data,SendData,k);
      index1 = get_index(vals.index[0]-d->local_0_start,vals.index[1],vals.index[2], d->np);
      for (i=0; i<3; i++){
        d->u[i][index1][0] = 1.;
        d->u[i][index1][1] = vals.img[i];
      }
    }
  }
}

/* function to generate a random symmetric velocity field
 * in the fourier space */
static void generate_vel_field ( GfsDomain * domain, InitSpectraData * d )
{
  gint dim,i,j,k,index1,index2;

  int rank = -1; 
#ifdef HAVE_MPI
  MPI_Comm_size (MPI_COMM_WORLD, &rank);
  d->seed += rank;
#endif

  /* init kz > 0 */
  for (i=d->local_0_start; i < d->local_0_start + d->local_n0; i++) {
    for (j=0; j < d->np; j++) {
      for (k=1; k < ( d->np / 2 + 1 ); k++) {
        index1 = get_index(i-d->local_0_start,j,k,d->np);
        d->u[1][index1][0] = 1.;
        for (dim=0; dim < 3; dim++){
          d->u[dim][index1][0] = 1.;
          srand(d->seed);
          d->u[dim][index1][1] = 100.*(0.5 - (rand() / ((double)RAND_MAX + 1))) ;
        }
      }
    }
  }

  gint rank_dest;
  SendData sd;
  /* init kz = 0 kx > 0 ky > 0*/
  gint istart = MAX(1,d->local_0_start);
  for (i=istart; (i < ( d->np / 2 + 1 ) && i < d->local_0_start + d->local_n0); i++) {
    for (j=1; j <  d->np; j++) {
        index1 = get_index(i-d->local_0_start,j,0,d->np);
        sd.index[0] = d->np-i;
        sd.index[1] = d->np-j;
        sd.index[2] = 0;
        for (dim=0; dim < 3; dim++){
          d->u[dim][index1][0] = 1.;
          srand(d->seed);
          d->u[dim][index1][1] = 100.*(0.5 - (rand() / ((double)RAND_MAX + 1))) ;
          sd.img[dim] = -d->u[dim][index1][1] ;
        }
        rank_dest = get_rank_dest( d->domain_ranges, d->np-i, d->td->comm_size );
        data_dest (d->td, rank_dest, sd);
    }
  }

  /* init kz = 0 kx = 0 ky > 0*/
  if ( d->td->rank == 0 ) {
    for (j=1; j < ( d->np / 2 + 1 ); j++) {
      index1 = get_index(0,j,0,d->np); 
      index2 = get_index(0,d->np-j,0,d->np);
      for (dim=0; dim < 3; dim++){
        d->u[dim][index1][0] = d->u[dim][index2][0] = 1.;
        srand(d->seed);
        d->u[dim][index1][1] = 100.*(0.5 - (rand() / ((double)RAND_MAX + 1))) ;
        d->u[dim][index2][1] = -d->u[dim][index1][1] ;
      }
    }
  }

  /* init kz = 0 kx > 0 ky = 0*/
  for (i=istart; i < ( d->np / 2 + 1 ) && (i < d->local_0_start + d->local_n0) ; i++) {
    index1 = get_index(i-d->local_0_start,0,0,d->np);
    sd.index[0] = d->np-i;
    sd.index[1] = 0;
    sd.index[2] = 0;
    for (dim=0; dim < 3; dim++){
      d->u[dim][index1][0] = 1.;
      srand(d->seed);
      d->u[dim][index1][1] = 100.*(0.5 - (rand() / ((double)RAND_MAX + 1))) ;
      sd.img[dim] = -d->u[dim][index1][1] ;
    }
    rank_dest = get_rank_dest( d->domain_ranges, d->np-i, d->td->comm_size );
    data_dest (d->td, rank_dest, sd);
  }

  send_recv_data (domain, d->td);
  assign_data_parallel (d);
}

/* Solenoidal velocity field in fourier space */
static void solenoidal_vel_field ( InitSpectraData * d, fftw_complex ** usolenoid )
{
  gint dim,i,j,k,index1;
  gdouble kx,ky,kz, kmod2;

  for (i=d->local_0_start; i < d->local_0_start + d->local_n0; i++) {
    if ( i < ( d->np / 2 + 1 ) ) kx = i*d->deltak;
    else kx = ( i - d->np ) * d->deltak;
    for (j=0; j < d->np; j++) {
      if ( j < ( d->np / 2 + 1 ) ) ky = j*d->deltak;
      else ky = ( j - d->np ) * d->deltak;
      for (k=0; k < ( d->np / 2 + 1 ); k++) {
        kz = k*d->deltak;
        kmod2 = pow(kx,2) + pow(ky,2) + pow(kz,2);
        index1 = get_index(i-d->local_0_start,j,k,d->np);
        if ( kmod2 != 0. ) {
          usolenoid[0][index1][0] = (1. - pow(kx,2)/kmod2)* d->u[0][index1][0] 
            - kx*ky/kmod2* d->u[1][index1][0] 
            - kx*kz/kmod2* d->u[2][index1][0];
          usolenoid[0][index1][1] = (1. - pow(kx,2)/kmod2)* d->u[0][index1][1] 
            - kx*ky/kmod2* d->u[1][index1][1] 
            - kx*kz/kmod2* d->u[2][index1][1];

          usolenoid[1][index1][0] =           - ky*kx/kmod2*d->u[0][index1][0] 
            + (1. - pow(ky,2)/kmod2)*d->u[1][index1][0] 
            - ky*kz/kmod2*d->u[2][index1][0];
          usolenoid[1][index1][1] =           - ky*kx/kmod2*d->u[0][index1][1] 
            + (1. - pow(ky,2)/kmod2)*d->u[1][index1][1] 
            - ky*kz/kmod2*d->u[2][index1][1];

          usolenoid[2][index1][0] =           - kz*kx/kmod2*d->u[0][index1][0] 
            - ky*kz/kmod2*d->u[1][index1][0] 
            + (1. - pow(kz,2)/kmod2)*d->u[2][index1][0];
          usolenoid[2][index1][1] =           - kz*kx/kmod2*d->u[0][index1][1] 
            - ky*kz/kmod2*d->u[1][index1][1] 
            + (1. - pow(kz,2)/kmod2)*d->u[2][index1][1];
        }
        else
          for (dim=0; dim < 3; dim++)
            usolenoid[dim][index1][0] =  usolenoid[dim][index1][1] = 0;
      }
    }
  }
}

typedef struct {
  gint local_0_start,local_n0,np;
  fftw_complex *u;
} SpectraEnergyData;

static void spectral_energy ( GfsDomain * domain, SpectraEnergyData * d, fftw_complex ** usolenoid, GfsInitSpectra * ds, gdouble deltak)
{
  gint i,j,k,index1,index2,dim;
  gint knx, kny;
    
  gint nk = 3 * pow( d->np / 2 + 1, 2);
  gdouble * Ek  = g_malloc0 ( nk * sizeof ( gdouble ) );

  for (dim=0; dim < 3; dim++) {
  for (i=d->local_0_start; i < d->local_0_start + d->local_n0; i++) {
    if ( i < ( d->np / 2 + 1 ) ) knx = i;
    else knx = ( d->np -i );
    for (j=0; j < d->np; j++) {
      if ( j < ( d->np / 2 + 1 ) ) kny = j;
      else kny = ( d->np - j );
      index1 =  pow(knx,2) + pow(kny,2);
      index2 = get_index(i-d->local_0_start,j,0,d->np);
      Ek[index1] += 0.5*(pow(usolenoid[dim][index2][0],2) + pow(usolenoid[dim][index2][1],2));
      for (k=1; k < ( d->np / 2 + 1 ); k++) {
        index1 = pow(knx,2) + pow(kny,2) + pow(k,2);
        index2 = get_index(i-d->local_0_start,j,k,d->np);
        Ek[index1] += pow(usolenoid[dim][index2][0],2) + pow(usolenoid[dim][index2][1],2);
      }
    }
  }
  }

    for (i=0; i < nk; i++) 
      gfs_all_reduce (domain, Ek[i], MPI_DOUBLE, MPI_SUM);

    gdouble * cscale =  g_malloc0 ( nk*sizeof ( gdouble ));
    gdouble Ekspectra=0; 
    gdouble kwave, Ei;
    gdouble Lint = pow(ds->totE,3./2.)/ds->epsilon;
    for (i=1; i < nk; i++) {
      kwave = deltak * sqrt((gdouble)i);
      if (Ek[i] != 0.){
        if (ds->ReL != 0) {
        gdouble fl = pow( Lint*kwave/sqrt( pow(Lint*kwave,2) + ds->c1 ),11./3.);
        gdouble feta = exp(- ds->c2*(pow(pow(Lint*kwave*pow(ds->ReL,-3./4.),4) + pow(ds->c3,4) ,0.25)- ds->c3 ));
        if (kwave < ds->kmax )
          Ei = ds->alpha*pow(ds->epsilon,2./3.)*pow(kwave,-5./3.)*fl*feta;
        else
          Ei = 0.;
        }
        else {
          if (kwave < ds->kmax )
            Ei = pow(kwave,2.);
          else
            Ei = 0.;
        }
        cscale[i] = sqrt(Ei/Ek[i]);
        Ekspectra += Ei;
      }
    }

    gdouble cscale2 = sqrt(ds->totE/Ekspectra);

    for (i=d->local_0_start; i < d->local_0_start + d->local_n0; i++) {
      if ( i < ( d->np / 2 + 1 ) ) knx = i;
      else knx = ( d->np - i );
      for (j=0; j < d->np; j++) {
        if ( j < ( d->np / 2 + 1 ) ) kny = j;
        else kny = ( d->np - j );
        for (k=0; k < ( d->np / 2 ) + 1 ; k++) {
          index1 = pow(knx,2) + pow(kny,2) + pow(k,2);
          index2 = get_index(i-d->local_0_start,j,k,d->np);
          for (dim=0; dim < 3; dim++) {
            usolenoid[dim][index2][0] *= cscale2*cscale[index1];
            usolenoid[dim][index2][1] *= cscale2*cscale[index1];
          }
        }
      }
    }

    g_free (cscale);
    g_free(Ek);
}

/*fixme: For now it only works in 3D and if the fftw_mpi.h libraries are installed */
static gboolean gfs_init_spectra_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_spectra_class ())->parent_class)->event) 
      (event, sim)) {

    fftw_mpi_init();

    GfsDomain * domain = GFS_DOMAIN (sim);

    GfsInitSpectra * d = GFS_INIT_SPECTRA(event);
    gdouble Lcorr = d->L;
    gdouble deltak = 2.*M_PI / Lcorr; 
    gint np = pow(2,d->level) ;

    TransferData td;
    td.ParallelList = g_array_new(FALSE, TRUE, sizeof (gpointer) );
    td.ndest = 0;
    td.CommData =  comm_data_new (gfs_comm_data_class ());
    MPI_Comm_rank(MPI_COMM_WORLD,&(td.rank));
    MPI_Comm_size (MPI_COMM_WORLD, &(td.comm_size)); 
    srand(td.rank+sim->time.i);

    /* memory allocation */
    ptrdiff_t alloc_local;
    ptrdiff_t local_n0, local_0_start;

    ptrdiff_t * np_array = g_malloc0( 3 * sizeof(ptrdiff_t) );
    np_array[0] = np;
    np_array[1] = np;
    np_array[2] = (np / 2) + 1;
    ptrdiff_t block0 = np / td.comm_size;
    alloc_local = fftw_mpi_local_size_many(3, np_array, 3, block0,
                                         MPI_COMM_WORLD, &local_n0, &local_0_start);

    ptrdiff_t nalloc = alloc_local*2;

    fftw_complex ** ucomplex = g_malloc ( sizeof(  fftw_complex * ) * 3 );

    gint dim;
    for (dim=0; dim < 3; dim++)
      ucomplex[dim] = fftw_alloc_complex( alloc_local );

    DomainData * domain_ranges =  g_malloc0 ( sizeof( DomainData ) * td.comm_size );
    get_domain_ranges ( td.comm_size, local_0_start, domain_ranges, td.rank);

    /* random velocity field in fourier space */
    InitSpectraData isd = { &td, ucomplex, local_0_start, local_n0, np, domain_ranges, deltak, d->seed };

    generate_vel_field ( domain , &isd );

    fftw_complex ** usolenoid = g_malloc0 ( sizeof(  fftw_complex * ) * 3 );
    for (dim=0; dim < 3; dim++)
      usolenoid[dim] = fftw_alloc_complex( alloc_local );

    solenoidal_vel_field ( &isd , usolenoid );

    /* free initial random field */
    for (dim=0; dim < 3; dim++)
      fftw_free ( ucomplex[dim] );
    g_free( ucomplex );

    /* Fit and rescale velocity field to spectra */
    SpectraEnergyData sed = { local_0_start,local_n0,np };
    spectral_energy ( domain, &sed, usolenoid, d, deltak );

    gint i,j,k,l;
    MPI_Request r;
    MPI_Status s;
    gint * npnodes = g_malloc0( sizeof ( gint ) * td.comm_size );

    /*initialize cgd data */
    SpectraInitData sid;
    sid.cgd = gfs_cartesian_grid_new (gfs_cartesian_grid_class ());
    sid.cgd->N = 3;
    sid.cgd->n = g_malloc0 (sid.cgd->N*sizeof (guint));
    sid.cgd->n[0] = sid.cgd->n[1] = sid.cgd->n[2] = np;
    sid.cgd->v = g_malloc0( sizeof ( gdouble ) * np * np * np  );
    gint maxsize = nalloc;
    gint nlocal;
    for (i = 0; i < td.comm_size; i++) {
      MPI_Isend (&nalloc, 1,  MPI_INT, i, 0, MPI_COMM_WORLD, &r);
      MPI_Recv(&nlocal, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &s);
      npnodes[i] = nlocal;
      maxsize = MAX ( maxsize, npnodes[i]);
    }
    gdouble * out = g_malloc0( sizeof ( gdouble ) * nalloc );
    gdouble * recv  = g_malloc0( sizeof ( gdouble ) * maxsize  );
    sid.cgd->name = g_malloc0 (sid.cgd->N*sizeof (char *));
    sid.cgd->name[0] = g_strdup ("x");
    sid.cgd->name[1] = g_strdup ("y");
    sid.cgd->name[2] = g_strdup ("z");
    sid.cgd->x = g_malloc0 (3*sizeof (gdouble *));
    for (i = 0; i < 3; i++) {
      sid.cgd->x[i] = g_malloc0 (sid.cgd->n[i]*sizeof (gdouble));
      for (j = 0; j < sid.cgd->n[i]; j++)
        sid.cgd->x[i][j] = (&(d->pos.x))[i] + Lcorr * ( (gdouble) j/(sid.cgd->n[i]-1) - 0.5 ) ;
    }

    fftw_plan p;

    gdouble Lbox = gfs_object_simulation (d)->physical_params.L;
    gint local_n0_rcv;

    for (dim=0; dim < FTT_DIMENSION; dim++) {
      p = fftw_mpi_plan_dft_c2r_3d(np, np, np, usolenoid[dim], out,
                               MPI_COMM_WORLD, FFTW_ESTIMATE);
      fftw_execute(p); 

      for (l = 0; l < td.comm_size; l++) {
        MPI_Isend (out, nalloc,  MPI_DOUBLE, l, 0, MPI_COMM_WORLD, &r);
        MPI_Recv (recv,npnodes[l], MPI_DOUBLE, l, 0, MPI_COMM_WORLD, &s);
        MPI_Isend (&local_n0, 1,  MPI_INT, l, 0, MPI_COMM_WORLD, &r);
        MPI_Recv (&local_n0_rcv, 1, MPI_INT, l, 0, MPI_COMM_WORLD, &s);
        for (i = 0; i < local_n0_rcv; i++) 
          for (j=0; j < np; j++) 
            for (k=0; k < np; k++) {
              sid.cgd->v[k + np * ( j + np * ( i + domain_ranges[l].ini ) ) ] = 
                                    recv[k + 2* (np/2 + 1) * ( j + np * i ) ]/Lbox;
            }
      }
      sid.u = d->v[dim];
      gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
                                (FttCellTraverseFunc) get_cell_values, &sid);

    fftw_destroy_plan(p);
    }
    
    fftw_mpi_cleanup();

    /* free memory allocation */
    g_free(npnodes);
    g_free(out); 
    g_free(recv);
    gts_object_destroy (GTS_OBJECT (sid.cgd));
    for (dim=0; dim < 3; dim++)
      fftw_free ( usolenoid[dim] );
    g_free( usolenoid );
    g_free(domain_ranges);
    free_ParallelSpectraData_List (td.ParallelList, td.ndest);

    g_free(np_array);

    return TRUE;
  }
  return FALSE;
}

static void gfs_init_spectra_class_init (GfsGenericInitClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_spectra_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_init_spectra_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_spectra_write;
  /*should I include a destroy method? For variables in v->v[dim]*/
}

GfsGenericInitClass * gfs_init_spectra_class (void)
{
  static GfsGenericInitClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_spectra_info = {
      "GfsInitSpectra",
      sizeof (GfsInitSpectra),
      sizeof (GfsGenericInitClass),
      (GtsObjectClassInitFunc) gfs_init_spectra_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
        &gfs_init_spectra_info);
  }

  return klass;
}

/** \endobject{GfsInitSpectra} */

/* GfsVariableTurbulentViscosity: header */

typedef struct _GfsVariableTurbulentViscosity               GfsVariableTurbulentViscosity;

struct _GfsVariableTurbulentViscosity {
  /*< private >*/
  GfsVariable parent;

  /*< public >*/
  gdouble Cs;
  gint model_type;
};

#define GFS_VARIABLE_TURBULENT_VISCOSITY(obj)            GTS_OBJECT_CAST (obj,\
                                                   GfsVariableTurbulentViscosity,\
                                                   gfs_variable_tracer_class ())
#define GFS_IS_VARIABLE_TURBULENT_VISCOSITY(obj)         (gts_object_is_from_class (obj,\
                                             gfs_variable_turbulent_viscosity_class ()))

GfsVariableClass * gfs_variable_turbulent_viscosity_class  (void);

/**
 * \beginobject{GfsVariableTurbulentViscosity}
 */

typedef struct {
  GfsVariableTurbulentViscosity * v;
  GfsVariable ** u;
} TVData;


static void get_smagorinsky_viscosity (FttCell * cell, TVData * tvd )
{
  gint i, j;
  gdouble s = 0;
  GfsVariable * v = GFS_VARIABLE(GFS_EVENT(tvd->v));
  FttVector g[FTT_DIMENSION];
  gdouble h = ftt_cell_size (cell);
  for (i=0; i<FTT_DIMENSION; i++)
    gfs_cm_gradient ( cell, tvd->u[i], &g[i]);
  for (i=0; i<FTT_DIMENSION; i++) {
  for (j=0; j<FTT_DIMENSION; j++) {
    s += pow(0.5*((&(g[i].x))[j] + (&(g[j].x))[i])/h,2);
  }
  }
  s = sqrt(2.*s);
  GFS_VALUE (cell, v) = pow(tvd->v->Cs*h,2)*s;
}

static void get_sigma_viscosity (FttCell * cell, TVData * tvd )
{
  gint i, j, k;
  GfsVariable * v = GFS_VARIABLE(GFS_EVENT(tvd->v));
  FttVector g[FTT_DIMENSION], g2[FTT_DIMENSION], g22[FTT_DIMENSION];
  gdouble h = ftt_cell_size (cell);

  for (i=0; i<FTT_DIMENSION; i++)
    gfs_cm_gradient ( cell, tvd->u[i], &g[i]);

  for (i=0; i<FTT_DIMENSION; i++) {
  for (j=0; j<FTT_DIMENSION; j++) {
    (&(g2[i].x))[j] = 0;
    for (k=0; k<FTT_DIMENSION; k++) 
      (&(g2[i].x))[j] += (&(g[k].x))[i]*(&(g[k].x))[j]/pow(h,2);
  }
  }

  for (i=0; i<FTT_DIMENSION; i++) {
  for (j=0; j<FTT_DIMENSION; j++) {
    (&(g22[i].x))[j] = 0;
    for (k=0; k<FTT_DIMENSION; k++) 
      (&(g22[i].x))[j] += (&(g2[i].x))[k]*(&(g2[k].x))[j];
  }
  }

  gdouble inv1 = 0., inv2 = 0., inv3;
  for (i=0; i<FTT_DIMENSION; i++){
    inv1 +=  (&(g2[i].x))[i];
    inv2 +=  (&(g22[i].x))[i];
  }

  inv2 = ( pow(inv1,2) - inv2 ) / 2;

#if FTT_2D
  inv3 = g2[0].x*g2[1].y -  g2[0].y*g2[1].x;
#else
  inv3  = (&(g2[0].x))[0]*(&(g2[1].x))[1]*(&(g2[2].x))[2];
  inv3 += (&(g2[0].x))[1]*(&(g2[1].x))[2]*(&(g2[2].x))[0];
  inv3 += (&(g2[0].x))[2]*(&(g2[1].x))[0]*(&(g2[2].x))[1];
  inv3 -= (&(g2[0].x))[2]*(&(g2[1].x))[1]*(&(g2[2].x))[0];
  inv3 -= (&(g2[0].x))[1]*(&(g2[1].x))[0]*(&(g2[2].x))[2];
  inv3 -= (&(g2[0].x))[0]*(&(g2[1].x))[2]*(&(g2[2].x))[1];
#endif

  gdouble alpha1 = pow(inv1,2) / 9. - inv2 / 3.;
  gdouble alpha2 = pow(inv1,3) / 27. - inv1*inv2/6. + inv3/2.;
  if ( alpha1 <= 0 ) {
    GFS_VALUE (cell, v) = 0.;
    return;
  }
  else if ( alpha2 >= pow(alpha1,3./2.) ) {
    GFS_VALUE (cell, v) = 0.;
    return;
  }
  gdouble alpha3 = 1./3.*acos(alpha2/pow(alpha1,3./2.));

  gdouble sigma1 = sqrt(inv1/3. + 2*sqrt(alpha1)*cos(alpha3));
  gdouble sigma2 = sqrt(inv1/3. - 2*sqrt(alpha1)*cos(M_PI/3. + alpha3));
  gdouble sigma3 = sqrt(inv1/3. - 2*sqrt(alpha1)*cos(M_PI/3. - alpha3));

  gdouble Dsigma;
  if (sigma1 != 0.)
    Dsigma = sigma3*(sigma1 - sigma2)*(sigma2 -sigma3)/pow(sigma1,2);
  else
    Dsigma = 0.;

  GFS_VALUE (cell, v) = pow(tvd->v->Cs*h,2)*Dsigma;
}

static gboolean variable_turbulent_viscosity_event (GfsEvent * event,  GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_variable_turbulent_viscosity_class ())->parent_class)->event) (event, sim)) {

  GfsDomain * domain = GFS_DOMAIN (sim);
  TVData tvd = { GFS_VARIABLE_TURBULENT_VISCOSITY(event), gfs_domain_velocity (domain) };
  if (GFS_VARIABLE_TURBULENT_VISCOSITY(event)->model_type == 1)
      gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                               (FttCellTraverseFunc) get_smagorinsky_viscosity, &tvd);
  else
      gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                               (FttCellTraverseFunc) get_sigma_viscosity, &tvd);
  return TRUE;
}
return FALSE;
}

static void variable_turbulent_viscosity_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_turbulent_viscosity_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number");
    return;
  }
  GFS_VARIABLE_TURBULENT_VISCOSITY(*o)->Cs = atof (fp->token->str);
  gts_file_next_token (fp);

}

static void variable_turbulent_viscosity_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_variable_turbulent_viscosity_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %g ", GFS_VARIABLE_TURBULENT_VISCOSITY(o)->Cs);
}

static void variable_turbulent_viscosity_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = variable_turbulent_viscosity_event;
  klass->read = variable_turbulent_viscosity_read;
  klass->write = variable_turbulent_viscosity_write;
}

static void variable_turbulent_viscosity_init (GfsVariableTurbulentViscosity * v)
{
  GFS_VARIABLE (GFS_EVENT(v))->description = g_strdup ("TurbulentViscosity");
  v->model_type = 1;  
}

GfsVariableClass * gfs_variable_turbulent_viscosity_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_variable_turbulent_viscosity_info = {
      "GfsVariableTurbulentViscosity",
      sizeof (GfsVariableTracer),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) variable_turbulent_viscosity_class_init,
      (GtsObjectInitFunc) variable_turbulent_viscosity_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL 
    };   
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()), 
                                  &gfs_variable_turbulent_viscosity_info);
  }

  return klass;
}

/** \endobject{GfsTurbulentViscosity} */

/* Initialize module */

const gchar gfs_module_name[] = "turbulence";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{ 
  gfs_variable_turbulent_viscosity_class ();
  gfs_parallel_spectra_data_class ();
  gfs_init_spectra_class ();
  return NULL; 
}

