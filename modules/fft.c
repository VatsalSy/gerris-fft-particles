/* Gerris - The GNU Flow Solver                       (-*-C-*-)
 * Copyright (C) 2011 Daniel Fuster
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

#include "fft.h"

/**
 * Storage for Parallel commmunication data.
 * \beginobject{GfsCommData}
 */

static void comm_data_read (GtsObject ** o, GtsFile * fp)
{
  GtsObjectClass * klass;
  if (GTS_OBJECT_CLASS (gfs_comm_data_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_comm_data_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsCommData)");
    return;
  }
  klass = gfs_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_comm_data_class ())) {
    gts_file_error (fp, "`%s' is not a GfsCommData", fp->token->str);
    return;
  }

  gts_file_next_token (fp);
  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting an integer (N)");
    return;
  }
  gint n = atoi (fp->token->str);
  GfsCommData * cd = GFS_COMM_DATA(*o);
  CommListData * cld;
  gint i;
  for (i = 0; i < n; i++) {
    gts_file_next_token (fp);
    cld = g_malloc0( sizeof(CommListData) );
    if (fp->type != GTS_INT) {
      gts_file_error (fp, "expecting an integer (N)");
      return;
    }
    cld->node = atoi (fp->token->str);

    gts_file_next_token (fp);
    if (fp->type != GTS_INT) {
      gts_file_error (fp, "expecting an integer (N)");
      return;
    }
    cld->src = atoi (fp->token->str);

    gts_file_next_token (fp);
    if (fp->type != GTS_INT) {
      gts_file_error (fp, "expecting an integer (N)");
      return;
    }
    cld->np = atoi (fp->token->str);
    cd->data = g_slist_prepend (cd->data, cld);
  }

}

static void comm_data_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_comm_data_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_comm_data_class ())->parent_class->write) 
      (o, fp);

  GfsCommData * cd = GFS_COMM_DATA(o);
  GSList * list = cd->data;
  fprintf (fp, "%s %i ", o->klass->info.name, g_slist_length(list));
  CommListData * cld;
  while(list) {
    cld = list->data;
    fprintf (fp, "%i %i %i ", cld->node, cld->src, cld->np);
    list = list->next;
  }
}

static void comm_data_destroy (GtsObject * o)
{
  GfsCommData * cd = GFS_COMM_DATA(o);
  GSList * list = cd->data;
  while(list){
    CommListData * cld = list->data;
    g_free (cld);
    list = list->next;
  }
  g_slist_free(cd->data);

  (* GTS_OBJECT_CLASS (gfs_comm_data_class ())->parent_class->destroy) (o); 
}


static void comm_data_class_init (GtsObjectClass * klass)
{
  klass->read  = comm_data_read;
  klass->write = comm_data_write;
  klass->destroy = comm_data_destroy;
}

GtsObjectClass * gfs_comm_data_class (void)
{
  static GtsObjectClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_comm_data_info = {
      "GfsCommData",
      sizeof (GfsCommData),
      sizeof (GtsObjectClass),
      (GtsObjectClassInitFunc) comm_data_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
        &gfs_comm_data_info);
  }

  return klass;
}

GfsCommData * comm_data_new (GtsObjectClass * klass)
{
  GfsCommData * object;

  object = GFS_COMM_DATA (gts_object_new (GTS_OBJECT_CLASS (klass)));

  return object;
}

/** \endobject{GfsCommData} */

/* GfsParallelData: header */

typedef struct _GfsParallelData                     GfsParallelData;

typedef struct {
    gint i,j,k;
    gdouble val;
} SendData;

struct _GfsParallelData {
  /*< private >*/
  GtsObject parent;

  /*< public >*/
  GArray * data; 
  gint node,ndata;
  FttVector pos_min, pos_max;
};

#define GFS_PARALLEL_DATA(obj)                 GTS_OBJECT_CAST (obj,\
                                                                GfsParallelData, \
                                                                gfs_parallel_data_class ())
#define GFS_IS_PARALLEL_DATA(obj)             (gts_object_is_from_class (obj,\
                                                                         gfs_parallel_data_class ()))

GtsObjectClass * gfs_parallel_data_class  (void);

/** \beginobject{GfsParallelData} */

static void gfs_parallel_data_read (GtsObject ** o, GtsFile * fp)
{
  GtsObjectClass * klass;
  if (GTS_OBJECT_CLASS (gfs_parallel_data_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_parallel_data_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsParallelData * pd = GFS_PARALLEL_DATA (*o);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsParallelData)");
    return;
  }
  klass = gfs_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_parallel_data_class ())) {
    gts_file_error (fp, "`%s' is not a GfsParallelData", fp->token->str);
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

  while (fp->type == '\n') 
    gts_file_next_token (fp);

  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x0", TRUE, &pd->pos_min.x},
    {GTS_DOUBLE, "y0", TRUE, &pd->pos_min.y},
    {GTS_DOUBLE, "z0", TRUE, &pd->pos_min.z},
    {GTS_DOUBLE, "x1", TRUE, &pd->pos_max.x},
    {GTS_DOUBLE, "y1", TRUE, &pd->pos_max.y},
    {GTS_DOUBLE, "z1", TRUE, &pd->pos_max.z},
    {GTS_NONE}
  };
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  gint i,j;
  SendData vals; 
  for (i = 0; i < pd->ndata; i++) {
    if (fp->type == '\n')
      gts_file_next_token (fp);
    for (j = 0; j < 3; j++) {
      if (fp->type != GTS_FLOAT && fp->type != GTS_INT) {
        gts_file_error (fp, "expecting a number at %i %i \n", i, j);
        return;
      }
      (&(vals.i))[j] = atoi (fp->token->str);
      gts_file_next_token (fp);
    }
    if (fp->type != GTS_FLOAT && fp->type != GTS_INT) {
        gts_file_error (fp, "expecting a number %i \n", i);
        return;
    }
    vals.val =  atof (fp->token->str);
    g_array_append_val(pd->data, vals);
    gts_file_next_token (fp);
  }
}

static void gfs_parallel_data_write (GtsObject * o, FILE * fp)
{
  if (GTS_OBJECT_CLASS (gfs_parallel_data_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_parallel_data_class ())->parent_class->write) 
      (o, fp);

  GfsParallelData * pd = GFS_PARALLEL_DATA (o);

  fprintf(fp, "%s %i %i \n", o->klass->info.name, pd->ndata, pd->node);
  fprintf(fp, " { x0 = %g y0 = %g z0 = %g x1 = %g y1 = %g z1 = %g } \n", 
      pd->pos_min.x, pd->pos_min.y, pd->pos_min.z, pd->pos_max.x, pd->pos_max.y, pd->pos_max.z );
  gint i;
   for (i = 0; i < pd->ndata; i++) {
    SendData vals = g_array_index(pd->data,SendData,i);
    fprintf(fp, "%i %i %i %g \n", vals.i, vals.j, vals.k, vals.val);
  } 
}

static void gfs_parallel_data_destroy ( GtsObject * o )
{
  g_array_free(GFS_PARALLEL_DATA (o)->data, TRUE);

  (* GTS_OBJECT_CLASS (gfs_parallel_data_class ())->parent_class->destroy) (o);
}

static void gfs_parallel_data_init ( GtsObject * o )
{
    GFS_PARALLEL_DATA (o)->data = g_array_new(FALSE, TRUE, sizeof (SendData) );
    GFS_PARALLEL_DATA (o)->ndata = 0;
}

static void gfs_parallel_data_class_init (GtsObjectClass * klass)
{
  klass->read  = gfs_parallel_data_read;
  klass->write = gfs_parallel_data_write;
  klass->destroy = gfs_parallel_data_destroy; 
}

GtsObjectClass * gfs_parallel_data_class (void)
{
  static GtsObjectClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_parallel_data_info = {
      "GfsParallelData",
      sizeof (GfsParallelData),
      sizeof (GtsObjectClass),
      (GtsObjectClassInitFunc) gfs_parallel_data_class_init,
      (GtsObjectInitFunc) gfs_parallel_data_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
        &gfs_parallel_data_info);
  }

  return klass;
}

static GfsParallelData * parallel_data_new (GtsObjectClass * klass)
{
  GfsParallelData * object;

  object = GFS_PARALLEL_DATA (gts_object_new (GTS_OBJECT_CLASS (klass)));

  return object;
}

static void free_ParallelData_List ( GArray * ParallelList, gint ndest )
{
    gint i;
    for (i = 0; i < ndest; i++) {
        GfsParallelData * pd = ( GfsParallelData * ) g_array_index(ParallelList,gpointer,i);
        gts_object_destroy (GTS_OBJECT (pd));
    } 
    g_array_free( ParallelList, TRUE);
}

/** \endobject{GfsParallelData} */

/** \beginobject{GfsOutputSpectra} */

static gboolean inside_domain (FttCell * cell,  gpointer data)
{
  SpectraData * d = data;
  FttVector pos; 
  ftt_cell_pos (cell, &pos);
  gint i;
  gboolean inside = TRUE;
  for (i = 0; i < d->Ndim; i++) 
    inside = inside && ((&(pos.x))[i] >= (&(d->pos_min_global.x))[i]) 
                    && ((&(pos.x))[i] <= (&(d->pos_max_global.x))[i]);
  
  if ( inside ) return TRUE;

  return FALSE;
}

static void array_min ( FttVector * pos1, FttVector * pos2 )
{
  pos1->x = MIN ( pos1->x, pos2->x) ;
  pos1->y = MIN ( pos1->y, pos2->y) ;
  pos1->z = MIN ( pos1->z, pos2->z) ;
}

static void array_max ( FttVector * pos1, FttVector * pos2 )
{
  pos1->x = MAX ( pos1->x, pos2->x) ;
  pos1->y = MAX ( pos1->y, pos2->y) ;
  pos1->z = MAX ( pos1->z, pos2->z) ;
}

static void coarsest_leaf (FttCell * cell,  SpectraData * d)
{
  d->levelmax = MAX( d->levelmax, ftt_cell_level(cell));
  FttVector pos; 
  ftt_cell_pos (cell, &pos);
  array_min (&(d->pos_min), &pos);
  array_max (&(d->pos_max), &pos);
}

static gint get_coord ( gdouble v1, gdouble v2, gdouble dx )
{
  gint i = (gint) (( v1 - v2 ) / dx);
  return i;
}

static gint get_index_matrix ( gint i, gint j, gint k, SpectraData * d)
{
  if (d->Ndim == 1 )
    return k+d->dirdata[2].npaux*(i*d->dirdata[1].npaux+j);
  else  if (d->Ndim == 2 )
    return k+d->dirdata[2].npaux*(i*2*d->dirdata[1].npaux+j);
    else
      return k+2*d->dirdata[2].npaux*(i*d->dirdata[1].npaux+j);
}

static void get_data (FttCell * cell,  SpectraData * d)
{
  if (ftt_cell_level(cell) == d->levelmax) { 
  FttVector pos; 
  ftt_cell_pos (cell, &pos);

  gint i = get_coord((&(pos.x))[d->dirdata[0].coord],
                      (&(d->pos_min.x))[d->dirdata[0].coord], d->dx);
  gint j = get_coord((&(pos.x))[d->dirdata[1].coord],
                      (&(d->pos_min.x))[d->dirdata[1].coord], d->dx);
  gint k = get_coord((&(pos.x))[d->dirdata[2].coord],
                      (&(d->pos_min.x))[d->dirdata[2].coord], d->dx);

  gint ntot = d->dirdata[2].np*d->dirdata[1].np*d->dirdata[0].np;
  gint index = get_index_matrix (i,j,k,d);
  d->cgd->v[index] = GFS_VALUE(cell,d->u)/ntot;
  }
}

static void get_data_parallel ( SpectraData * d )
{
  gint i,j;
  for (j=0; j < d->td.ndest; j++) {
    GfsParallelData * pd = (GfsParallelData *) g_array_index(d->td.ParallelList,gpointer,j);
    for (i=0; i < pd->ndata; i++){
      SendData vals = g_array_index(pd->data,SendData,i);;

      gint i2 = (&(vals.i))[d->dirdata[0].coord]-d->local_0_start;
      gint j2 = (&(vals.i))[d->dirdata[1].coord];
      gint k2 = (&(vals.i))[d->dirdata[2].coord];

      gint ntot;

      ntot = d->dirdata[2].np*d->dirdata[1].np*d->dirdata[0].np;
      gint index = get_index_matrix (i2,j2,k2,d);
      d->cgd->v[index] = vals.val/ntot; 
    }
  }
}

static void get_deep_level( GfsDomain * domain, SpectraData * d)
{
  d->pos_max = d->pos_min_global;
  d->pos_min = d->pos_max_global;
  gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                                      (FttCellTraverseFunc) coarsest_leaf, d, inside_domain, d);          
  return;
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

static void pack_data (FttCell * cell,  SpectraData * d, SendData * vals, FttVector * pos )
{

    GfsParallelData * pd;
    gint rank_dest = get_rank_dest (d->domain_ranges, vals->i, d->td.comm_size);

    gboolean found = FALSE;
    GfsCommData * cd = d->td.CommData;
    GSList * list  = cd->data;
    gint i;

    CommListData * nodes;
    while (list && !found) {
      nodes = list->data;
      if ( nodes->node == rank_dest ) found = TRUE;
      else list = list->next;
    }    

    if (found) {
      i=0;
      pd = (GfsParallelData *) g_array_index(d->td.ParallelList,gpointer,i);
      while(pd->node != rank_dest){
        i++;
        pd = (GfsParallelData *) g_array_index(d->td.ParallelList,gpointer,i);
      }    
      nodes->np++;
    }    
    else {
      pd = parallel_data_new ( gfs_parallel_data_class ()); 
      pd->pos_max = d->pos_min_global;
      pd->pos_min = d->pos_max_global;
      CommListData * newnodes = g_malloc0( sizeof ( CommListData ) ); 
      newnodes->node = pd->node = rank_dest;
      newnodes->src  = d->td.rank;
      newnodes->np = 1;
      cd->data = g_slist_prepend (cd->data, newnodes);
      g_array_append_val(d->td.ParallelList,pd);
      d->td.ndest++;
    }    
    g_array_append_val(pd->data, *vals);
    pd->ndata++;
    array_min (&(pd->pos_min), pos);
    array_max (&(pd->pos_max), pos);

}


static void data_dest (FttCell * cell,  SpectraData * d )
{
  if (ftt_cell_level(cell) == d->levelmax ) {
    FttVector pos; 
    ftt_cell_pos (cell, &pos);

    SendData vals;
    gint i;
    for ( i = 0; i < FTT_DIMENSION; i++)
      (&(vals.i))[i] = get_coord((&(pos.x))[d->dirdata[i].coord],
                                 (&(d->pos_min_global.x))[d->dirdata[i].coord],d->dx);
    vals.val = GFS_VALUE(cell,d->u);

    pack_data(cell, d, &vals, &pos);
  }
}

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
  GfsParallelData * pd;
  request = g_ptr_array_new ();
  for (i = 0; i < td->ndest; i++) {
    pd = (GfsParallelData *) g_array_index(td->ParallelList,gpointer,i);;
    if ( pd->node != td->rank) {
    tmp = NULL;
    tmp = g_slist_prepend (tmp, pd);
    g_ptr_array_add (request, gfs_send_objects (tmp, pd->node));
    g_slist_free(tmp);
  }
  }
  
  /*cleaning sent data*/
  free_ParallelData_List (td->ParallelList, td->ndest);
  td->ParallelList = g_array_new(FALSE, TRUE, sizeof (gpointer) );
  td->ndest = 0;

  for (i = 0; i < nreceived; i++) {
    gint node = g_array_index(rcv_nodes,gint,i);
    if ( node != td->rank) {
    rcv_list = gfs_receive_objects (domain, node);
    while (rcv_list){
      pd = GFS_PARALLEL_DATA(rcv_list->data);
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


static void node_communication ( GfsDomain * domain, SpectraData * d )
{

  MPI_Comm_rank(MPI_COMM_WORLD,&(d->td.rank));
  MPI_Comm_size (MPI_COMM_WORLD, &(d->td.comm_size)); 

  /*Creating sending structures*/
  d->td.CommData = comm_data_new (gfs_comm_data_class ());

  d->domain_ranges = g_malloc ( sizeof( DomainData ) * d->td.comm_size );
  get_domain_ranges(d->td.comm_size, d->local_0_start, d->domain_ranges, d->td.rank);

  gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL,  -1,
                                      (FttCellTraverseFunc) data_dest, d, inside_domain, d);          

  send_recv_data( domain, &(d->td));
}

static gdouble position (FttCell * cell, GfsVariable * v, FttVector * posref) 
{
  FttVector p;
  FttCell * parent;
  gdouble height=0.;

  if (gfs_vof_center (cell, GFS_VARIABLE_TRACER_VOF (v), &p)) {
//    if (FTT_CELL_IS_LEAF (cell))
      height = (&p.x)[2];
/*    else {
      FttCell * cellaux = NULL;

      FttCellChildren child;
      ftt_cell_children(cell,&child);

      gint i;
      gdouble mindist = G_MAXDOUBLE;
      for (i = 0; i < FTT_CELLS; i++) {
        if (child.c[i]) {
          FttVector poschild;
          ftt_cell_pos(child.c[i],&poschild);
          gdouble dist = ftt_vector_distance(&poschild,posref);
          if (GFS_VALUE (child.c[i], v) > 0. && GFS_VALUE (child.c[i], v) < 1. && dist < mindist) {
            mindist = dist; 
            cellaux = child.c[i];
          }
        }
      } 

      if (cellaux)
        height = position (cellaux, v, posref);
      else
        height = (&p.x)[2];
    }*/
  }
  else {
    parent = ftt_cell_parent(cell);
    if (parent)
      height = position (parent, v, posref);
  }

  return height;
}


static void fill_interpolated_parallel_matrix ( GfsDomain * domain, SpectraData * d )
{
  gint i,j,k;

  FttVector pos;
  FttCell * cell;

  d->cgd->x = g_malloc0 (3*sizeof (gdouble *));
   for (i = 0; i < 3; i++) {
     d->cgd->x[i] = g_malloc (d->cgd->n[i]*sizeof (gdouble));
     for (j = 0; j < d->cgd->n[i]; j++) 
        d->cgd->x[i][j] = (&(d->pos_min_global.x))[i] + j*d->dx;
  }   

   gdouble val;

  SendData vals;
  gdouble avg = 0.;
  gint np = 0;

  for (i = 0; i < d->cgd->n[0]; i++)
    for (j = 0; j < d->cgd->n[1]; j++)
      for (k = 0; k < d->cgd->n[2]; k++) {
        pos.x = d->cgd->x[0][i];
        pos.y = d->cgd->x[1][j];
        pos.z = d->cgd->x[2][k];

        cell = gfs_domain_locate (domain, pos, -1, NULL);
        if (cell) {
          if (d->plane_vof)
            val = position(cell,d->u,&pos);
          else
            val = GFS_VALUE(cell,d->u);
          avg += val;
          np++;
        }
      }

  gfs_all_reduce (domain, avg, MPI_DOUBLE, MPI_SUM);
  gfs_all_reduce (domain, np, MPI_INT, MPI_SUM);
  avg /= np;

  for (i = 0; i < d->cgd->n[0]; i++)
    for (j = 0; j < d->cgd->n[1]; j++)
      for (k = 0; k < d->cgd->n[2]; k++) {
        pos.x = d->cgd->x[0][i];
        pos.y = d->cgd->x[1][j];
        pos.z = d->cgd->x[2][k];

        cell = gfs_domain_locate (domain, pos, -1, NULL);
        if (cell) {
          if (d->plane_vof)
            val = position(cell,d->u,&pos)-avg;
          else
            val = GFS_VALUE(cell,d->u)-avg;
          vals.i = i;
          vals.j = j;
          vals.k = k;
          vals.val = val;
          pack_data(cell, d, &vals, &pos);
        }
      }
}

static void node_communication_from_matrix ( GfsDomain * domain, SpectraData * d )
{

  MPI_Comm_rank(MPI_COMM_WORLD,&(d->td.rank));
  MPI_Comm_size (MPI_COMM_WORLD, &(d->td.comm_size)); 

  /*Creating sending structures*/
  d->td.CommData = comm_data_new (gfs_comm_data_class ());

  d->domain_ranges = g_malloc ( sizeof( DomainData ) * d->td.comm_size );
  get_domain_ranges(d->td.comm_size, d->local_0_start, d->domain_ranges, d->td.rank);
  fill_interpolated_parallel_matrix ( domain, d );
  send_recv_data( domain, &(d->td));
}

static gint substract ( gpointer a, gpointer b)
{
  return ((DirOrderData *)b)->np - ((DirOrderData *)a)->np;
}

static void order_array ( SpectraData * d )
{
  DirOrderData dirs[3];
  gint i;
  GArray * dir_order = g_array_new(FALSE, FALSE, sizeof (DirOrderData) );
  for (i = 0; i < 3; i++) {
    dirs[i].coord = i;
    dirs[i].np    =  ABS((&(d->pos_max_global.x))[i] - (&(d->pos_min_global.x))[i])/d->dx + 1;
    dirs[i].npaux = dirs[i].np;
    g_array_append_val( dir_order,  dirs[i]);
  }
  g_array_sort ( dir_order, (GCompareFunc)substract);
  d->dirdata[0] = g_array_index( dir_order, DirOrderData, 0);
  d->dirdata[1] = g_array_index( dir_order, DirOrderData, 1);
  d->dirdata[2] = g_array_index( dir_order, DirOrderData, 2);
  g_array_free ( dir_order, TRUE );
  /*correcting aux number of point for allocation purposes*/
  d->dirdata[d->Ndim-1].npaux =  ( d->dirdata[d->Ndim-1].np / 2 ) + 1;

  return;
}

static void fill_interpolated_cartesian_matrix ( GfsDomain * domain, SpectraData * d )
{
  gint i,j,k, index;

  FttVector pos;
  FttCell * cell;

  d->cgd->x = g_malloc0 (3*sizeof (gdouble *));
   for (i = 0; i < 3; i++) {
     d->cgd->x[i] = g_malloc (d->cgd->n[i]*sizeof (gdouble));
     for (j = 0; j < d->cgd->n[i]; j++) 
        d->cgd->x[i][j] = (&(d->pos_min_global.x))[i] + j*d->dx;
  }   

   gdouble val;
   gdouble avg = 0.;
   gint np = 0;

  for (i = 0; i < d->cgd->n[0]; i++)
    for (j = 0; j < d->cgd->n[1]; j++)
      for (k = 0; k < d->cgd->n[2]; k++) {
        pos.x = d->cgd->x[0][i];
        pos.y = d->cgd->x[1][j];
        pos.z = d->cgd->x[2][k];

        cell = gfs_domain_locate (domain, pos, -1, NULL);
        if (cell) {
          if (d->plane_vof)
            val = position(cell,d->u,&pos);
          else
            val = GFS_VALUE(cell,d->u);
          index = get_index_matrix (i,j,k,d);
          d->cgd->v[index] = val;
          avg += val;
          np++;
        }
      }

  gfs_all_reduce (domain, avg, MPI_DOUBLE, MPI_SUM);
  avg /= np;

  for (i = 0; i < d->cgd->n[0]; i++)
    for (j = 0; j < d->cgd->n[1]; j++)
      for (k = 0; k < d->cgd->n[2]; k++) {
        pos.x = d->cgd->x[0][i];
        pos.y = d->cgd->x[1][j];
        pos.z = d->cgd->x[2][k];

        cell = gfs_domain_locate (domain, pos, -1, NULL);
        if (cell) {
          index = get_index_matrix (i,j,k,d);
          d->cgd->v[index] -= avg;
          d->cgd->v[index] /= np;
        }
      }
}

typedef struct {
  GfsDomain * domain;
  gdouble val, vol;
  GfsVariable * v, * u;
} AvgData;

static void add_data (FttCell * cell,  AvgData * d)
{
  gdouble vol = gfs_cell_volume (cell, d->domain);
  d->val += vol*GFS_VALUE (cell, d->v);
  d->vol += vol;
}

static void substract_avg (FttCell * cell,  AvgData * d)
{
  GFS_VALUE (cell, d->u) = GFS_VALUE (cell, d->v) - d->val;
}

static void substract_average ( GfsDomain * domain, SpectraData * d )
{
  AvgData ad = {domain, 0, 0, d->v, d->u };
  gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
                                     (FttCellTraverseFunc) add_data, &ad, inside_domain, d);          
  gfs_all_reduce (domain, ad.val, MPI_DOUBLE, MPI_SUM);
  gfs_all_reduce (domain, ad.vol, MPI_DOUBLE, MPI_SUM);
  g_assert( ad.vol > 0. );
  ad.val = ad.val / ad.vol;
  gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
                                     (FttCellTraverseFunc) substract_avg, &ad, inside_domain, d);          
}

/*fixme: do something when there is no data on the processor*/
static void get_domain_limits ( GfsDomain * domain, SpectraData * d )
{
  gint i;
#ifdef HAVE_MPI
  gfs_all_reduce (domain, d->levelmax, MPI_INT, MPI_MIN);
#endif

  d->dx = ftt_level_size(d->levelmax)* gfs_object_simulation (d->v)->physical_params.L;

  for (i = 0; i < 3; i++) {
    (&(d->pos_max_global.x))[i] = (&(d->pos_max.x))[i];
    (&(d->pos_min_global.x))[i] = (&(d->pos_min.x))[i];
#ifdef HAVE_MPI  
    gfs_all_reduce (domain, (&(d->pos_max_global.x))[i], MPI_DOUBLE, MPI_MAX);
    gfs_all_reduce (domain, (&(d->pos_min_global.x))[i], MPI_DOUBLE, MPI_MIN);
#endif
  }
}

static void allocate_cartesian_matrix ( GfsDomain * domain, SpectraData * d )
{
  gint i;
  if (d->levelmax == 0) {
//    fill_interpolated_cartesian_matrix(domain, d);
    g_error("interpolation not implemented yet");
    return;
  }

  /* number of local points in each direction */
  for (i = 0; i < 3; i++) {
    if ( ( (&(d->pos_max.x))[i] -  (&(d->pos_min.x))[i] ) == 0 )
      d->cgd->n[i] = 1;
    else
      d->cgd->n[i] = (gint)(( (&(d->pos_max.x))[i] -  (&(d->pos_min.x))[i] )/d->dx) + 1;
  }

  order_array(d);

/* allocate local data size */
  ptrdiff_t alloc_local;
#ifdef HAVE_MPI  
  alloc_local = fftw_mpi_local_size_3d(d->dirdata[0].npaux, d->dirdata[1].npaux,
                                       d->dirdata[2].npaux, MPI_COMM_WORLD, 
                                       &(d->local_n0), &(d->local_0_start));
   
  d->out = g_malloc( sizeof( fftw_complex ) * alloc_local );
  d->cgd->v = g_malloc0( sizeof ( gdouble ) * 2 * alloc_local );
#else
  d->local_0_start = d->local_n0 = 0;
  alloc_local = d->dirdata[0].np * d->dirdata[1].np * d->dirdata[2].np;
  d->out = fftw_alloc_complex ( alloc_local ) ;
  d->cgd->v = fftw_alloc_real ( 2 * alloc_local ) ; 
#endif
}

static void fill_cartesian_matrix ( GfsDomain * domain, SpectraData * d )
{
  d->u = gfs_temporary_variable (domain);
  substract_average(domain, d);
  /*obtain global data*/
#ifdef HAVE_MPI  

  d->td.ParallelList = g_array_new(FALSE, TRUE, sizeof (gpointer) );
  d->td.ndest = 0;
  node_communication ( domain, d);
  
  if (d->td.comm_size > 1) {
    GfsParallelData * pd;
    d->pos_min = d->pos_max_global;
    d->pos_max = d->pos_min_global;
    gint i;
    for (i = 0; i < d->td.ndest; i++) {
      pd = ( GfsParallelData * ) g_array_index(d->td.ParallelList,gpointer,i);
      array_min (&(d->pos_min), &(pd->pos_min));
      array_max (&(d->pos_max), &(pd->pos_max));
  }
    get_data_parallel (d);
  } else
    gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
                                        (FttCellTraverseFunc) get_data, d, inside_domain, d);          

  free_ParallelData_List (d->td.ParallelList, d->td.ndest);

#else
  gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
                                     (FttCellTraverseFunc) get_data, d, inside_domain, d);          
#endif
    
  gts_object_destroy (GTS_OBJECT (d->u));

}

static void fill_cartesian_matrix_plane ( GfsDomain * domain, SpectraData * d )
{
  /*obtain global data*/
#ifdef HAVE_MPI  

  d->td.ParallelList = g_array_new(FALSE, TRUE, sizeof (gpointer) );
  d->td.ndest = 0;
  node_communication_from_matrix ( domain, d);
  
  if (d->td.comm_size > 1) {
    GfsParallelData * pd;
    d->pos_min = d->pos_max_global;
    d->pos_max = d->pos_min_global;
    gint i;
    for (i = 0; i < d->td.ndest; i++) {
      pd = ( GfsParallelData * ) g_array_index(d->td.ParallelList,gpointer,i);
      array_min (&(d->pos_min), &(pd->pos_min));
      array_max (&(d->pos_max), &(pd->pos_max));
    }
    get_data_parallel (d);
  } else 
    fill_interpolated_cartesian_matrix (domain, d);


  free_ParallelData_List (d->td.ParallelList, d->td.ndest);
#else
  fill_interpolated_cartesian_matrix (domain, d);
#endif
}

static FttVector init_kmax ( SpectraData * d )
{
  FttVector kmax;
  gint i;

  for (i = 0; i < 3; i++) {
    gdouble L = ABS((&(d->pos_max_global.x))[i] -  (&(d->pos_min_global.x))[i]);
    if (L != 0) 
      (&(kmax.x))[i] = 2.*M_PI/L; 
    else
      (&(kmax.x))[i] = 0;
  }

  return kmax;
}

static void write_spectra ( FILE * fp, SpectraData * d )
{
  gint i,j,l; 
  gint aux;
  FttVector kmax =  init_kmax ( d );
  FttVector k;
  k.x = k.y = k.z = 0;
  gint np = 1;
  for ( i = 0; i < d->Ndim; i++ )
    np *= d->dirdata[i].np;
  fprintf(fp, "# %i \n", np);
  fputs ("# 1:kx 2:ky 3:kz 4:real 5:img\n", fp);

  for ( i = d->local_0_start; i < (d->local_0_start + d->local_n0); i++ ) {
    if ( i < d->dirdata[0].np/2 + 1 )  (&(k.x))[d->dirdata[0].coord] = (&(kmax.x))[d->dirdata[0].coord]*i;
    else { 
      aux = i-d->dirdata[0].np;
      (&(k.x))[d->dirdata[0].coord] = (&(kmax.x))[d->dirdata[0].coord]*aux;
    }
    for ( j = 0; j < d->dirdata[1].np; j++ ) {
      if ( j < (d->dirdata[1].np/2) + 1 ) (&(k.x))[d->dirdata[1].coord] = (&(kmax.x))[d->dirdata[1].coord]*j;
      else {
        aux = j- d->dirdata[1].np;
        (&(k.x))[d->dirdata[1].coord] = (&(kmax.x))[d->dirdata[1].coord]*aux;
      }
      for ( l = 0; l < d->dirdata[2].npaux; l++ ) {
        (&(k.x))[d->dirdata[2].coord] = (&(kmax.x))[d->dirdata[2].coord]*l;
        fprintf (fp, "%g %g %g %g %g\n", 
            k.x, k.y, k.z , 
            d->out[l+d->dirdata[2].npaux*((i-d->local_0_start)*d->dirdata[1].np+j)][0]
            *gfs_object_simulation (d->v)->physical_params.L, 
            d->out[l+d->dirdata[2].npaux*((i-d->local_0_start)*d->dirdata[1].np+j)][1]
            *gfs_object_simulation (d->v)->physical_params.L);
      }
    }
  }
}

static fftw_plan get_fftw_plan ( SpectraData * d )
{
  fftw_plan p;

#ifdef HAVE_MPI  
  p = fftw_mpi_plan_dft_r2c_3d( d->dirdata[0].np, d->dirdata[1].np,  d->dirdata[2].np, d->cgd->v, d->out, MPI_COMM_WORLD, FFTW_ESTIMATE);
#else
  p = fftw_plan_dft_r2c_3d(d->dirdata[0].np, d->dirdata[1].np, d->dirdata[2].np, d->cgd->v, d->out, FFTW_ESTIMATE);
#endif
      
  return p;
}


static gboolean output_spectra_event (GfsEvent * event, 
    GfsSimulation * sim) 
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_spectra_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsOutputSpectra * v = GFS_OUTPUT_SPECTRA (event);
    fftw_plan p;

  SpectraData d;
  d.v = v->v;
  d.u = v->v;
  d.levelmax = 0;

  d.cgd = gfs_cartesian_grid_new (gfs_cartesian_grid_class ()); 

  /* number of dims of the fft */
  d.cgd->N = 3;
  d.cgd->n = g_malloc ( 3 * sizeof ( guint ) );  

  d.Ndim = 3;
  gint realdim = 0;
  gint i,k=1;
  for (i = 0; i < 3; i++) {
    if (((&(v->pos_max.x))[i] - (&(v->pos_min.x))[i]) != 0) {
      realdim++;
      k++;
    }
  }
  g_assert ( realdim > 0);

  d.pos_min_global = v->pos_min;
  d.pos_max_global = v->pos_max;
  get_deep_level (domain, &d);
  d.levelmax = MIN( d.levelmax, v->level);

#ifdef HAVE_MPI  
    fftw_mpi_init();
#endif

    get_domain_limits( domain, &d);
    allocate_cartesian_matrix( domain, &d);
    d.plane_vof = FALSE;
    if (realdim == 3)
      fill_cartesian_matrix( domain, &d);
    else if (realdim == 2)
      fill_cartesian_matrix_plane( domain, &d);
    else
      g_error("FFT in %i dimension is not implemented yet", realdim);
    p = get_fftw_plan (&d);
    fftw_execute(p); 
    write_spectra ( GFS_OUTPUT (event)->file->fp, &d);

    fftw_destroy_plan(p); 
    g_free(d.domain_ranges);
    g_free(d.out);
    gts_object_destroy (GTS_OBJECT (d.cgd));

#ifdef HAVE_MPI  
    fftw_mpi_cleanup();
#endif

    return TRUE;
  }
  return FALSE;
}

static void output_spectra_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_spectra_class ())->parent_class->read) (o, fp); 
  if (fp->type == GTS_ERROR)
    return;
  GfsOutputSpectra * v = GFS_OUTPUT_SPECTRA (*o);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (v)");
    return;
  }
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(v->v = gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x0", TRUE, &v->pos_min.x},
    {GTS_DOUBLE, "y0", TRUE, &v->pos_min.y},
    {GTS_DOUBLE, "z0", TRUE, &v->pos_min.z},
    {GTS_DOUBLE, "x1", TRUE, &v->pos_max.x},
    {GTS_DOUBLE, "y1", TRUE, &v->pos_max.y},
    {GTS_DOUBLE, "z1", TRUE, &v->pos_max.z},
    {GTS_NONE}
  };
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type == GTS_INT) {
    v->level = atoi (fp->token->str);
    gts_file_next_token (fp);
  }
  else
    v->level = gfs_domain_depth (domain);

}

static void output_spectra_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_spectra_class ())->parent_class->write) (o, fp); 

  GfsOutputSpectra * v = GFS_OUTPUT_SPECTRA (o);
  fprintf (fp, " %s { x0 = %g y0 = %g z0 = %g x1 = %g y1 = %g z1 = %g } %d ",
      v->v->name, v->pos_min.x, v->pos_min.y, v->pos_min.z, v->pos_max.x, v->pos_max.y, 
      v->pos_max.z, v->level);
}

static void output_spectra_init ( GtsObject * o )
{
  GfsOutputSpectra * v = GFS_OUTPUT_SPECTRA (o);
  v->pos_min.x = v->pos_min.y = v->pos_min.z = 0.;
  v->pos_max.x = v->pos_max.y = v->pos_max.z = 0.;
}

static void output_spectra_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = output_spectra_event;
  klass->read =  output_spectra_read;
  klass->write = output_spectra_write;
}

GfsOutputClass * gfs_output_spectra_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsOutputSpectra",
      sizeof (GfsOutputSpectra),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) output_spectra_class_init,
      (GtsObjectInitFunc) output_spectra_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsOutputSpectra} */

/** \beginobject{GfsOutputSpectraInterface} */

static gboolean output_spectra_interface_event (GfsEvent * event,  GfsSimulation * sim)
{

  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_spectra_class ())->parent_class)->event)
      (event, sim)) {

    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsOutputSpectra * v = GFS_OUTPUT_SPECTRA (event);
    fftw_plan p;
  SpectraData d;

  d.u = v->v;
  d.v = v->v;
  d.levelmax = 0;

  d.cgd = gfs_cartesian_grid_new (gfs_cartesian_grid_class ()); 

  /* number of dims of the fft */
   d.cgd->N = 3;
   d.cgd->n = g_malloc ( 3 * sizeof ( guint ) );  

   d.Ndim = 3;

   d.pos_min_global = v->pos_min;
   d.pos_max_global = v->pos_max;
   d.pos_min = d.pos_min_global;
   d.pos_max = d.pos_max_global;
   d.levelmax = v->level;

#ifdef HAVE_MPI  
    fftw_mpi_init();
#endif

   d.dx = ftt_level_size(d.levelmax)*sim->physical_params.L;
   allocate_cartesian_matrix( domain, &d);
   d.plane_vof = TRUE;
   fill_cartesian_matrix_plane( domain, &d);
   p = get_fftw_plan (&d);
   fftw_execute(p); 
   write_spectra ( GFS_OUTPUT (event)->file->fp, &d);
   fftw_destroy_plan(p); 

   g_free(d.domain_ranges);
   g_free ( d.out );

#ifdef HAVE_MPI  
    fftw_mpi_cleanup();
#endif
    gts_object_destroy (GTS_OBJECT (d.cgd));

    return TRUE;
  }
}

static void output_spectra_interface_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = output_spectra_interface_event;
}


GfsOutputClass * gfs_output_spectra_interface_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsOutputSpectraInterface",
      sizeof (GfsOutputSpectraInterface),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) output_spectra_interface_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_spectra_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsOutputSpectraInterface} */

/** \beginobject{GfsOutputEnergySpectra} */

static void write_energy_spectra ( FILE * fp, gdouble * Ek, gdouble deltak, gint nk, gdouble Etot)
{
  gint i;
  fprintf(fp, "# Total energy = %g \n", Etot);
  fputs ("# 1:k 2:Ek \n", fp);
  for (i=1; i < nk; i++) 
    fprintf (fp, "%g %g \n", deltak * sqrt((gdouble)i), Ek[i] );
  
}

gint get_index ( gint i, gint j, gint k, gint np)
{
#if FTT_2D
  return j + ( np / 2 + 1 ) * i ;
#else /* 3D */
  return k + ( np / 2 + 1 ) * ( j + np * i );
#endif /* 3D */
}


static gboolean output_energy_spectra_event (GfsEvent * event,  GfsSimulation * sim)
{

  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_energy_spectra_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsOutputEnergySpectra * v = GFS_OUTPUT_ENERGY_SPECTRA (event);
    fftw_plan p;
    gint i,j,k,np,nk,dim;
    gint index1,index2;
    gint knx, kny;
    gdouble * Ek;

    SpectraData d;
    d.levelmax = 0;

    d.cgd = gfs_cartesian_grid_new (gfs_cartesian_grid_class ()); 

    /* number of dims of the fft */
    d.cgd->N = (gint) FTT_DIMENSION + 1;
    d.cgd->n = g_malloc ( d.cgd->N * sizeof ( guint ) );  

    d.Ndim = 0;
    k=1;
    for (i = 0; i < 3; i++) {
      if (((&(v->pos_max.x))[i] - (&(v->pos_min.x))[i]) != 0) {
        d.Ndim++;
        k++;
      }
    }
    g_assert ( d.Ndim > 0);

    d.pos_min_global = v->pos_min;
    d.pos_max_global = v->pos_max;
    get_deep_level (domain, &d);
    d.levelmax = MIN( d.levelmax, v->level);

    fftw_mpi_init();

    GfsVariable ** u = gfs_domain_velocity (domain);
    d.v = u[0];

    /*fixme: it does not work for geometries with different numbers of points */
    get_domain_limits( domain, &d);
    allocate_cartesian_matrix( domain, &d);
    np = d.dirdata[0].np; 
    nk = d.cgd->N * pow( np / 2 + 1, 2);
    Ek = g_malloc0 ( nk * sizeof ( gdouble ) );

    for (dim=0; dim < FTT_DIMENSION; dim++) {
      d.v = u[dim];
      fill_cartesian_matrix( domain, &d);

      p = get_fftw_plan (&d); 
      fftw_execute(p); 
 
#if FTT_2D  
      for (i=d.local_0_start; i < d.local_0_start + d.local_n0; i++) {
        if ( i < ( np / 2 + 1 ) ) knx = i;
        else knx = ( np - i );
        index1 =  pow(knx,2);
        index2 = get_index(i-d.local_0_start,0,0,np);
        Ek[index1] += 0.5*(pow(d.out[index2][0],2) + pow(d.out[index2][1],2));

        for (j=0; j < ( np / 2 ) + 1; j++) {
          index1 =  pow(knx,2) + pow(j,2);
          index2 = get_index(i-d.local_0_start,j,0,np);
          Ek[index1] += (pow(d.out[index2][0],2) + pow(d.out[index2][1],2));
        }
      }
#else /* 3D */
      for (i=d.local_0_start; i < d.local_0_start + d.local_n0; i++) {
        if ( i < ( np / 2 + 1 ) ) knx = i;
        else knx = ( np - i );
        for (j=0; j < np; j++) {
          if ( j < ( np / 2 + 1 ) ) kny = j;
          else kny = ( np - j );
          index1 =  pow(knx,2) + pow(kny,2);
          index2 = get_index(i-d.local_0_start,j,0,np);
          Ek[index1] += 0.5*(pow(d.out[index2][0],2) + pow(d.out[index2][1],2));
            for (k=1; k < ( np / 2 ) + 1 ; k++) {
              index1 = pow(knx,2) + pow(kny,2) + pow(k,2);
              index2 = get_index(i-d.local_0_start,j,k,np);
              Ek[index1] += (pow(d.out[index2][0],2) + pow(d.out[index2][1],2));
            }
        }
      }
#endif /* 3D */

      fftw_destroy_plan(p);
    }

    g_free(d.domain_ranges);
    g_free ( d.out );
    gts_object_destroy (GTS_OBJECT (d.cgd));

    gdouble Etot = 0.;
    for (j=0; j < nk; j++) { 
      gfs_all_reduce (domain, Ek[j], MPI_DOUBLE, MPI_SUM);
      Etot += Ek[j];
    }

    /* fixme: it only works for cubes */
    if (d.td.rank == 0) {
      gdouble deltak = 2*M_PI / ( d.pos_max_global.x -  d.pos_min_global.x ); 
      write_energy_spectra ( GFS_OUTPUT (event)->file->fp, Ek, deltak, nk, Etot);
    }

    g_free(Ek);

    fftw_mpi_cleanup();
    return TRUE;
  }
  return FALSE;
}

static void output_energy_spectra_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_energy_spectra_class ())->parent_class->read) (o, fp); 
  if (fp->type == GTS_ERROR)
    return;
  GfsOutputEnergySpectra * v = GFS_OUTPUT_ENERGY_SPECTRA (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));

  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x0", TRUE, &v->pos_min.x},
    {GTS_DOUBLE, "y0", TRUE, &v->pos_min.y},
    {GTS_DOUBLE, "z0", TRUE, &v->pos_min.z},
    {GTS_DOUBLE, "x1", TRUE, &v->pos_max.x},
    {GTS_DOUBLE, "y1", TRUE, &v->pos_max.y},
    {GTS_DOUBLE, "z1", TRUE, &v->pos_max.z},
    {GTS_NONE}
  };
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type == GTS_INT) {
    v->level = atoi (fp->token->str);
    gts_file_next_token (fp);
  }
  else
    v->level = gfs_domain_depth (domain);

}

static void output_energy_spectra_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_energy_spectra_class ())->parent_class->write) (o, fp); 

  GfsOutputEnergySpectra * v = GFS_OUTPUT_ENERGY_SPECTRA (o);
  fprintf (fp, " { x0 = %g y0 = %g z0 = %g x1 = %g y1 = %g z1 = %g } %d",
      v->pos_min.x, v->pos_min.y, v->pos_min.z, v->pos_max.x, v->pos_max.y, v->pos_max.z, v->level);
}

static void output_energy_spectra_init ( GtsObject * o )
{
  GfsOutputEnergySpectra * v = GFS_OUTPUT_ENERGY_SPECTRA (o);
  v->pos_min.x = v->pos_min.y = v->pos_min.z = 0.;
  v->pos_max.x = v->pos_max.y = v->pos_max.z = 0.;
}

static void output_energy_spectra_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = output_energy_spectra_event;
  klass->read =  output_energy_spectra_read;
  klass->write = output_energy_spectra_write;
}

GfsOutputClass * gfs_output_energy_spectra_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsOutputEnergySpectra",
      sizeof (GfsOutputEnergySpectra),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc)  output_energy_spectra_class_init,
      (GtsObjectInitFunc) output_energy_spectra_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsOutputEnergySpectra} */

/* Initialize module */

const gchar gfs_module_name[] = "fft";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{ 
  gfs_comm_data_class ();
  gfs_parallel_data_class ();
  gfs_output_spectra_class ();
  gfs_output_spectra_interface_class ();
  gfs_output_energy_spectra_class ();
  return NULL; 
}
