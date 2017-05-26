//
//  morphy.h
//  MorPhy2
//
//  Created by mbrazeau on 23/04/2017.
//  Copyright © 2017 brazeaulab. All rights reserved.
//

#ifndef morphy_h
#define morphy_h

#ifdef DEBUG
#include <stdio.h>
#define dbg_printf(...) printf(__VA_ARGS__)
#else
#define dbg_printf(...)
#endif

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <float.h>
#include <math.h>
//#include <glib.h>

//typedef struct mplarray_s {
//    int nelems;
//    int maxelems;
//    size_t elemsize;
//    void **data;
//} MPLarray;

/* Function prototypes */

Morphyp         mpl_new_Morphy_t(void);
void*           mpl_get_from_matrix(const int row, const int col, const int ncol, const size_t size, const void* data);
int             mpl_set_numtaxa(const int ntax, Morphyp m);
int             mpl_set_num_charac(const int ncharac, Morphyp m);
int             mpl_get_gaphandl(Morphyp handl);
int             mpl_check_data_loaded(Morphyp m);
char            mpl_get_gap_symbol(Morphyp handl);
bool            mpl_almost_equal(double a, double b);
bool            mpl_isreal(const double n);
void            mpl_set_new_weight_public(const double wt, const int char_id, Morphyp handl);
void            mpl_scale_all_intweights(Morphyp handl);
MPLchtype*      mpl_get_charac_types(Morphyp handl);
int             mpl_assign_partition_fxns(MPLpartition* part);
int             mpl_fetch_parsim_fxn_setter (void(**pars_assign)(MPLpartition*), MPLchtype chtype);
int             mpl_extend_intarray(int** array, size_t size);
int             mpl_part_push_index(int newint, MPLpartition* part);
int             mpl_part_remove_index(int index, MPLpartition* part);
int             mpl_delete_partition(MPLpartition* part);
MPLpartition*   mpl_new_partition(const MPLchtype chtype, const bool hasNA);
int             mpl_count_gaps_in_columns(Morphyp handl);
int             mpl_setup_partitions(Morphyp handle);
int             mpl_get_numparts(Morphyp handl);
MPLndsets*      mpl_alloc_stateset(int numchars);
void            mpl_free_stateset(MPLndsets* statesets);
int             mpl_delete_all_partitions(Morphyp handl);
int             mpl_setup_statesets(Morphyp handl);
int             mpl_destroy_statesets(Morphyp handl);
int             mpl_copy_data_into_tips(Morphyp handl);
int             mpl_assign_intwts_to_partitions(Morphyp handl);
int             mpl_update_root(MPLndsets* lower, MPLndsets* upper, MPLpartition* part);
int             mpl_update_NA_root(MPLndsets* lower, MPLndsets* upper, MPLpartition* part);

//MPLarray*   mpl_new_array(size_t elemsize);
//void        mpl_destroy_array(MPLarray* arr);
#endif /* morphy_h */
