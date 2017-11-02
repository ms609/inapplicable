//
//  fitch.c
//  MorPhy2
//
//  Created by mbrazeau on 02/05/2017.
//  Copyright © 2017 brazeaulab. All rights reserved.
//
#include "mpl.h"
#include "morphydefs.h"
#include "morphy.h"
#include "mplerror.h"
#include "fitch.h"
#include "statedata.h"

static inline int mpl_check_up_NA_steps
(MPLstate ndset, MPLstate lactive, MPLstate ractive);

static inline int mpl_check_down_NA_steps
(MPLstate left, MPLstate right, MPLstate lactive, MPLstate ractive);


/**/
int mpl_fitch_downpass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int steps = 0;
    const int* indices  = part->charindices;
    int nchars          = part->ncharsinpart;
    MPLstate* left      = lset->downpass1;
    MPLstate* right     = rset->downpass1;
    MPLstate* n         = nset->downpass1;
    
    unsigned long* weights = part->intwts;
    
    for (i = 0; i < nchars; ++i) {
        j = indices[i];
        
        n[j] = left[j] & right[j];
        
        if (n[j] == 0) {
            n[j] = left[j] | right[j];
            steps += weights[i];
        }
    }
    
    // TODO: rewrite for updated stateset checks.
    
    return steps;
}


int mpl_fitch_uppass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset,
 MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    const int* indices  = part->charindices;
    int nchars      = part->ncharsinpart;
    MPLstate* left  = lset->downpass1;
    MPLstate* right = rset->downpass1;
    MPLstate* npre  = nset->downpass1;
    MPLstate* nfin  = nset->uppass1;
    MPLstate* anc   = ancset->uppass1;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        nfin[j] = anc[j] & npre[j];
        
        if (nfin[j] != anc[j]) {
            
            if (left[j] & right[j]) {
                nfin[j] = (npre[j] | (anc[j] & (left[j] | right[j])));
            }
            else {
                nfin[j] = npre[j] | anc[j];
            }
        }
#ifdef DEBUG
        assert(nfin[j]);
#endif
    }
    
    return 0;
}


int mpl_fitch_local_reopt
(MPLndsets* srcset, MPLndsets* tgt1set, MPLndsets* tgt2set, MPLpartition* part,
 int maxlen, bool domaxlen)
{
   
    int i     = 0;
    int j     = 0;
    int steps = 0;
    const int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    MPLstate* tgt1  = tgt1set->uppass1;
    MPLstate* tgt2  = tgt2set->uppass1;
    MPLstate* src   = srcset->downpass1;
    
    unsigned long* weights = part->intwts;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (!(src[j] & (tgt1[j] | tgt2[j]))) {
            
            steps += weights[i];
            
            if (steps > maxlen && domaxlen == true)
            {
                return steps;
            }
        }
    }
    
    return steps;
}


int mpl_NA_fitch_first_downpass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    const int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    MPLstate* left  = lset->downpass1;
    MPLstate* right = rset->downpass1;
    MPLstate* n     = nset->downpass1;
    MPLstate* nt    = nset->temp_downpass1;
    
    for (i = 0; i < nchars; ++i) {
        j = indices[i];
        
        
        n[j] = (left[j] & right[j]);
        
        if (n[j] == 0) {
            n[j] = (left[j] | right[j]);
            
            if ((left[j] & ISAPPLIC) && (right[j] & ISAPPLIC)) {
                n[j] = n[j] & ISAPPLIC;
            }
        }
        else {
            if (n[j] == NA) {
                if ((left[j] & ISAPPLIC) && (right[j] & ISAPPLIC)) {
                    n[j] = (left[j] | right[j]);
                }
            }
        }
        
        nt[j] = n[j];   // Store a copy for partially reoptimising the subtree
#ifdef DEBUG
        assert(n[j]);
#endif
    }
    
    return 0;
}


int mpl_NA_fitch_first_update_downpass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part)
{
    /*------------------------------------------------------------------------*
     |  This function is for doing a partial downpass when proposing a        |
     |  subtree reinsertion during branchswapping. Its purpose is to          |
     |  (partially) correct any character state sets that are affected by     |
     |  the proposed reinsertion. It is nearly identical to its original-     |
     |  pass counterpart except that it does not overwrite the temp state     |
     |  storage.                                                              |
     *------------------------------------------------------------------------*/
    int        i        = 0;
    int        j        = 0;
    const int* indices  = part->update_NA_indices;
    int        nchars   = part->nNAtoupdate;
    MPLstate*  left     = lset->downpass1;
    MPLstate*  right    = rset->downpass1;
    MPLstate*  n        = nset->downpass1;
    MPLstate*  ntemp    = nset->temp_downpass1;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        n[j] = (left[j] & right[j]);
        
        if (n[j] == 0) {
            n[j] = (left[j] | right[j]);
            
            if ((left[j] & ISAPPLIC) && (right[j] & ISAPPLIC)) {
                n[j] = n[j] & ISAPPLIC;
            }
        }
        else {
            if (n[j] == NA) {
                if ((left[j] & ISAPPLIC) && (right[j] & ISAPPLIC)) {
                    n[j] = (left[j] | right[j]);
                }
            }
        }
        
        if (n[j] != ntemp[j]) {
            nset->updated = true;
        }
        
#ifdef DEBUG
        assert(n[j]);
#endif
    }
    
    return 0;
}


int mpl_NA_fitch_first_uppass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset,
 MPLpartition* part)
{
    int         i       = 0;
    int         j       = 0;
    const int*  indices = part->charindices;
    int         nchars  = part->ncharsinpart;
    MPLstate*   left    = lset->downpass1;
    MPLstate*   right   = rset->downpass1;
    MPLstate*   npre    = nset->downpass1;
    MPLstate*   nifin   = nset->uppass1;
    MPLstate*   anc     = ancset->uppass1;
    MPLstate*   nfint   = nset->temp_uppass1;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (npre[j] & NA) {
            if (npre[j] & ISAPPLIC) {
                if (anc[j] == NA) {
                    nifin[j] = NA;
                }
                else {
                    nifin[j] = npre[j] & ISAPPLIC;
                }
            }
            else {
                if (anc[j] == NA) {
                    nifin[j] = NA;
                }
                else {
                    if ((left[j] | right[j]) & ISAPPLIC
                        && (left[j] | right[j]) != MISSING) { // TODO: This possibly isn't quite safe or right.
                        nifin[j] = ((left[j] | right[j]) & ISAPPLIC);
                    }
                    else {
                        nifin[j] = NA;
                    }
                }
            }
        }
        else {
            nifin[j] = npre[j];
        }
        
        // Store the set for restoration during tree searches.
        nfint[j] = nifin[j];
        
#ifdef DEBUG
        assert(nifin[j]);
#endif
    }
    
    return 0;
}


int mpl_NA_fitch_first_update_uppass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset,
 MPLpartition* part)
{
    /*------------------------------------------------------------------------*
     |  This function is for doing a partial uppsass when proposing a subtree |
     |  reinsertion during branchswapping. Its purpose is to (partially)      |
     |  correct any character state sets that are affected by the proposed    |
     |  reinsertion.                                                          |
     *------------------------------------------------------------------------*/
    int         i       = 0;
    int         j       = 0;
    const int*  indices = part->update_NA_indices;
    int         nchars  = part->nNAtoupdate;
    MPLstate*   left    = lset->downpass1;
    MPLstate*   right   = rset->downpass1;
    MPLstate*   npre    = nset->downpass1;
    MPLstate*   nifin   = nset->uppass1;
    MPLstate*   anc     = ancset->uppass1;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (npre[j] & NA) {
            if (npre[j] & ISAPPLIC) {
                if (anc[j] == NA) {
                    nifin[j] = NA;
                }
                else {
                    nifin[j] = npre[j] & ISAPPLIC;
                }
            }
            else {
                if (anc[j] == NA) {
                    nifin[j] = NA;
                }
                else {
                    if ((left[j] | right[j]) & ISAPPLIC
                        && (left[j] | right[j]) != MISSING) { // TODO: This possibly isn't quite safe or right.
                        nifin[j] = ((left[j] | right[j]) & ISAPPLIC);
                    }
                    else {
                        nifin[j] = NA;
                    }
                }
            }
        }
        else {
            nifin[j] = npre[j];
        }
        
#ifdef DEBUG
        assert(nifin[j]);
#endif
    }
    
    return 0;
}


int mpl_NA_fitch_second_downpass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part)
{
    int             i       = 0;
    int             j       = 0;
    int             steps   = 0;
    const int*      indices = part->charindices;
    int             nchars  = part->ncharsinpart;
    MPLstate*       left    = lset->downpass2;
    MPLstate*       right   = rset->downpass2;
    MPLstate*       nifin   = nset->uppass1;
    MPLstate*       npre    = nset->downpass2;
    MPLstate*       npret   = nset->temp_downpass2;
    MPLstate*       stacts  = nset->subtree_actives;
    MPLstate*       tstatcs = nset->temp_subtr_actives;
    MPLstate*       lacts   = lset->subtree_actives;
    MPLstate*       racts   = rset->subtree_actives;
    MPLstate        temp    = 0;
    unsigned long*  weights = part->intwts;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (nifin[j] & ISAPPLIC) {
            if ((temp = (left[j] & right[j]))) {
                if (temp & ISAPPLIC) {
                    npre[j] = temp & ISAPPLIC;
                } else {
                    npre[j] = temp;
                }
            }
            else {
                npre[j] = (left[j] | right[j]) & ISAPPLIC;
                
                if (left[j] & ISAPPLIC && right[j] & ISAPPLIC) {
                    steps += weights[i];
                } else if (lacts[j] && racts[j]) {
                    steps += weights[i];
                }
            }
        }
        else {
            npre[j] = nifin[j];
        }
        
        npret[j] = npre[j]; // Storage for temporary updates.
        stacts[j] = (lacts[j] | racts[j]) & ISAPPLIC;
        tstatcs[j] = stacts[j]; // Storage for temporary updates.
    
#ifdef DEBUG
        assert(npre[j]);
#endif
    }
    
    return steps;
}


static inline int mpl_check_down_NA_steps
(MPLstate left, MPLstate right, MPLstate lactive, MPLstate ractive)
{
    int steps = 0;

#ifdef DEBUG    
    assert(left && right);
#endif
    
    if (!(left & right)) {
        if (left & ISAPPLIC && right & ISAPPLIC) {
            ++steps;
        }
        else if (lactive && ractive) {
            ++steps;
        }
    }
    
    return steps;
}


int mpl_NA_fitch_second_update_downpass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part)
{
    /*------------------------------------------------------------------------*
     |  This function is for doing a partial downpass when proposing a        |
     |  subtree reinsertion during branchswapping. Its purpose is to          |
     |  (partially) correct any character state sets that are affected by     |
     |  the proposed reinsertion. It is nearly identical to its original-     |
     |  pass counterpart except that it does not overwrite the temp state     |
     |  storage.                                                              |
     *------------------------------------------------------------------------*/
    int             i           = 0;
    int             j           = 0;
    int             steps       = 0;
    int             step_recall = 0;
    const int*      indices     = part->update_NA_indices;
    int             nchars      = part->nNAtoupdate;
    MPLstate*       left        = lset->downpass2;
    MPLstate*       right       = rset->downpass2;
    const MPLstate* tleft       = lset->temp_downpass2;
    const MPLstate* tright      = rset->temp_downpass2;
    MPLstate*       nifin       = nset->uppass1;
    MPLstate*       npre        = nset->downpass2;
    const MPLstate* npret       = nset->temp_downpass2;
    MPLstate*       stacts      = nset->subtree_actives;
    MPLstate*       lacts       = lset->subtree_actives;
    MPLstate*       racts       = rset->subtree_actives;
    MPLstate*       tlacts      = lset->temp_subtr_actives;
    MPLstate*       tracts      = rset->temp_subtr_actives;
    MPLstate        temp        = 0;
    unsigned long*  weights     = part->intwts;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (nifin[j] & ISAPPLIC) {
            if ((temp = (left[j] & right[j]))) {
                if (temp & ISAPPLIC) {
                    npre[j] = temp & ISAPPLIC;
                } else {
                    npre[j] = temp;
                }
            }
            else {
                npre[j] = (left[j] | right[j]) & ISAPPLIC;
                
                if (left[j] & ISAPPLIC && right[j] & ISAPPLIC) {
                    steps += weights[i];
                } else if (lacts[j] && racts[j]) {
                    steps += weights[i];
                }
            }
        }
        else {
            npre[j] = nifin[j];
        }
        
        stacts[j] = (lacts[j] | racts[j]) & ISAPPLIC;
        
        /* Count whether any steps need to be taken back at this node */
        
        /* Flag as updated if current set is different from previous */
        if (npre[j] != npret[j]) {
            nset->updated = true;
//            if (nset->temp_uppass1[j] & ISAPPLIC) {
//                int rec = 0;
//                rec = mpl_check_down_NA_steps(tleft[j], tright[j], tlacts[j], tracts[j]);
//                step_recall += (weights[j] * rec);
//            }
        }
#ifdef DEBUG
        assert(npre[j]);
#endif
    }
    
    nset->steps_to_recall = step_recall;
    
    return steps;
}


int mpl_NA_fitch_second_uppass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset,
 MPLpartition* part)
{
    int             i       = 0;
    int             j       = 0;
    int             steps   = 0;
    const int*      indices = part->charindices;
    int             nchars  = part->ncharsinpart;
    MPLstate*       left    = lset->downpass2;
    MPLstate*       right   = rset->downpass2;
    MPLstate*       npre    = nset->downpass2;
    MPLstate*       nfin    = nset->uppass2;
    MPLstate*       nfint   = nset->temp_uppass2;
    MPLstate*       anc     = ancset->uppass2;
    MPLstate*       lacts   = lset->subtree_actives;
    MPLstate*       racts   = rset->subtree_actives;
    unsigned long*  weights = part->intwts;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (npre[j] & ISAPPLIC) {
            if (anc[j] & ISAPPLIC) {
                if ((anc[j] & npre[j]) == anc[j]) {
                    nfin[j] = anc[j] & npre[j];
                } else {
                    if (left[j] & right[j]) {
                        nfin[j] = (npre[j] | (anc[j] & (left[j] | right[j])));
                    }
                    else {
                        if ((left[j] | right[j]) & NA) {
                            if ((left[j] | right[j]) & anc[j]) {
                                nfin[j] = anc[j];
                            } else {
                                nfin[j] = (left[j] | right[j] | anc[j]) & ISAPPLIC;
                            }
                        } else {
                            nfin[j] = npre[j] | anc[j];
//                            if ((anc[j] & nfin[j]) == anc[j]) {
//                                nfin[j] = anc[j] & nfin[j];
//                            }
                        }
                    }
                }
            }
            else {
                nfin[j] = npre[j];
            }
        }
        else {
            nfin[j] = npre[j];
            
            if (lacts[j] && racts[j]) {
                steps += weights[i];
            }
        }
        
        nfint[j] = nfin[j]; // Storage of states for undoing temp updates
#ifdef DEBUG
        assert(nfin[j]);
#endif
    }
    
    return steps;
}


static inline int mpl_check_up_NA_steps
(MPLstate ndset, MPLstate lactive, MPLstate ractive)
{
    int steps = 0;
    
    if (ndset == NA) {
        if (lactive != 0 && ractive != 0) {
            ++steps;
        }
    }
    
    return steps;
}

int mpl_NA_fitch_second_update_uppass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset,
 MPLpartition* part)
{
    int             i           = 0;
    int             j           = 0;
    int             steps       = 0;
    int             step_recall = 0;
    const int*      indices     = part->update_NA_indices;
    int             nchars      = part->nNAtoupdate;
    MPLstate*       left        = lset->downpass2;
    MPLstate*       right       = rset->downpass2;
    MPLstate*       npre        = nset->downpass2;
    MPLstate*       nfin        = nset->uppass2;
    MPLstate*       nfint       = nset->temp_uppass2;
    MPLstate*       anc         = ancset->uppass2;
    MPLstate*       lacts       = lset->subtree_actives;
    MPLstate*       racts       = rset->subtree_actives;
    MPLstate*       tlacts      = lset->temp_subtr_actives;
    MPLstate*       tracts      = lset->temp_subtr_actives;
    unsigned long*  weights     = part->intwts;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (npre[j] & ISAPPLIC) {
            if (anc[j] & ISAPPLIC) {
                if ((anc[j] & npre[j]) == anc[j]) {
                    nfin[j] = anc[j] & npre[j];
                } else {
                    if (left[j] & right[j]) {
                        nfin[j] = (npre[j] | (anc[j] & (left[j] | right[j])));
                    }
                    else {
                        if ((left[j] | right[j]) & NA) {
                            if ((left[j] | right[j]) & anc[j]) {
                                nfin[j] = anc[j];
                            } else {
                                nfin[j] = (left[j] | right[j] | anc[j]) & ISAPPLIC;
                            }
                        } else {
                            nfin[j] = npre[j] | anc[j];
//                            if ((anc[j] & nfin[j]) == anc[j]) {
//                                nfin[j] = anc[j] & nfin[j];
//                            }
                        }
                    }
                }
            }
            else {
                nfin[j] = npre[j];
            }
        }
        else {
            nfin[j] = npre[j];
            
            if (lacts[j] && racts[j]) {
                steps += weights[i];
            }
        }
        
        if (nfint[j] != nfin[j]) {
            nset->updated = true;
            int rec = 0;
            rec = mpl_check_up_NA_steps(nfint[j], tlacts[j], tracts[j]);
            step_recall += (rec * weights[i]);
        }
        
#ifdef DEBUG
        assert(nfin[j]);
#endif
    }
    
    nset->steps_to_recall = step_recall;
    
    return steps;
}


int mpl_fitch_NA_local_reopt
(MPLndsets* srcset, MPLndsets* tgt1set, MPLndsets* tgt2set, MPLpartition* part,
 int maxlen, bool domaxlen)
{
    
    part->ntoupdate = 0; // V. important: resets the record of characters needing updates
    
    int i           = 0;
    int j           = 0;
    int need_update = 0;
    int steps       = 0;
    const int* indices  = part->charindices;
    int nchars          = part->ncharsinpart;
    MPLstate* tgt1d1    = tgt1set->downpass1;
//    MPLstate* tgt1u1    = tgt1set->uppass1;
    MPLstate* tgt2d1    = tgt2set->downpass1;
//    MPLstate* tgt2u1    = tgt2set->uppass1;
    MPLstate* tgt1f     = tgt1set->uppass2;
    MPLstate* tgt2f     = tgt2set->downpass2;
    MPLstate* src       = srcset->downpass1; // TODO: Verify this.
    
    unsigned long* weights = part->intwts;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (!(src[j] & (tgt1f[j] | tgt2f[j]))) {
            
            if (src[j] & ISAPPLIC) {
                if ((tgt1f[j] | tgt2f[j]) & ISAPPLIC) {
                    steps += weights[i];
                }
                else {
                    
                    /* NOTE: This will be written simply at first, but there are
                     * possible additional checks on tgt preliminary sets that 
                     * could reduce the number of characters that need to be 
                     * updated. */
                    
                    part->update_NA_indices[need_update] = j;
                    ++need_update;
                }
            }
            else {
                if (src[j] & (tgt1d1[j] | tgt2d1[j])) {
                    part->update_NA_indices[need_update] = j;
                    ++need_update;
                }
            }
        }
    }
    
    part->nNAtoupdate = need_update;
    
    return steps;
}


int mpl_fitch_tip_update(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    // TODO: Check these!!!!!!!
    MPLstate* tprelim = tset->downpass1;
    MPLstate* tfinal  = tset->uppass1;
    MPLstate* ttfinal = tset->temp_uppass1;
    MPLstate* astates = ancset->uppass1;
    
    for (i = 0; i < nchars; ++i) {
        j = indices[i];
        if (tprelim[j] & astates[j]) {
            tfinal[j] = tprelim[j] & astates[j];
        }
        else {
            tfinal[j] = tprelim[j];
        }
        ttfinal[j] = tfinal[j];
#ifdef DEBUG
        assert(tfinal[j]);
#endif
    }
    return 0;
}

int mpl_fitch_one_branch
(MPLndsets* tipanc, MPLndsets* node, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices     = part->charindices;
    int nchars       = part->ncharsinpart;
    MPLstate* tipset = tipanc->downpass1;
    MPLstate* tipfin = tipanc->uppass1;
    MPLstate* ndset  = node->downpass1;
    MPLstate temp    = 0;
    unsigned long* weights = part->intwts;
    int length = 0;
    
    for (i = 0; i < nchars; ++i) {
        j = indices[i];
        
        temp = tipset[j] & ndset[j];

        if (temp == 0) {
            tipfin[j] = tipset[j];
            length += weights[i];
            node->uppass1[j] = ndset[j];
        }
        else {
            tipfin[j] = temp;
            node->uppass1[j] = temp;
        }
    }
    
    return length;
}

int mpl_fitch_NA_one_branch
(MPLndsets* tipanc, MPLndsets* node, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices     = part->charindices;
    int nchars       = part->ncharsinpart;
    MPLstate* tipset = tipanc->downpass1;
    MPLstate* tipifin = tipanc->uppass1;
    MPLstate* tipfin = tipanc->uppass2;
    MPLstate* ndset  = node->downpass1;
    MPLstate* ndacts = node->subtree_actives;
    MPLstate  temp   = 0;
    unsigned long* weights = part->intwts;
    int length = 0;
    
    for (i = 0; i < nchars; ++i) {
        j = indices[i];
        
        temp = tipset[j] & ndset[j];
        
        if (temp == 0) {
            if (tipset[j] & ISAPPLIC) {
                if (ndset[j] & ISAPPLIC) {
                    length += weights[i];
                }
                else {
                    if (ndacts[j]) {
                        length += weights[i];
                    }
                }
            }
            
            tipifin[j]        = tipset[j];
            //node->uppass2[j] = ndset[j];
        }
        else {
            tipifin[j]        = temp;
            //node->uppass1[j] = temp;
            //node->uppass2[j] = temp;
        }
    }
    
    return length;
}

int mpl_fitch_NA_tip_update
(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices        = part->charindices;
    int nchars          = part->ncharsinpart;
    
    MPLstate* tpass1    = tset->downpass1;
    MPLstate* tpass2    = tset->uppass1;
    MPLstate* tpass3    = tset->downpass2;
    MPLstate* ttpass1   = tset->temp_downpass1;
    MPLstate* ttpass2   = tset->temp_uppass1;
    MPLstate* ttpass3   = tset->temp_downpass2;
    MPLstate* astates   = ancset->uppass1;
    MPLstate* stacts    = tset->subtree_actives;
    MPLstate* tstatcs   = tset->temp_subtr_actives;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (tpass1[j] & astates[j]) {
            stacts[j] = (tpass1[j] & astates[j] & ISAPPLIC);
        }
        else {
            stacts[j] |= tpass1[j] & ISAPPLIC;
        }

        tpass2[j] = tpass1[j];
        
        if (tpass2[j] & astates[j]) {
            if (astates[j] & ISAPPLIC) {
                tpass2[j] &= ISAPPLIC;
            }
        }
        
        tpass3[j]  = tpass2[j];
        
        // Store the temp sets for restoring after temporary updates
        ttpass1[j] = tpass1[j];
        ttpass2[j] = tpass2[j];
        ttpass3[j] = tpass3[j];
        tstatcs[j] = stacts[j];
#ifdef DEBUG   
        assert(tpass3[j]);
        assert(tpass2[j]);
#endif
    }
    
    return 0;
}

int mpl_fitch_NA_tip_finalize
(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    MPLstate* tpass1    = tset->downpass1;
    MPLstate* tfinal    = tset->uppass2;
    MPLstate* ttfinal   = tset->temp_uppass2;
    MPLstate* astates   = ancset->uppass2;
    MPLstate* stacts    = tset->subtree_actives;
    MPLstate* tstacts   = tset->temp_subtr_actives;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (tpass1[j] & astates[j]) {
            tfinal[j] = tpass1[j] & astates[j];
        }
        else {
            tfinal[j] = tpass1[j];
        }
        
        stacts[j] = tfinal[j] & ISAPPLIC;
        
        // Store the temp buffers:
        ttfinal[j] = tfinal[j];
        tstacts[j] = stacts[j];
#ifdef DEBUG
        assert(tfinal[j]);
#endif
    }
    
    return 0;
}
