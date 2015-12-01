/*!
 ***************************************************************************
 *
 * \file spatial_domain_coding.h
 *
 * \brief
 *    prototypes of spatial domain coding functions
 *
 * \date
 *    18.05.2005
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    Matthias Narroschke  narrosch@tnt.uni-hannover.de
 **************************************************************************/

#ifndef _SPATIAL_DOMAIN_CODING_H_
#define _SPATIAL_DOMAIN_CODING_H_



int  de_quantizer(int value);
int  quantizer_index(int value, float adaptive_f);

void encode_in_spatial_domain    (int block8x8,int64  *cbp_blk, int *cbp, int *coeff_cost);
void encode_in_spatial_domain_8x8(int block8x8,int64  *cbp_blk, int *cbp, int *coeff_cost);
void eliminate_expensive_samples    (int bx, int by, int *SSD_best, int *ctr);
void eliminate_expensive_samples_8x8(int bx, int by, int *SSD_best, int *ctr);

#endif
