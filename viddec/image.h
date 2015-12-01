
/*!
 ************************************************************************
 * \file image.h
 *
 * \brief
 *    prototypes for image.c
 *
 ************************************************************************
 */

#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "mbuffer.h"

#include "global.h"

extern StorablePicture *dec_picture;

void find_snr(struct snr_par *snr, StorablePicture *p, int p_ref);
void get_block(int ref_frame, StorablePicture **list, int x_pos, int y_pos, struct img_par *img, int block[BLOCK_SIZE][BLOCK_SIZE]);
int  picture_order(struct img_par *img);


#ifdef EIGHTH_PEL
void average_block(int block[BLOCK_SIZE][BLOCK_SIZE], int block2[BLOCK_SIZE][BLOCK_SIZE]); 
void average_block4(int block[BLOCK_SIZE][BLOCK_SIZE], int block2[BLOCK_SIZE][BLOCK_SIZE], int block3[BLOCK_SIZE][BLOCK_SIZE], int block4[BLOCK_SIZE][BLOCK_SIZE]); 
void average_block_2(int block[MB_BLOCK_SIZE][MB_BLOCK_SIZE], int block2[MB_BLOCK_SIZE][MB_BLOCK_SIZE]); 
void average_block4_2(int block[MB_BLOCK_SIZE][MB_BLOCK_SIZE], int block2[MB_BLOCK_SIZE][MB_BLOCK_SIZE], int block3[MB_BLOCK_SIZE][MB_BLOCK_SIZE], int block4[MB_BLOCK_SIZE][MB_BLOCK_SIZE]);
#endif
#endif

