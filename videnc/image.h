
/*!
 ************************************************************************
 * \file image.h
 *
 * \brief
 *    headers for image processing
 *
 * \author
 *  Inge Lille-Langoy               <inge.lille-langoy@telenor.com>
 *  Copyright (C) 1999  Telenor Satellite Services, Norway
 ************************************************************************
 */
#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "mbuffer.h"
#include "global.h"

extern StorablePicture *enc_picture;
extern StorablePicture *enc_frame_picture;
extern StorablePicture *enc_frame_picture2;
extern StorablePicture *enc_frame_picture3;
extern StorablePicture *enc_top_picture;
extern StorablePicture *enc_bottom_picture;
#ifdef ADAPTIVE_FILTER
extern StorablePicture *enc_frame_picture_aif;
#endif
int encode_one_frame (void);
void report_frame_statistic(void);
Boolean dummy_slice_too_big(int bits_slice);
void copy_rdopt_data (int field_type);       // For MB level field/frame coding tools
#ifdef ADAPTIVE_FD_SD_CODING
void store_rdopt_data_for_FDSD_coding ();
void restore_rdopt_data_for_FDSD_coding ();
void store_rdopt_data_for_FDSD_coding_interlace ();
void restore_rdopt_data_for_FDSD_coding_interlace ();
#endif

void UnifiedOneForthPix (StorablePicture *s);
#ifdef EIGHTH_PEL
void oneeighthpix (StorablePicture *s);
#endif
#ifdef USE_HP_FILTER
void UnifiedOneForthPixWithHighPrecisionH264Filter_16 (StorablePicture *s);
void UnifiedOneForthPixEnhancedFDIF(StorablePicture *s);
#endif
#ifdef USE_NEW_OFFSET
void SetFrameOffset(); 
void CalcFrameOffset();
void ResetFrameOffset(); 
#endif 
#endif

