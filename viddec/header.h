
/*!
*************************************************************************************
* \file header.h
* 
* \brief
*    Prototypes for header.c
*************************************************************************************
*/

#ifndef _HEADER_H_
#define _HEADER_H_

int FirstPartOfSliceHeader();
int RestOfSliceHeader();

void dec_ref_pic_marking(Bitstream *currStream);

void decode_poc(struct img_par *img);
int dumppoc(struct img_par *img);

#ifdef DIRECTIONAL_FILTER
void readFilterCoefs(int filterID,Bitstream *currStream);

#ifdef E_DAIF
void readFilterCoefs_EDAIF(int filterID,Bitstream *currStream);
void readAIFInteger(Bitstream *bitstream);
void readAIFOffset(int sub_pos, Bitstream *bitstream);
#endif  // E_DAIF
#ifdef EAIF
void readFilterCoefs_EAIF(int filterID,Bitstream *currStream);
void readAIFSubpel(int sub_pos, Bitstream *bitstream);
#endif
#endif  // DIRECTIONAL_FILTER


#endif

