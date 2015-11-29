
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

int SliceHeader();
int Partition_BC_Header();

int  writeERPS(SyntaxElement *sym, DataPartition *partition);
// int  SequenceHeader(FILE *outf);
void write_terminating_bit (short);
int nBitsAIF(void);

#ifdef DIRECTIONAL_FILTER
int sendCoefsAIF(Bitstream *bitstream);

#endif  // DIRECTIONAL_FILTER

#ifdef E_DAIF
int sendCoefs_EDAIF(Bitstream *bitstream);
int estimateCostOfCoefs_EDAIF(int sub_pos);
int sendAIFInteger(Bitstream *bitstream);
int sendAIFOffset(int sub_pos, Bitstream *bitstream);       // Used for both full and subpel offsets
#endif  // E_DAIF

#ifdef EAIF
int sendAIFSubpel(int sub_pos, Bitstream *bitstream);
int sendCoefs_EAIF(Bitstream *bitstream);
int estimateCostOfCoefs_EAIF(int sub_pos);
#endif

#endif

