
/*!
 ************************************************************************
 * \file  memalloc.h
 *
 * \brief
 *    Memory allocation and free helper funtions
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 ************************************************************************
 */

#ifndef _MEMALLOC_H_
#define _MEMALLOC_H_

int  get_mem2D(byte ***array2D, int rows, int columns);
int  get_mem3D(byte ****array2D, int frames, int rows, int columns);

int  get_mem2Dint(int ***array2D, int rows, int columns);
int  get_mem3Dint(int ****array3D, int frames, int rows, int columns);
int  get_mem4Dint(int *****array4D, int idx, int frames, int rows, int columns );
int  get_mem5Dint(int ******array5D, int refs, int blocktype, int rows, int columns, int component);

int  get_mem2Dint64(int64 ***array2D, int rows, int columns);
int  get_mem3Dint64(int64 ****array3D, int frames, int rows, int columns);

int  get_mem2Dshort(short ***array2D, int rows, int columns);
int  get_mem3Dshort(short ****array3D, int frames, int rows, int columns);
int  get_mem4Dshort(short *****array4D, int idx, int frames, int rows, int columns );
int  get_mem5Dshort(short ******array5D, int refs, int blocktype, int rows, int columns, int component);
int  get_mem6Dshort(short *******array6D, int list, int refs, int blocktype, int rows, int columns, int component);

int get_mem2Dpel(imgpel ***array2D, int rows, int columns);
int get_mem3Dpel(imgpel ****array3D, int frames, int rows, int columns);

int get_mem2Ddouble(double ***array2D, int rows, int columns);
int get_mem2Ddb_offset(double ***array2D, int rows, int columns, int offset);
int get_mem2Dint_offset(int ***array2D, int rows, int columns, int offset);

#ifdef ADAPTIVE_QUANTIZATION
int  get_mem6Dint(int *******array6D, int len, int refs, int blocktype, int rows, int columns, int component);
void free_mem6Dint(int ******array6D, int len, int refs, int blocktype, int rows);
#endif

void free_mem2D(byte **array2D);
void free_mem3D(byte ***array2D, int frames);

void free_mem2Dint(int **array2D);
void free_mem3Dint(int ***array3D, int frames);
void free_mem4Dint(int ****array4D, int idx, int frames);
void free_mem5Dint(int *****array5D, int refs, int blocktype, int rows);

void free_mem2Dint64(int64 **array2D);
void free_mem3Dint64(int64 ***array3D64, int frames);

void free_mem2Dshort(short **array2D);
void free_mem3Dshort(short ***array3D, int frames);
void free_mem4Dshort(short ****array4D, int idx, int frames);
void free_mem5Dshort(short *****array5D, int refs, int blocktype, int height);
void free_mem6Dshort(short ******array5D, int list, int refs, int blocktype, int height);

void free_mem2Dpel(imgpel **array2D);
void free_mem3Dpel(imgpel ***array3D, int frames);

void free_mem2Ddouble(double **array2D);
void free_mem2Ddb_offset(double **array2D, int offset);
void free_mem2Dint_offset(int **array2D, int offset);

int init_top_bot_planes(imgpel **imgFrame, int rows, int columns, imgpel ***imgTopField, imgpel ***imgBotField);
void free_top_bot_planes(imgpel **imgTopField, imgpel **imgBotField);


void no_mem_exit(char *where);

#ifdef ADAPTIVE_FILTER
int  get_mem1Ddouble(double **array1D, int rows);
int  get_mem2Ddouble(double ***array2D, int rows, int columns);
int  get_mem3Ddouble(double ****array3D, int frames, int rows, int columns);
void free_mem1Ddouble(double *array1D);
void free_mem2Ddouble(double **array2D);
void free_mem3Ddouble(double ***array3D, int frames);
#endif

#ifdef ADAPTIVE_LOOP_FILTER
int  get_mem1Dint(int **array1D, int rows);
void free_mem1Dint(int *array1D);
#endif

#endif
