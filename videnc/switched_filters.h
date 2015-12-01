//
//
//
#ifndef SWITCHED_FILTERS_H
#define SWITCHED_FILTERS_H

#define NUM_SIFO                        3             // Number of filters
#define NORM_FACTOR                     0.25          // 0.5

#define NO_SAMPLES                      {  1,  2, 16, 64} // noSamples
#define BL_SIZE                         16            // Block size
#define NO_DC_VAL                       64            // Number of DC levels
#define USE_DCMINMAX
#define THRESHOLD_FRAME_OFFSETS

#define HPF_SIFO                        3             // input->UseHPFilter
#define HPF_SIFO_FPO                    4             // input->UseHPFilter

// Defines for img->filterParam
#define SIFO_FIRST_PASS                 0             // Basic first pass with sequence filters and zero offsets
#define SIFO_SEQ_FILTER                 1             // Use sequence filters and no offsets
#define SIFO_FRAME_FILTER_WITH_OFFSET   2             // Use frame filters with offsets
#define SIFO_FIRST_PASS_FPO             3             // Use sequence filters and first pass offsets 
#define SIFO_INTRA                      4             // P slice converted to I
#define SIFO_FRAME_FILTER               5             // Use frame filters with no offsets

#define SIFO_FILTERS_ONLY               0
#define SIFO_USE_OFFSETS                1
#define NO_BEST_FILTER                  0

void SwapFilteredFrames(void);                        // slice.c
void UnifiedOneForthPixFiltSel(int filterParam);      // slice.c
void rdPictureCodingFilterP(void);                    // image.c
void rdPictureCodingFilterB(void);                    // image.c
void rdPictureCodingFilterP_IPass(void);              // image.c

void initSIFOFilters();                               // mbuffer.c
int Clip3Fun(int low, int high, int val);             // epzs.c
void setFirstPassSubpelOffset(int list);              // slice.c

int ComputeFiltersAndOffsets(void);                   // image.c
void UpdateSequenceFilters_B(void);                   // image.c
void UpdateSequenceFilters_P(void);                   // image.c

void resetFirstPassSubpelOffset(int list);            // slice.c
int getNonzero(int list);                             // header.c
#endif  // SWITCHED_FILTERS_H

