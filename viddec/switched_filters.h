//
//
//
#ifndef SWITCHED_FILTERS_H
#define SWITCHED_FILTERS_H

#include "mbuffer.h"

#define NUM_SIFO                        3             // Number of filters

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

void readSIFOHeader();
int getblock2DFilt_quarter_pel(int ref_idx, StorablePicture **list, int x_pos, int y_pos,
                               int mvx, int mvy, int img_width, int img_height, 
                               int block[BLOCK_SIZE][BLOCK_SIZE]);
int getblock2DFilt(int ref_idx, StorablePicture **list, int x_pos, int y_pos,
                   int mvx, int mvy, int img_width, int img_height, int block[BLOCK_SIZE][BLOCK_SIZE]);


#endif  // SWITCHED_FILTERS_H
