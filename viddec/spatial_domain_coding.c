/*!
 ***************************************************************************
 * \file spatial_domain_coding.c
 *
 * \brief
 *    coding of 4x4 blocks in spatial domain
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    Matthias Narroschke  narrosch@tnt.uni-hannover.de
 * \date
 *    18. May 2005
 **************************************************************************
 */

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "global.h"
#include "spatial_domain_coding.h"

#ifdef ADAPTIVE_FD_SD_CODING



static const int rep_values_qp00_04[256]={0,  1,  2,   3,   4 ,  5 ,  6 ,  7  ,   8 , 9 , 10 , 11 , 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 , 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 , 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 , 36 , 37 , 38 , 39 , 40 , 41 , 42 , 43 , 44 , 45 , 46 , 47 , 48 , 49 , 50 , 51 , 52 , 53 , 54 , 55 , 56 , 57 , 58 , 59 , 60 , 61 , 62 , 63 , 64 , 65 , 66 , 67 , 68 , 69 , 70 , 71 , 72 , 73 , 74 , 75 , 76 , 77 , 78 , 79 , 80 , 81 , 82 , 83 , 84 , 85 , 86 , 87 , 88 , 89 , 90 , 91 , 92 , 93 , 94 , 95 , 96 , 97 , 98 , 99 , 100 , 101 , 102 , 103 , 104 , 105 , 106 , 107 , 108 , 109 , 110 , 111 , 112 , 113 , 114 , 115 , 116 , 117 , 118 , 119 , 120 , 121 , 122 , 123 , 124 , 125 , 126 , 127 , 128 , 129 , 130 , 131 , 132 , 133 , 134 , 135 , 136 , 137 , 138 , 139 , 140 , 141 , 142 , 143 , 144 , 145 , 146 , 147 , 148 , 149 , 150 , 151 , 152 , 153 , 154 , 155 , 156 , 157 , 158 , 159 , 160 , 161 , 162 , 163 , 164 , 165 , 166 , 167 , 168 , 169 , 170 , 171 , 172 , 173 , 174 , 175 , 176 , 177 , 178 , 179 , 180 , 181 , 182 , 183 , 184 , 185 , 186 , 187 , 188 , 189 , 190 , 191 , 192 , 193 , 194 , 195 , 196 , 197 , 198 , 199 , 200 , 201 , 202 , 203 , 204 , 205 , 206 , 207 , 208 , 209 , 210 , 211 , 212 , 213 , 214 , 215 , 216 , 217 , 218 , 219 , 220 , 221 , 222 , 223 , 224 , 225 , 226 , 227 , 228 , 229 , 230 , 231 , 232 , 233 , 234 , 235 , 236 , 237 , 238 , 239 , 240 , 241 , 242 , 243 , 244 , 245 , 246 , 247 , 248 , 249 , 250 , 251 , 252 , 253 , 254 , 255};
static const int rep_values_qp05_09[128]={0,  2,  4,   6,   8 , 10 , 12 , 14  ,  16 , 18 , 20 , 22 , 24 , 26 , 28 , 30 , 32 , 34 , 36 , 38 , 40 , 42 , 44 , 46 , 48 , 50 , 52 , 54 , 56 , 58 , 60 , 62 , 64 , 66 , 68 , 70 , 72 , 74 , 76 , 78, 80 , 82 , 84 , 86 , 88 , 90 , 92 , 94 , 96 , 98 , 100 , 102 , 104 , 106 , 108 , 110 , 112 , 114 , 116 , 118 , 120 , 122 , 124 , 126 , 128 , 130 , 132 , 134 , 136 , 138 , 140 , 142 , 144 , 146 , 148 , 150 , 152 , 154 , 156 , 158 , 160 , 162 , 164 , 166 , 168 , 170 , 172 , 174 , 176 , 178 , 180 , 182 , 184 , 186 , 188 , 190 , 192 , 194 , 196 , 198 , 200 , 202 , 204 , 206 , 208 , 210 , 212 , 214 , 216 , 218 , 220 , 222 , 224 , 226 , 228 , 230 , 232 , 234 , 236 , 238 , 240 , 242 , 244 , 246 , 248 , 250 , 252 , 254};
static const int rep_values_qp10[87]={0,  3,  6,   9,  12 , 15 , 18 , 21  ,  24 , 27 , 30 , 33 , 36 , 39 , 42 , 45 , 48 , 51 , 54 , 57 , 60 , 63 , 66 , 69 , 72 , 75 , 78 , 81 , 84 , 87 , 90 , 93 , 96 , 99 , 102 , 105 , 108 , 111 , 114 , 117 , 120 , 123 , 126 , 129 , 132 , 135 , 138 , 141 , 144 , 147 , 150 , 153 , 156 , 159 , 162 , 165 , 168 , 171 , 174 , 177 , 180 , 183 , 186 , 189 , 192 , 195 , 198 , 201 , 204 , 207 , 210 , 213 , 216 , 219 , 222 , 225 , 228 , 231 , 234 , 237 , 240 , 243 , 246 , 249 , 252 , 255};
static const int rep_values_qp11[65]={0,  3,  7,  11,  15 , 19 , 23 , 27  ,  31  ,  35 , 39 , 43 , 47 , 51 , 55 , 59 , 63 , 67 , 71 , 75 , 79 , 83 , 87 , 91 , 95 , 99 , 103 , 107 , 111 , 115 , 119 , 123 , 127 , 131 , 135 , 139 , 143 , 147 , 151 , 155 , 159 , 163 , 167 , 171 , 175 , 179 , 183 , 187 , 191 , 195 , 199 , 203 , 207 , 211 , 215 , 219 , 223 , 227 , 231 , 235 , 239 , 243 , 247 , 251 , 255};
static const int rep_values_qp12[65]={0,  4,  8,  12,  16 , 20 , 24 , 28  ,  31 , 35 , 39 , 43 , 47 , 51 , 55 , 59 , 63 , 67 , 71 , 75 , 79 , 83 , 87 , 91 , 95 , 99 , 103 , 107 , 111 , 115 , 119 , 123 , 127 , 131 , 135 , 139 , 143 , 147 , 151 , 155 , 159 , 163 , 167 , 171 , 175 , 179 , 183 , 187 , 191 , 195 , 199 , 203 , 207 , 211 , 215 , 219 , 223 , 227 , 231 , 235 , 239 , 243 , 247 , 251 , 255};
static const int rep_values_qp13[44]={0,  4,  8,  12,  18 , 24 , 30 , 36  ,  42 , 48 , 54 , 60 , 66 , 72 , 78 , 84 , 90 , 96 , 102 , 108 , 114 , 120 , 126 , 132 , 138 , 144 , 150 , 156 , 162 , 168 , 174 , 180 , 186 , 192 , 198 , 204 , 210 , 216 , 222 , 228 , 234 , 240 , 246 , 252};
static const int rep_values_qp14[38]={0,  4,  8,  13,  19 , 25 , 32 , 39  ,  46 , 53 , 60 , 67 , 74 , 81 , 88 , 95 , 102 , 109 , 116 , 123 , 130 , 137 , 144 , 151 , 158 , 165 , 172 , 179 , 186 , 193 , 200 , 207 , 214 , 221 , 228 , 235 , 242 , 249};
static const int rep_values_qp15[34]={0,  4,  8,  14,  20 , 26 , 34 , 42  ,  50  , 58 , 66 , 74 , 82 , 90 , 98 , 106 , 114 , 122 , 130 , 138 , 146 , 154 , 162 , 170 , 178 , 186 , 194 , 202 , 210 , 218 , 226 , 234 , 242 , 250};
static const int rep_values_qp16[34]={0,  5, 10,  17,  25 , 32 , 39 , 46  ,  53 , 61 , 69 , 77 , 85 , 93 , 101 , 109 , 117 , 125 , 133 , 141 , 149 , 157 , 165 , 173 , 181 , 189 , 197 , 205 , 213 , 221 , 229 , 237 , 245 , 253};
static const int rep_values_qp17[27]={0,  6, 11,  19,  28 , 36 , 46 , 56  ,  66 , 76 , 86 , 96 , 106 , 116 , 126 , 136 , 146 , 156 , 166 , 176 , 186 , 196 , 206 , 216 , 226 , 236 , 246};
static const int rep_values_qp18[25]={0,  6, 12,  21,  30 , 39 , 49 , 60  ,  71 , 82 , 93 , 104 , 115 , 126 , 137 , 148 , 159 , 170 , 181 , 192 , 203 , 214 , 225 , 236 , 247};
static const int rep_values_qp19[23]={0,  7, 15,  26,  36 , 48 , 60 , 72  ,  84 , 96 , 108 , 120 , 132 , 144 , 156 , 168 , 180 , 192 , 204 , 216 , 228 , 240 , 252};
static const int rep_values_qp20[20]={0,  8, 16,  28,  40 , 54 , 68 , 82  ,  96 , 110 , 124 , 138 , 152 , 166 , 180 , 194 , 208 , 222 , 236 , 250};
static const int rep_values_qp21[18]={0,  8, 20,  32,  46 , 60 , 76 , 92  , 108 , 124 , 140 , 156 , 172 , 188 , 204 , 220 , 236 , 252};
static const int rep_values_qp22[17]={0,  9, 22,  36,  52 , 68 , 84 , 101 , 118 , 135 , 152 , 169 , 186 , 203 , 220 , 237 , 254};
static const int rep_values_qp23[13]={0, 11, 28,  46,  66 , 96 , 116 , 136 , 156 , 176 , 196 , 216 , 236};
static const int rep_values_qp24[8]={0, 12, 30,  50,  80 , 130 , 180 , 230};
static const int rep_values_qp25[8]={0, 13, 33,  54,  95 , 145 , 200 , 255};
static const int rep_values_qp26[7]={0, 14, 36,  58, 110 , 170 , 230};
static const int rep_values_qp27[6]={0, 16, 42,  69, 122 , 200};
static const int rep_values_qp28[6]={0, 18, 48,  80, 134 , 210};
static const int rep_values_qp29[6]={0, 20, 54,  92, 148 , 250};
static const int rep_values_qp30[5]={0, 22, 62, 105, 172};
static const int rep_values_qp31[5]={0, 25, 69, 118, 196};
static const int rep_values_qp32[5]={0, 28, 76, 130, 220};
static const int rep_values_qp33[5]={0, 32, 84, 138, 244};
static const int rep_values_qp34[5]={0, 36, 96, 160, 255};
static const int rep_values_qp35[4]={0, 40,108, 184};
static const int rep_values_qp36[4]={0, 44,124, 210};
static const int rep_values_qp37[4]={0, 50,128, 236};
static const int rep_values_qp38[4]={0, 56,152, 255};
static const int rep_values_qp39[4]={0, 64,168, 255};
static const int rep_values_qp40[3]={0, 72,192};
static const int rep_values_qp41[3]={0, 80,216};
static const int rep_values_qp42[3]={0, 88,248};
static const int rep_values_qp43[3]={0,100,255};
static const int rep_values_qp44[3]={0,112,255};
static const int rep_values_qp45[2]={0,128};
static const int rep_values_qp46[2]={0,144};
static const int rep_values_qp47[2]={0,160};
static const int rep_values_qp48[2]={0,176};
static const int rep_values_qp49[2]={0,200};
static const int rep_values_qp50[2]={0,230};
static const int rep_values_qp51[2]={0,255};

static const int *rep_value_ptr[52]={
rep_values_qp00_04, rep_values_qp00_04, rep_values_qp00_04, rep_values_qp00_04, rep_values_qp00_04,
rep_values_qp05_09, rep_values_qp05_09, rep_values_qp05_09, rep_values_qp05_09, rep_values_qp05_09,
rep_values_qp10, rep_values_qp11, rep_values_qp12, rep_values_qp13, rep_values_qp14, rep_values_qp15,
rep_values_qp16, rep_values_qp17, rep_values_qp18, rep_values_qp19, rep_values_qp20, rep_values_qp21,
rep_values_qp22, rep_values_qp23, rep_values_qp24, rep_values_qp25, rep_values_qp26, rep_values_qp27,
rep_values_qp28, rep_values_qp29, rep_values_qp30, rep_values_qp31, rep_values_qp32, rep_values_qp33,
rep_values_qp34, rep_values_qp35, rep_values_qp36, rep_values_qp37, rep_values_qp38, rep_values_qp39,
rep_values_qp40, rep_values_qp41, rep_values_qp42, rep_values_qp43, rep_values_qp44, rep_values_qp45,
rep_values_qp46, rep_values_qp47, rep_values_qp48, rep_values_qp49, rep_values_qp50, rep_values_qp51};


static const int SD_RD_quant_s_for_all_qps_rec_int[52]={
  256,  256,  256,  256,  256,  256,  256,  256,  256,  256,  768,  768,  768,  768,  768,
  998, 1050, 1280, 1510, 1536, 1613, 1997, 2202, 2483, 2790, 3277, 3661, 4070, 4634, 4890,
 5658, 6426, 6938, 8218, 8730,10010,11341,12595,14106,16179,17357,21683,26496,35558,39245,
39731,40141,40960,43520,46080,48640,51200};

static const int SD_RD_quant_szd_for_all_qps_rec_int[52]={
  384,  384,   384,   384,   384,   384,   384,   384,   512,   512,   998,   998,   998,   998,  1254,
 1261, 1167,  1510,  1674,  2048,  2029,  2296,  2811,  2792,  3068,  3564,  4073,  4590,  5114,  5627,
 6395, 7419,  8175,  9468, 11007, 12530, 14025, 16102, 18676, 21443, 24466, 27964, 32128, 35857, 40696,
46812,46244, 49382, 52070, 54758, 55014, 60134};

int de_quantizer(int value)
{
  if (img->SD_Quantizer==0)
  {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
    return (value>0?1:-1)*(*(rep_value_ptr[Clip3(0,51,img->qp+img->bitdepth_luma_qp_scale)]+(abs(value))));
#else
    return (value>0?1:-1)*(*(rep_value_ptr[img->qp]+(abs(value))));
#endif
  }
  else
  {
    if (value==0)     return 0;
    else if (value>0)
    {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
      return   min(((1<<img->bitdepth_luma)-1),max(0,((((value>0?value:-value)-1)*SD_RD_quant_s_for_all_qps_rec_int[Clip3(0,51,img->qp+img->bitdepth_luma_qp_scale)]+SD_RD_quant_szd_for_all_qps_rec_int[Clip3(0,51,img->qp+img->bitdepth_luma_qp_scale)])>>8)));
#else
      return   min(255,max(0,((((value>0?value:-value)-1)*SD_RD_quant_s_for_all_qps_rec_int[img->qp]+SD_RD_quant_szd_for_all_qps_rec_int[img->qp])>>8)));
#endif
    }
    else
    {
#ifdef  INTERNAL_BIT_DEPTH_INCREASE
      return  -min(((1<<img->bitdepth_luma)-1),max(0,((((value>0?value:-value)-1)*SD_RD_quant_s_for_all_qps_rec_int[Clip3(0,51,img->qp+img->bitdepth_luma_qp_scale)]+SD_RD_quant_szd_for_all_qps_rec_int[Clip3(0,51,img->qp+img->bitdepth_luma_qp_scale)])>>8)));
#else
      return  -min(255,max(0,((((value>0?value:-value)-1)*SD_RD_quant_s_for_all_qps_rec_int[img->qp]+SD_RD_quant_szd_for_all_qps_rec_int[img->qp])>>8)));
#endif
    }
  }
}
#endif
