
/*!
 ************************************************************************
 * \file refbuf.c
 *
 * \brief
 *    Declarations of teh reference frame buffer types and functions
 ************************************************************************
 */

#include <stdlib.h>
#include <memory.h>
#include "global.h"

#include "refbuf.h"

/*!
 ************************************************************************
 * \brief
 *    Reference buffer write routines
 ************************************************************************
 */
void PutPel_14 (pel_t **Pic, int y, int x, pel_t val)
{
  Pic [y][x] = val;
}

void PutPel_11 (pel_t *Pic, int pel_pos, pel_t val)
{
  Pic [pel_pos] = val;
}

/*!
 ************************************************************************
 * \note
 *    The following functions returning line are NOT reentrant!  Use a buffer
 *    provided by the caller to change that (but it costs a memcpy()...
 ************************************************************************
 */
static pel_t line[16];


#ifdef EIGHTH_PEL
int i, maxx, xfloor, yfloor, xx, yy, xq, yq;
static int buf_switch;
pel_t *out;
#endif

pel_t *FastLine16Y_11 (pel_t *Pic, int y, int x, int height, int width)
{
  return &Pic [y*width+x];
}


pel_t *UMVLine16Y_11 (pel_t *Pic, int y, int x, int height, int width)
{
  int i, maxx;
  pel_t *Picy;

  Picy = &Pic [max(0,min(height-1,y)) * width];

  if (x < 0) 
  {                    // Left edge ?

    maxx = min(0,x+16);
    for (i = x; i < maxx; i++)
      line[i-x] = Picy [0];       // Replicate left edge pixel

    maxx = x+16;
    for (i = 0; i < maxx; i++)    // Copy non-edge pixels
      line[i-x] = Picy [i];
  }
  else if (x > width-16)  
  {  // Right edge ?

    maxx = width;
    for (i = x; i < maxx; i++)
      line[i-x] = Picy [i];       // Copy non-edge pixels

    maxx = x+16;
    for (i = max(width,x); i < maxx; i++)
      line[i-x] = Picy [width-1];  // Replicate right edge pixel
  }
  else                            // No edge
    return &Picy [x];

  return line;
}


pel_t *FastLineX (int dummy, pel_t* Pic, int y, int x, int height, int width)
{
  return Pic + y*width + x;
}

static pel_t umv_line[2][16];

pel_t *UMVLineX_impl (int size, pel_t* Pic, int y, int x, int height, int width, int line_idx)
{
  int i;
  pel_t *Picy;

  Picy = Pic + max(0,min(height-1,y)) * width;

  if (x < 0)                            // Left edge
  {   
    for (i = x; i < min(0,x+size); i++)
    {
      umv_line[line_idx][i-x] = Picy [0];             // Replicate left edge pixel
    }

    memcpy(&umv_line[line_idx][min(-x,15)],Picy,max(x+size,0) * sizeof(pel_t)); // Copy non-edge pixels
  }
  else if (x > width-size)         // Right edge
  {
    memcpy(umv_line[line_idx],&Picy[x], max((width - x),0) * sizeof(pel_t)); // Copy non-edge pixels
    for (i = max(width,x); i < x+size; i++)
    {
      umv_line[line_idx][i-x] = Picy [width-1];  // Replicate right edge pixel
    }
  }
  else                                  // No edge
  {
    return Picy + x;
  }

  return umv_line[line_idx];
}

pel_t *UMVLineX (int size, pel_t* Pic, int y, int x, int height, int width)
{
  return UMVLineX_impl (size, Pic, y, x, height, width, 0);
}

pel_t *UMVLineX2 (int size, pel_t* Pic, int y, int x, int height, int width)
{
  return UMVLineX_impl (size, Pic, y, x, height, width, 1);
}


/*!
 ************************************************************************
 * \brief
 *    Reference buffer, 1/4 pel
 ************************************************************************
 */
pel_t UMVPelY_14 (pel_t **Pic, int y, int x, int height4, int width4)
{
  int ypos = (y < 0 ? y&3 : (y >  height4 ? height4+(y&3) : y));
  int xpos = (x < 0 ? x&3 : (x >  width4 ? width4+(x&3) : x));

  return Pic [ypos][xpos];
}

pel_t FastPelY_14 (pel_t **Pic, int y, int x, int height, int width)
{
  return Pic [y][x];
}



pel_t *FastLine4X (pel_t **Pic, int y, int x, int height, int width)
{
  return &Pic [y][x];
}

static pel_t line4[2][64];
#ifdef EIGHTH_PEL
static pel_t line41[64];
#endif
/*!
 ************************************************************************
 * \brief
 *    Reference buffer, 1/4 pel
 ************************************************************************
 */
pel_t *UMVLine4X_impl (pel_t **Pic, int y, int x, int height4, int width4, int line_idx)
{
  int i, xx;

  int xpos, ypos = (y < 0 ? y&3 : (y >  height4 ? height4+(y&3) : y));

  for (i=0; i< 64; i+=4)
  {
    xx = x + i;
    xpos = (xx < 0 ? (xx) &3 : ((xx) >  width4 ? width4+((xx) &3) : (xx) ));
    line4[line_idx][i] = Pic [ypos][xpos];
  }
  return line4[line_idx];
}

pel_t *UMVLine4X (pel_t **Pic, int y, int x, int height4, int width4)
{
  return UMVLine4X_impl (Pic, y, x, height4, width4, 0);
}

pel_t *UMVLine4X2 (pel_t **Pic, int y, int x, int height4, int width4)
{
  return UMVLine4X_impl (Pic, y, x, height4, width4, 1);
}

#ifdef EIGHTH_PEL

pel_t *UMVLine4X8 (pel_t **Pic, int y, int x, int height4, int width4)
{
  /* Maximum values (padding included) */
  int maxx8 = ( width4 - 4 ) << 1;
  int maxy8 = (height4 - 4 ) << 1;

  // quick hack to use two different buffers for b-frames
  buf_switch = ( buf_switch == 0 ) ? 16 : 0;
  out = line41 + buf_switch;


  if (y < 0)
    yfloor = (y&7) >> 1;
  else if (y > maxy8)
    yfloor = ( maxy8 + (y&7)) >> 1;
  else
    yfloor = y>>1;

  yy = y&1;
  xx = x&1;

  for( i=0; i < 4; i++, out+=4, x+=8)
  {
    if (x < 0)
      xfloor = (x&7) >> 1;
    else if (x > maxx8)
      xfloor = (maxx8 + (x&7)) >> 1;
    else
      xfloor = x>>1;

    if( xx && yy ) // are x and y on an eighth-pel position?
      *out = ( Pic[yfloor][xfloor] + Pic[yfloor  ][xfloor+1] + Pic[yfloor+1][xfloor] + Pic[yfloor+1][xfloor+1]) >> 2;
    else if( xx ) // only x is on an eighth-pel position
      *out = ( Pic[yfloor][xfloor] + Pic[yfloor  ][xfloor+1] ) >> 1;
    else if( yy ) // only y is on an eighth-pel position
      *out = ( Pic[yfloor][xfloor] + Pic[yfloor+1][xfloor  ] ) >> 1;
    else // quarter-pel position
      *out = Pic[yfloor][xfloor];
  }
  return out-16;
}

pel_t *FastLine4X8 (pel_t **Pic, int y, int x, int height4, int width4)
{
  xfloor = x>>1; yfloor = y>>1; // quarter-pel units
  xx = x&1; yy = y&1; // is eighth-pel position?

  xq = x&7; // eighth-pel sub position
  yq = y&7;

  // quick hack to use 2 different buffers for b-frames
  buf_switch = ( buf_switch == 0 ) ? 16 : 0;
  out = line41 + buf_switch;

  for( i=0; i < 4; i++, out+=4, x+=8, xfloor+=4)
  {
    if( xx && yy ) // are x and y on an eighth-pel position?
      *out = ( Pic[yfloor][xfloor] + Pic[yfloor  ][xfloor+1] + Pic[yfloor+1][xfloor] + Pic[yfloor+1][xfloor+1]) >> 2;
    else if( xx ) // only x is on an eighth-pel position
      *out = ( Pic[yfloor][xfloor] + Pic[yfloor  ][xfloor+1] ) >> 1;
    else if( yy ) // only y is on an eighth-pel position
      *out = ( Pic[yfloor][xfloor] + Pic[yfloor+1][xfloor  ] ) >> 1;
    else // quarter-pel position
      *out = Pic[yfloor][xfloor];
  }
  return out-16;
}
#endif
