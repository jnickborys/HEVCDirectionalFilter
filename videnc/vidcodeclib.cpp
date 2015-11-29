#include <iostream>
#include <fstream>
using namespace std;

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "codeclib.h"
#include "arith.h"
#include "bits.h"
#include "xform.h"
#include "quant.h"

//-----------------------------------
// ICodec class
//-----------------------------------
ICodec::ICodec(unsigned w, unsigned h)
  :_w(w), _h(h), _w2(w / 2), _h2(h / 2),
  m_iMBNumW (w / MBSIZE), m_iMBNumH (h / MBSIZE)
{
  unsigned i;

  //Init Image buffer: the (y, x)-th pixel of the image can be accessed  by _ydata[y][x].
  m_pfFrameBuf[0] = new float[_w * _h * 3 / 2];
  m_pfFrameBuf[1] = new float[_w * _h * 3 / 2];

  _ydata = new float *[_h];
  _udata = new float *[_h2];
  _vdata = new float *[_h2];

  _ydataRef = new float *[_h];
  _udataRef = new float *[_h2];
  _vdataRef = new float *[_h2];

  // Upsampled Data
  _ydataUp = new float *[_h * 16];
  _udataUp = new float *[_h2 * 16];
  _vdataUp = new float *[_h2 * 16];
  
  SetFrameBufPointer(m_pfFrameBuf[0], m_pfFrameBuf[1]);

  //MV buffer
  m_iMVy = new int *[m_iMBNumH];
  m_iMVy[0] = new int[m_iMBNumW * m_iMBNumH];
  for(i = 1; i < m_iMBNumH; i++) {
    m_iMVy[i] = m_iMVy[i-1] + m_iMBNumW;
  }

  m_iMVx = new int *[m_iMBNumH];
  m_iMVx[0] = new int[m_iMBNumW * m_iMBNumH];
  for(i = 1; i < m_iMBNumH; i++) {
    m_iMVx[i] = m_iMVx[i-1] + m_iMBNumW;
  }

  //context prob models
  for (i = 0; i < MAX_CONTEXT_NUM; i++) {
      _context[i] = new binModel(1, 1, _w); 
  }

}

ICodec::~ICodec(void)
{
  int i;

  delete m_pfFrameBuf[0];
  delete m_pfFrameBuf[1];

  delete _ydata;
  delete _udata;
  delete _vdata;

  delete _ydataRef;
  delete _udataRef;
  delete _vdataRef;
 
  delete m_iMVy[0];
  delete m_iMVy;
  delete m_iMVx[0];
  delete m_iMVx;

  for (i = 0; i < MAX_CONTEXT_NUM; i++) {
      delete _context[i];
  }
}

//Two set of buffers are needed for video coding:
//_ydata, _udata, _vdata hold YUV of the current frame
//_ydataRef, _udataRef, _vdataRef hold the reference frame
void ICodec::SetFrameBufPointer(
    float *pfBuf1,
    float *pfBuf2
)
{
  unsigned y;

  //_d, _udata, _vdata hold the current frame
  _ydata[0] = pfBuf1;
  for(y = 1; y < _h; y++) {
    _ydata[y] = _ydata[y-1] + _w;
  }

  _udata[0] = _ydata[0] + _w * _h;
  for(y = 1; y < _h2; y++) {
    _udata[y] = _udata[y-1] + _w2;
  }

  _vdata[0] = _udata[0] + _w2 * _h2;
  for(y = 1; y < _h2; y++) {
    _vdata[y] = _vdata[y-1] + _w2;
  }

  //_ydataRef, _udataRef, _vdataRef hold the reference frame
  _ydataRef[0] = pfBuf2;
  for(y = 1; y < _h; y++) {
    _ydataRef[y] = _ydataRef[y-1] + _w;
  }

  _udataRef[0] = _ydataRef[0] + _w * _h;
  for(y = 1; y < _h2; y++) {
    _udataRef[y] = _udataRef[y-1] + _w2;
  }

  _vdataRef[0] = _udataRef[0] + _w2 * _h2;
  for(y = 1; y < _h2; y++) {
    _vdataRef[y] = _vdataRef[y-1] + _w2;
  }

}

//Create a reconstruced frame of the current frame,
//so that the next frame can use it as reference frame.
//Current frame already contains the reconstructed prediction error.
void ICodec::GetReconstructedFrame()
{
    unsigned y, x;

    for(y = 0; y < _h; y += BLOCKSIZE) {
        for(x = 0; x < _w; x += BLOCKSIZE) {
            GetBlockRecon(_ydata, _ydataRef, y, x, m_iMVy[y / MBSIZE][x / MBSIZE], m_iMVx[y / MBSIZE][x / MBSIZE]);
        }
    }

    for(y = 0; y < _h2; y += BLOCKSIZE) {
        for(x = 0; x < _w2; x += BLOCKSIZE) {
            GetBlockRecon(_udata, _udataRef, y, x, m_iMVy[y / MBSIZEUV][x / MBSIZEUV] / 2, m_iMVx[y / MBSIZEUV][x / MBSIZEUV] / 2);
        }
    }

    for(y = 0; y < _h2; y += BLOCKSIZE) {
        for(x = 0; x < _w2; x += BLOCKSIZE) {
            GetBlockRecon(_vdata, _vdataRef, y, x, m_iMVy[y / MBSIZEUV][x / MBSIZEUV] / 2, m_iMVx[y / MBSIZEUV][x / MBSIZEUV] / 2);
        }
    }

}




void ICodec::GetBlockRecon(
    float **pfCurrFrame,    
    float **pfRefFrame,     
    int y0,     //y0, x0 give the upper-left corner of the block
    int x0, 
    int iMVy,   //iMVy, iMVx: MV.
    int iMVx)
{
    int i, j;

    for (i = 0; i < BLOCKSIZE; i++) {
        for (j = 0; j < BLOCKSIZE; j++) {
            pfCurrFrame[y0 + i][x0 + j] = pfCurrFrame[y0 + i][x0 + j] + pfRefFrame[y0 + i + iMVy][x0 + j + iMVx];
        }
    }
}


// Map signed integer to unsigned index for Golomb Rice encoding
// Positive numbers are mapped to odd indices.
// Negative numbers are mapped to even indices.
unsigned ICodec::ForwardGolombRiceIndexMapping(
    int iNumber)
{
    //Convert signed int to unsigned index:
    //iNumber:   0     1   -1  2   -2  3   -3 ...
    //uiIndex:   0     1    2  3    4  5    6 ...
    unsigned uiIndex;
    if (iNumber > 0) {
        uiIndex = 2 * iNumber - 1;
    } else {
        uiIndex = (unsigned) (-2 * iNumber);
    }
    return uiIndex;
}

// Map unsigned index back to signed integer for Golomb Rice decoding
// Odd indices are mapped to positive numbers.
// Even indices are mapped to negative numbers.
int ICodec::InverseGolombRiceIndexMapping(
    unsigned uiIndex)
{
    //Convert signed int to unsigned index:
    //uiIndex: 0     1    2  3    4  5    6 ...
    //iNumber: 0     1   -1  2   -2  3   -3 ...

    int iNumber;
    if (uiIndex & 1) {
        iNumber = (uiIndex + 1) >> 1;
    } else {
        iNumber = -1 * ((int) (uiIndex >> 1));
    }
    return iNumber;
}


//Copy image to the internal buffer
//and subtract 128.
void ICodec::SetImage(
	unsigned char *pcBuf)
{
  unsigned len = _w * _h, i;
  for(i = 0; i < len; i++) {
    _ydata[0][i] = (float) (pcBuf[i]);
  }

  pcBuf += len;
  for(i = 0; i < len / 4; i++) {
    _udata[0][i] = (float) (pcBuf[i]);
  }

  pcBuf += len / 4;
  for(i = 0; i < len / 4; i++) {
    _vdata[0][i] = (float) (pcBuf[i]);
  }
}

//Read internal image to external buffer
//After SetFrameBufPointer in CodeImage, the recently decoded image is in Ref buffer.
void ICodec::GetImage(
	unsigned char *pcBuf)
{
  int tmp;
  unsigned len = _w * _h, i;

  for(i = 0; i < len; i++) {
      //round to integer
      tmp = (int) (_ydataRef[0][i] + 0.5);

      //clipping to [0, 255]
      pcBuf[i] = (unsigned char) (tmp < 0 ? 0 : (tmp > 255 ? 255 : tmp));
  }

  pcBuf += len;
  for(i = 0; i < len / 4; i++) {
      //round to integer
      tmp = (int) (_udataRef[0][i] + 0.5);

      //clipping to [0, 255]
      pcBuf[i] = (unsigned char) (tmp < 0 ? 0 : (tmp > 255 ? 255 : tmp));
  }

  pcBuf += len / 4;
  for(i = 0; i < len / 4; i++) {
      //round to integer
      tmp = (int) (_vdataRef[0][i] + 0.5);

      //clipping to [0, 255]
      pcBuf[i] = (unsigned char) (tmp < 0 ? 0 : (tmp > 255 ? 255 : tmp));
  }

}


//-----------------------------------
// IEncoder class
//-----------------------------------

IEncoder::IEncoder(unsigned w,unsigned h)
  :ICodec(w, h)
{
  _ace = new ACEncoder();
  _out = new OFlow();
}

IEncoder::~IEncoder(void)
{
  delete _ace;
  delete _out;
}

// Perform unary code and binary arithmetic code:
// The first bit use Conext 0, and the rest bits use Context 1.
// Similar approach to H.264, since the first bit has more prob of 0.
void IEncoder::EncodeUnary(
    unsigned uiIndex)   //index to be encoded
{    
    if (uiIndex == 0) {
        //Codeword is 0: code with Context 0.
        _ace->codeSymbol(0, _context[0], _out);
    } else {
        //Other codeword has the fomat of 111...10
        //Code the first bit with Context 0.
        _ace->codeSymbol(1, _context[0], _out);
        uiIndex--;

        //Code other bits with Context 1.
        while (uiIndex > 0) {
            _ace->codeSymbol(1, _context[1], _out);
            uiIndex--;
        }

        //the last bit is 0, still use context 1
        _ace->codeSymbol(0, _context[1], _out);
    }    
}


// Perform Golomb-Rice code and binary arithmetic code:
// The first bit use Conext 0, and the rest bits use Context 1.
// Similar approach to H.264, since the first bit has more prob of 0.
void IEncoder::EncodeGolombRice(
    unsigned uiIndex,           //index to be encoded
    int iGRPara)                //Number of Golomb-Rice remainder bits
{    
    bool bBit;

    //Encode Group ID (with 1 or 2 contexts)
    unsigned uiGroup = uiIndex >> iGRPara;
    EncodeUnary(uiGroup);

    //encode remainder bits with Context 1
    for (int i = iGRPara - 1; i >= 0; i--) {
        bBit = (uiIndex >> i) & 1;
        _ace->codeSymbol(bBit, _context[1], _out);
    }
}


// Perform Exp-Golomb code and binary arithmetic code:
// The first bit use Conext 0, and the rest bits use Context 1.
// Similar approach to H.264, since the first bit has more prob of 0.
void IEncoder::EncodeExpGolomb(
    unsigned uiIndex)           //index to be encoded
{    
    bool bBit;

    if (uiIndex == 0) {
        EncodeUnary(0);
    } else {
        //Get Group ID
        unsigned uiGroup = 2;
        while ((int) uiIndex > (1 << uiGroup) - 2) {
            uiGroup++;
        }
        uiGroup = uiGroup - 1;

        //Encode Group ID (with 1 or 2 contexts)
        EncodeUnary(uiGroup);

        //encode offset within each group bits with Context 1
        unsigned uiOffset = uiIndex - (1 << uiGroup) + 1;
        for (int i = uiGroup - 1; i >= 0; i--) {
            bBit = (uiOffset >> i) & 1;
            _ace->codeSymbol(bBit, _context[1], _out);
        }
    }
}

// Main function to encode an image
int IEncoder::codeImage(
    bool bIsIFrame,
	unsigned char *pcBitstreamBuf,   //output buffer
    float fQstep)                    //quantization step size,
{
  unsigned x, y;

  //Initialize output buffer
  _out->reset(pcBitstreamBuf);

  //Initialize arithmetic coder
  _ace->start();

  if (!bIsIFrame) {
      //P frames: motion est, find prediction error

	  //AdaptiveInterpolationFilter();

	  MotionEst();

      EncodeMV();

      GetPredError();

	  //call custom DAIF function here
	  //pretty sure it is called for P and B frames
	  //DAIF()

  }


  //encode Y component: _d contains reconstructed result after codeBlock().
  for(y = 0; y < _h; y += BLOCKSIZE) {
      for(x = 0; x < _w; x += BLOCKSIZE) {
          codeBlock(_ydata, y, x, fQstep);
      }
  }

  //encode U component
  for(y = 0; y < _h2; y += BLOCKSIZE) {
      for(x = 0; x < _w2; x += BLOCKSIZE) {
          codeBlock(_udata, y, x, fQstep);
      }
  }

  //encode V component
  for(y = 0; y < _h2; y += BLOCKSIZE) {
      for(x = 0; x < _w2; x += BLOCKSIZE) {
          codeBlock(_vdata, y, x, fQstep);
      }
  }

  if (!bIsIFrame) {
      //Get reconstructed frame, will be used as reference for the next frame
      GetReconstructedFrame();
  }

  //swap pointers:
  //the next frame will be written into the ref frame of the current frame,
  //and the current frame will become the ref frame of the next frame.
  SetFrameBufPointer(_ydataRef[0], _ydata[0]);
  
  _ace->stop(_out);
  return _out->bytesUsed();
}

//Encode an image block
void IEncoder::codeBlock(
    float **buf,    //Y/U/V buffer header,
    unsigned y,     //y index of upper-left corner of the block
    unsigned x,     //x index of upper-left corner of the block
    float fQstep)   //quantization step size,
{
    unsigned i, j, uiIndex;
    bool skipblk = false;

    //update m_fBlkBuf so that it points to the current block
    for (i = 0; i < BLOCKSIZE; i++) {
        m_fBlkBuf[i] = buf[y + i] + x;
    }
  
    //Forward DCT
    if (BLOCKSIZE == 4) {
        Transform::FDCT4(m_fBlkBuf);
    } else {
        Transform::FDCT8(m_fBlkBuf);
    }

    //Quantization
    Quant::QuantMidtread(m_fBlkBuf, BLOCKSIZE, fQstep);

    //check whether there is any non-zero coeff, 
    //send 1 bit flag to signal this: 1: skip, 0: not skip, encode all coeffs after this.
    //this is a simplified version of CBP, but we send it at the beginning of each block.
    skipblk = true;
    for (i = 0; i < BLOCKSIZE; i++) {
        for (j = 0; j < BLOCKSIZE; j++) {
            if (m_fBlkBuf[i][j] != 0) {
                skipblk = false;
                break;
            }
        }
        if (!skipblk) {
            break;
        }
    }
    //encode the skip bit using context 2.
    _ace->codeSymbol(skipblk, _context[2], _out);

    //encode all coeffs of m_fBlkBuf, not efficient since we encode all zeros at the end.
    if (!skipblk) {
        for (i = 0; i < BLOCKSIZE; i++) {
            for (j = 0; j < BLOCKSIZE; j++) {
                //convert signed int to unsigned index for G-R code
                uiIndex = ForwardGolombRiceIndexMapping((int) m_fBlkBuf[i][j]);

                //encode with unary code and arithmetic code
                EncodeUnary(uiIndex);
            }
        }
    }

    //Get recon for prediction
    //Dequantization
    Quant::DequantMidtread(m_fBlkBuf, BLOCKSIZE, fQstep);

    //IDCT
    if (BLOCKSIZE == 4) {
        Transform::IDCT4(m_fBlkBuf);
    } else {
        Transform::IDCT8(m_fBlkBuf);
    }

}

//find MV of each macro-block
void IEncoder::MotionEst()
{
    unsigned y, x;

    for(y = 0; y < _h; y += MBSIZE) {
        for(x = 0; x < _w; x += MBSIZE) {
            MBMotionEst(_ydata, _ydataRef, y, x);
        }
    }


}


void IEncoder::MBMotionEst(
    float **pfCurrFrame,    
    float **pfRefFrame,
    int y,
    int x)
{
    int i, j, k, l, iSAD0, iSAD = 0;
    int iMVy = 0, iMVx = 0;

    //Get SAD at MV=(0, 0) first for early termination purpose.
    iSAD0 = GetSAD(pfCurrFrame, pfRefFrame, y, x, 0, 0, 65535);
    
    //------------------------
    // HW5: full-search motion estimation:
    // For each MV candidate, call GetSAD() below to get its prediction error (SAD).
    // Only update MV if the new SAD is less than certain ratio of previous minimum SAD (saved in iSAD0).
    // Skip invalid MV (any pixel in the reference block is out of the frame boundary).
    // Keep track of best candidate with iMVy and iMVx.
    //------------------------

	for(i = -18; i <= 18; i++) {     //iterate through search window
		for(j = -18; j <= 18; j++) { //iterate through search window

			//boundary conditions
		if( y+i < 0						||  
			(y+i+MBSIZE-1) >= (int) _h  ||  
			x+j < 0						|| 
			(x+j+MBSIZE-1) >= (int) _w  || 
			(i == 0 && j == 0)   ) { 
			continue; //go to next iteration of loop, just ignore this one
		}

					//compute SAD
					iSAD = GetSAD(pfCurrFrame, pfRefFrame, y, x, i, j, iSAD0);  

					//check if new SAD is lower than old SAD
					//if it is then update MV's
					if( iSAD < (0.925 * iSAD0) ) {
						iMVy = i;
						iMVx = j;

						//update to new SAD
						iSAD0 = iSAD;
					}
		}
	}

    m_iMVy[y / MBSIZE][x / MBSIZE] = iMVy;
    m_iMVx[y / MBSIZE][x / MBSIZE] = iMVx;
}

//return the sum of absolute difference (SAD) of two blocks with the given motion vectors
int IEncoder::GetSAD(
    float **pfCurrFrame,    //pointer to the current frame,
    float **pfRefFrame,     //pointer to the reference frame,
    int y,                  //(y, x) is the upper-left corner of the current block
    int x,                  
    int iMVy,               //(iMVy, iMVx) is the given MV
    int iMVx,
    int iSAD0)              //Minimum SAD so far, used for early termination.
{
    int iSAD = 0;
	int i, j;

    //------------------------
    //HW5: compute SAD, compare with the input SAD0 for early termination.
    //------------------------

		//implement SAD algorithm with early termination
		for (i = 0; i < MBSIZE && iSAD <= iSAD0; i++) {
			for (j = 0; j < MBSIZE && iSAD <= iSAD0; j++) {
				iSAD = iSAD + abs( pfCurrFrame[y + i][x + j] - pfRefFrame[y + i + iMVy][x + j + iMVx] );
			}
		}

    return iSAD;
}


/*
	Upsamples the Image and the applies a filter

	Image size is increased by 16x
*/
void IEncoder::AdaptiveInterpolationFilter()
{
	// Upsample by 2

	// apply filter h1
	// Upsample by 2
	// apply filter h1
	// upsample by 4
	// apply filter h2
	
	return ;
}

void ComputeWienerHopfEqn()
{
	// For a 6 Tap Filter
	for (int t = 0; t < 6; ++t)
	{

	}
}

void IEncoder::ComputeFilterCoefficients()
{
	return;
}

//custom DAIF (directional adaptive interpolation filter) function
void IEncoder::DAIF()
{
	//declare variables
	unsigned y, x;
	double		a, b, c, 
		     d, e, f, g, 
			 h, i, j, k, 
			 l, m, n, o;

	//go through Y frame
	//not sure if the same needs to be done for U and V frames
	for(y = 0; y < _h; y++) {     //iterate through frame in y-direction
		for(x = 0; x < _w; x++) { //iterate through frame in x-direction

			// for now we`re going to ignore the boundary conditions where pixels are closer
			// than 6 pixels to the edge, in which case the coefficients are going to
			// be determined a bit differently
			
			//_ydata means current frame, store sub-pixel locations:	
			a, c = _ydata[y][x-2] + _ydata[y][x-1] + _ydata[y][x] + _ydata[y][x+1] + _ydata[y][x+2] + _ydata[y][x+3];

			b = (_ydata[y][x-2] + _ydata[y][x+3]) + (_ydata[y][x-1] + _ydata[y][x+2]) + (_ydata[y][x] + _ydata[y][x+1]);

			d, l = _ydata[y-3][x] + _ydata[y-2][x] + _ydata[y-1][x] + _ydata[y][x] + _ydata[y+1][x] + _ydata[y+2][x];

			h = (_ydata[y-3][x] + _ydata[y+2][x]) + (_ydata[y-2][x] + _ydata[y+1][x]) + (_ydata[y-1][x] + _ydata[y][x]);

			e, g = _ydata[y+2][x-2] + _ydata[y+1][x-1] + _ydata[y][x] + _ydata[y-1][x+1] + _ydata[y-2][x+2] + _ydata[y-3][x+3];

			m, o = _ydata[y+2][x+3] + _ydata[y+1][x+2] + _ydata[y][x+1] + _ydata[y-1][x] + _ydata[y-2][x-1] + _ydata[y-3][x-2];

			j = (_ydata[y+2][x-2] + _ydata[y-3][x+3] + _ydata[y+2][x+3] + _ydata[y-3][x-2])  +
				(_ydata[y+1][x-1] + _ydata[y-2][x+2] + _ydata[y+1][x+2] + _ydata[y-2][x-1])  +
				(_ydata[y][x] + _ydata[y-1][x+1] + _ydata[y][x+1] + _ydata[y-1][x]) ;
			
			f, n = (_ydata[y+2][x-2] + _ydata[y+2][x+3])  + (_ydata[y+1][x-1] + _ydata[y+1][x+2])  +
				   (_ydata[y][x] + _ydata[y][x+1])  + (_ydata[y-1][x] + _ydata[y-1][x+1])  +
				   (_ydata[y-2][x-1] + _ydata[y-2][x+2])  + (_ydata[y-3][x-2] + _ydata[y-3][x+3]) ;

			i, k = (_ydata[y+2][x-2] + _ydata[y-3][x-2])  + (_ydata[y+1][x-1] + _ydata[y-2][x-1])  +
				   (_ydata[y][x] + _ydata[y-1][x])  + (_ydata[y][x+1] + _ydata[y-1][x+1])  +
				   (_ydata[y+1][x+2] + _ydata[y-2][x+2])  + (_ydata[y+2][x+3] + _ydata[y-3][x+3]) ;

		}//end inner for loop
	}//end outer for loop
}

//Save the MV information in a separate file for Matlab to plot.
void IEncoder::DumpMV(ofstream& DumpFile)
{
    DumpFile.write((const char *) m_iMVy[0], _w * _h / MBSIZE / MBSIZE * sizeof(int));
    DumpFile.write((const char *) m_iMVx[0], _w * _h / MBSIZE / MBSIZE * sizeof(int));
}

void IEncoder::EncodeMV()
{
    unsigned y, x;
    unsigned uiIndex;

    for(y = 0; y < m_iMBNumH; y++) {
        for(x = 0; x < m_iMBNumW; x++) {
            uiIndex = ForwardGolombRiceIndexMapping(m_iMVy[y][x]);
            EncodeUnary(uiIndex);
            uiIndex = ForwardGolombRiceIndexMapping(m_iMVx[y][x]);
            EncodeUnary(uiIndex);
         }
    }
}

//Called by encoder after motion est to get the prediction error.
//The error will then be encoded.
void IEncoder::GetPredError()
{
    unsigned y, x;

    for(y = 0; y < _h; y += BLOCKSIZE) {
        for(x = 0; x < _w; x += BLOCKSIZE) {
            GetBlockPredError(_ydata, _ydataRef, y, x, m_iMVy[y / MBSIZE][x / MBSIZE], m_iMVx[y / MBSIZE][x / MBSIZE]);
        }
    }

    for(y = 0; y < _h2; y += BLOCKSIZE) {
        for(x = 0; x < _w2; x += BLOCKSIZE) {
            GetBlockPredError(_udata, _udataRef, y, x, m_iMVy[y / MBSIZEUV][x / MBSIZEUV] / 2, m_iMVx[y / MBSIZEUV][x / MBSIZEUV] / 2);
        }
    }

    for(y = 0; y < _h2; y += BLOCKSIZE) {
        for(x = 0; x < _w2; x += BLOCKSIZE) {
            GetBlockPredError(_vdata, _vdataRef, y, x, m_iMVy[y / MBSIZEUV][x / MBSIZEUV] / 2, m_iMVx[y / MBSIZEUV][x / MBSIZEUV] / 2);
        }
    }

}


void IEncoder::GetBlockPredError(
    float **pfCurrFrame,    
    float **pfRefFrame,     
    int y0,     //y0, x0 give the upper-left corner of the block
    int x0, 
    int iMVy,   //iMVy, iMVx: MV.
    int iMVx)
{
    //HW4: Get prediction error of the current frame
    int i, j;

    for (i = 0; i < BLOCKSIZE; i++) {
        for (j = 0; j < BLOCKSIZE; j++) {
            pfCurrFrame[y0 + i][x0 + j] = pfCurrFrame[y0 + i][x0 + j] - pfRefFrame[y0 + i + iMVy][x0 + j + iMVx];
        }
    }
}

//-----------------------------------
// IDecoder class
//-----------------------------------
IDecoder::IDecoder(unsigned w,unsigned h)
  :ICodec(w,h)
{
  _acd=new ACDecoder();
  _in=new IFlow();
}

IDecoder::~IDecoder(void)
{
  delete _acd;
  delete _in;
}

// Decode unary code and binary arithmetic code:
// The first bit use Conext 0, and the rest bits use Context 1.
// Similar approach to H.264, since the first bit has more prob of 0.
unsigned IDecoder::DecodeUnary()
{
  unsigned uiIndex = 0;
  bool nextbit = _acd->decodeSymbol(_context[0], _in);

  if (nextbit == 1) {
    uiIndex = 0;
    do {
      nextbit = _acd->decodeSymbol(_context[1], _in);
      uiIndex++;
    } while (nextbit != 0);

  }

  return uiIndex;
}


// Decode Golomb-Rice code and binary arithmetic code:
// The first bit use Conext 0, and the rest bits use Context 1.
// Similar approach to H.264, since the first bit has more prob of 0.
unsigned IDecoder::DecodeGolombRice(
    int iGRPara)      //Number of Golomb-Rice remainder bits
{

    bool nextbit;

    //Decode Group ID
    unsigned uiIndex = DecodeUnary();

    //Decode the remainder bits
    for (unsigned char i = 0; i < iGRPara; i++) {
        nextbit = _acd->decodeSymbol(_context[1], _in);
        uiIndex = (uiIndex << 1) + (unsigned) nextbit;
    }

    return uiIndex;
}

// Decode Exp-Golomb code and binary arithmetic code:
// The first bit use Conext 0, and the rest bits use Context 1.
// Similar approach to H.264, since the first bit has more prob of 0.
unsigned IDecoder::DecodeExpGolomb()
{

    bool nextbit;

    //Decode Group ID
    unsigned uiGroup = DecodeUnary();

    //Decode the offset within each group
    unsigned uiIndex = 0;
    for (unsigned char i = 0; i < uiGroup; i++) {
        nextbit = _acd->decodeSymbol(_context[1], _in);
        uiIndex = (uiIndex << 1) + (unsigned) nextbit;
    }
    uiIndex += ((1 << uiGroup) - 1);

    return uiIndex;
}


// Main function to decode an image
int IDecoder::decodeImage(
    bool bIsIFrame,
	unsigned char *b,    //input buffer
    float fQstep)        //quantization step size,
{
  unsigned x, y;

  //Initialize _in to input buffer
  _in->reset(b);

  //Initialize arithmetic decoder
  _acd->start(_in);

  //HW4:
  if (!bIsIFrame) {
      //decode MV for P frames
      DecodeMV();
  }

  for(y = 0; y < _h; y += BLOCKSIZE) {
      for(x = 0; x < _w; x += BLOCKSIZE) {
          decodeBlock(_ydata, y, x, fQstep);
      }
  }

  for(y = 0; y < _h2; y += BLOCKSIZE) {
      for(x = 0; x < _w2; x += BLOCKSIZE) {
          decodeBlock(_udata, y, x, fQstep);
      }
  }

  for(y = 0; y < _h2; y += BLOCKSIZE) {
      for(x = 0; x < _w2; x += BLOCKSIZE) {
          decodeBlock(_vdata, y, x, fQstep);
      }
  }

  if (!bIsIFrame) {
      //add back prediction for P frames
      GetReconstructedFrame();
  }

  //swap pointers:
  //the next frame will be written into the ref frame of the current frame,
  //and the current frame will become the ref frame of the next frame.
  SetFrameBufPointer(_ydataRef[0], _ydata[0]);

  return _in->bytesUsed();
}


//Decode intra block
void IDecoder::decodeBlock(
    float **buf, 
    unsigned y,     //y index of upper-left corner of the block
    unsigned x,     //x index of upper-left corner of the block
    float fQstep)   //quantization step size,
{
    unsigned i, j, uiIndex;
    bool skipblk = false;

    //update m_fBlkBuf so that it points to the current block
    for (i = 0; i < BLOCKSIZE; i++) {
        m_fBlkBuf[i] = buf[y + i] + x;
    }

    //decode the skip bit using context 2.
    skipblk = _acd->decodeSymbol(_context[2], _in);

    if (skipblk) {
        //if skip bit = 1, reset block buffer, skip inverse quant and inverse DCT.
        for (i = 0; i < BLOCKSIZE; i++) {
            for (j = 0; j < BLOCKSIZE; j++) {
                m_fBlkBuf[i][j] = 0;
            }
        }
    } else {
        //decode all coeffs of the current block into m_fBlkBuf:
        for (i = 0; i < BLOCKSIZE; i++) {
            for (j = 0; j < BLOCKSIZE; j++) {
                uiIndex = DecodeUnary();

                //convert to signed int
                m_fBlkBuf[i][j] = (float) InverseGolombRiceIndexMapping(uiIndex);
            }
        }

        //Dequantization
        Quant::DequantMidtread(m_fBlkBuf, BLOCKSIZE, fQstep);

        //IDCT
        if (BLOCKSIZE == 4) {
            Transform::IDCT4(m_fBlkBuf);
        } else {
            Transform::IDCT8(m_fBlkBuf);
        }
    }
}


void IDecoder::DecodeMV()
{
    unsigned y, x;
    unsigned uiIndex;

    for(y = 0; y < m_iMBNumH; y++) {
        for(x = 0; x < m_iMBNumW; x++) {
            uiIndex = DecodeUnary();
            m_iMVy[y][x] = InverseGolombRiceIndexMapping(uiIndex);

            uiIndex = DecodeUnary();
            m_iMVx[y][x] = InverseGolombRiceIndexMapping(uiIndex);
         }
    }
}

