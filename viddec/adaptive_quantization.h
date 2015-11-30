/*!
************************************************************************
* \file adaptive_quantization.h
*
* \brief
*    adaptive quantization matrix selection routines
*
* \author
*  Akiyuki Tanizawa       <akiyuki.tanizawa@toshiba.co.jp>           \n
*  Takeshi Chujoh         <takeshi.chujoh@toshiba.co.jp>             \n
*
* This software may be used for research purposes only.
* TOSHIBA would also appreciate that any technical report, conference or
* journal paper which uses the software includes one of the following references:
*    - A. Tanizawa, T. Chujoh, "Adaptive Quantization Matrix Selection",
*                              VCEG contribution D.266, Geneva, November 2006.
*    - A. Tanizawa, T. Chujoh, "Simulation results of AQMS on KTA software",
*                              VCEG contribution VCEG-AC07, Klagenfurt, July 2006.
*    - A. Tanizawa, T. Chujoh, "Adaptive Quantization Matrix Selection on KTA software",
*                              VCEG contribution VCEG-AD06, Hangzhou, Oct. 2006.
*    - A. Tanizawa, T. Chujoh, "Simulation results of AQMS on KTA software version 1.3",
*                              VCEG contribution VCEG-AF08, San Jose, April 2007.
*    - A. Tanizawa, T. Chujoh, T. Yamakage, 
*                              "Encoding complexity reduction of AQMS",
*                              VCEG contribution C.403, Geneve, April 2008.
*    - A. Tanizawa, T. Chujoh, T. Yamakage, 
*                              "Improvement of Adaptive Quantization Matrix Selection",
*                              VCEG contribution VCEG-AI19, Berlin, July 2008.
************************************************************************
*/

/*
* CONFIGURATION:
*
* Configuration of the .Cfg file:
* In section KTA STUFF:
*
* // AQMS command:
* UseAdaptiveQuantMatrix    = 0 (0:disable, 1:enable (IAQMS))
*
* NOTE:
* IAQMS can perform on VCEG-AA10 conditions.
*    - TK. Tan, G. Sullivan, T. Wedi, "Recommended Simulation Common Conditions for Coding Efficiency Experiments",
*                              VCEG contribution VCEG-AA10, Nice, Oct. 2005.
*    - TK. Tan, G. Sullivan, T. Wedi, "Recommended Simulation Common Conditions for Coding Efficiency Experiments",
*                              VCEG contribution VCEG-AE10r1, Marrakech, Jan. 2007.
********************
* KNOWN LIMITATIONS
********************
* Not compatible with:                
*  - interlace coding
*  - RdOpt != 1   
*  - RateControl
*
********************
* VERSIONS
********************
* 1.0 KTA version 1.2  - First Implementation for AQMS
* 1.1 KTA version 1.3  - Modification of Compatibility with AIF, APEC, 
* 1.2 KTA version 1.9  - Addition of Fast AQMS and Modification of Compatibility with RDOQ, 
* 1.3 KTA version 2.1  - Replacement of IAQMS (Improved AQMS) and Deletion of Adaptive QP
*
*/

#ifndef _ADAPTIVE_QUANTIZATION_H_
#define _ADAPTIVE_QUANTIZATION_H_

#include "global.h"
#include "mbuffer.h"

#ifdef ADAPTIVE_QUANTIZATION

void AssignQuantParamForAQMS(pic_parameter_set_rbsp_t* pps);
void CalculateQuantParamForAQMS();
void CalculateQuant8ParamForAQMS();
void CreateModelingMatrixForIAQMS(void);
void SetAutoScalingListDataForIAQMS();
void CreateScalingList4x4ForIAQMS(int Gradient, int ListIdx);
void CreateScalingList8x8ForIAQMS(int Gradient, int ListIdx);
void InitScalingListSPS(seq_parameter_set_rbsp_t *sps);
void InitScalingListPPS(pic_parameter_set_rbsp_t *pps);
void SetAutoScalingListDataForMQM(int param0, int param1);
void CreateScalingList4x4ForMQM(int pA, int pB);
void CreateScalingList8x8ForMQM(int pA, int pB);

#endif

#endif //_ADAPTIVE_QUANTIZATION_H_

