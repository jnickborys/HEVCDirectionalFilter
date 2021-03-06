#ifdef EDAIF2
#include "aif_common.h"
#endif
#ifdef DIRECTIONAL_FILTER
/*
	Main idea: 1D 6-tap adaptive Wiener nterpolation filter with possible directionality
	Benefits: 
	1. Reduced encoding complexity.
	2. Reduced decoding complexity.
	3. Less bit-overhead.

	Version 1:
	Key 1: 
		Keep the 1D filter for a,b,c,d,h and l sub-pels
	Key 2: 
		Assign the 1D diagonal AIF for e,g,m,o,j (NW-SE)
	Key 3: 
		Keep 2D 6x6 filters for f,j,k,n

	Results: Good enough for container, glasgow, Foreman qic (+), foreman CIF(-0.03)
		Failes with mobile!


*/		

static int TwoDEquationPattern_v4[15][SQR_FILTER]  =		// get equation number for symmetric filter coefficients
{
	{ 0,  1,  2,  3,  4,  5,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,        // a_pos
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0  },
	{ 0,  1,  2,  2,  1,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,        // b_pos
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0  },
	{ 5,  4,  3,  2,  1,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,        // c_pos
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0  },
	{ 0,  1,  2,  3,  4,  5,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,        // d_pos
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0  },
	{ 0,  0,  0,  0,  0,  0,
	0,  1,  0,  0,  0,  0,
	0,  0,  2,  0,  0,  0,        // e_pos - diagonal NW-SE filter
	0,  0,  0,  3,  0,  0,
	0,  0,  0,  0,  4,  0,
	0,  0,  0,  0,  0,  5  },
	{ 0,  0,  0,  0,  0,  0,
	0,  1,  0,  0,  1,  0,
	0,  0,  2,  2,  0,  0,        //f_pos - 
	0,  0,  3,  3,  0,  0,
	0,  4,  0,  0,  4,  0,
	5,  0,  0,  0,  0,  5  },
	{ 0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  1,  0,
	0,  0,  0,  2,  0,  0,        //g_pos - diagonal SW-NE filter
	0,  0,  3,  0,  0,  0,
	0,  4,  0,  0,  0,  0,
	5,  0,  0,  0,  0,  0  },
	{ 0,  1,  2,  2,  1,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,        // h_pos
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0  },
	{ 0,  0,  0,  0,  0,  5,
	0,  1,  0,  0,  4,  0,
	0,  0,  2,  3,  0,  0,        // i_pos - diagonal NW-SE filter
	0,  0,  2,  3,  0,  0,
	0,  1,  0,  0,  4,  0,
	0,  0,  0,  0,  0,  5  },
	{ 0,  0,  0,  0,  0,  0,
	0,  1,  0,  0,  1,  0,
	0,  0,  2,  2,  0,  0,        // j_pos - diagonal NW-SE filter
	0,  0,  2,  2,  0,  0,
	0,  1,  0,  0,  1,  0,
	0,  0,  0,  0,  0,  0  },
	{ 5,  0,  0,  0,  0,  0,
	0,  4,  0,  0,  1,  0,
	0,  0,  3,  2,  0,  0,        // k_pos - diagonal NW-SE filter
	0,  0,  3,  2,  0,  0,
	0,  4,  0,  0,  1,  0,
	5,  0,  0,  0,  0,  0  },
	{ 5,  4,  3,  2,  1,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,         // l_pos
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0  },
	{ 0,  0,  0,  0,  0,  5,
	0,  0,  0,  0,  4,  0,
	0,  0,  0,  3,  0,  0,        // m_pos - diagonal SW-NE filter
	0,  0,  2,  0,  0,  0,
	0,  1,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0  },
	{ 5,  0,  0,  0,  0,  5,
	0,  4,  0,  0,  4,  0,
	0,  0,  3,  3,  0,  0,        // j_pos - diagonal NW-SE filter
	0,  0,  2,  2,  0,  0,
	0,  1,  0,  0,  1,  0,
	0,  0,  0,  0,  0,  0  },
	{ 5,  0,  0,  0,  0,  0,
	0,  4,  0,  0,  0,  0,
	0,  0,  3,  0,  0,  0,        // o_pos - diagonal SE-NW filter
	0,  0,  0,  2,  0,  0,
	0,  0,  0,  0,  1,  0,
	0,  0,  0,  0,  0,  0  },
};

static double STANDARD_2D_FILTER_v4[15][SQR_FILTER] = {
{ 0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	1.0/64.0  ,   -5.0/64.0  ,   52.0/64.0  ,   20.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0
}, // a_pos
{ 0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	1.0/32.0  ,   -5.0/32.0  ,   20.0/32.0  ,   20.0/32.0  ,   -5.0/32.0  ,   1.0/32.0  ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0
}, // b_pos
{ 0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	1.0/64.0  ,   -5.0/64.0  ,   20.0/64.0  ,   52.0/64.0  ,   -5.0/64.0  ,   1.0/64.0  ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0
}, // c_pos
{ 0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   52.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   20.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0
}, // d_pos
{ 1.0/64.0    ,    0.0       ,    0.0		,    0.0       ,		0.0       ,   0.0       ,
	0.0       ,  -5.0/64.0   ,    0.0		,    0.0       ,		0.0       ,   0.0       ,
	0.0		  ,    0.0		 ,   52.0/64.0  ,   0.0		   ,	    0.0		  ,   0.0		,
	0.0       ,    0.0       ,    0.0		,    20.0/64.0 ,		0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0		,    0.0       ,    -5.0/64.0     ,   0.0       ,
	0.0       ,    0.0       ,    0.0		,    0.0       ,		0.0       ,   1.0/64.0
	}, // e_pos
{ 1.0/128.0,			0.0	,		0.0		,		0.0		,		0.0		,		1.0/128.0,
	0.0		,   -5.0/128.0	,	 0.0		,		0.0		,   -5.0/128.0	,			0.0 ,
	0.0		, 0.0			,   52/128.0	,		52/128.0,		0.0		,			0.0 ,
	0.0		, 0.0			,   20/128.0	,		20/128.0,		0.0		,			0.0 ,
	0.0		,   -5.0/128.0	,		0.0		,		0.0		,	 -5.0/128.0	,			0.0 ,
  1.0/128.0	,   0.0 	   ,	0.0			,		0.0		,		0.0		,		1.0/128.0
}, // f_pos
{   0.0       ,    0.0		,    0.0        ,    0.0       ,	  0.0         ,	1.0/64.0    ,
	0.0       ,    0.0       ,    0.0		,    0.0       ,	-5.0/64.0	  ,   0.0       ,
	0.0		  ,    0.0		 ,    0.0	    ,	52.0/64.0  ,		0.0		  ,   0.0		,
	0.0       ,    0.0       ,   20.0/64.0  ,	 0.0	   ,    	0.0       ,   0.0       ,
	0.0       ,    -5.0/64.0 ,    0.0		,    0.0       ,		0.0       ,   0.0       ,
	1.0/64.0  ,    0.0       ,    0.0		,    0.0       ,		0.0       ,   0.0       
	}, // g_pos
{ 0.0       ,    0.0       ,    1.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   -5.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   20.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   20.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   -5.0/32.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    1.0/32.0  ,    0.0       ,    0.0       ,   0.0
}, // h_pos
{  1.0/128.0,			0.0	,		0.0		,		0.0		,		0.0		,		1.0/128.0,
	0.0		,   -5.0/128.0	,	 0.0		,		0.0		,   -5.0/128.0	,			0.0 ,
	0.0		, 0.0			,   52/128.0	,		20/128.0,		0.0		,			0.0 ,
	0.0		, 0.0			,   52/128.0	,		20/128.0,		0.0		,			0.0 ,
	0.0		,   -5.0/128.0	,		0.0		,		0.0		,	 -5.0/128.0	,			0.0 ,
  1.0/128.0	,   0.0 	   ,	0.0			,		0.0		,		0.0		,		1.0/128.0
	},// i_pos
{ 1.0/64.0    ,    0.0       ,    0.0		,    0.0       ,		0.0       ,   1.0/64.0  ,
	0.0       ,  -5.0/64.0   ,    0.0		,    0.0       ,	-5.0/64.0     ,   0.0       ,
	0.0		  ,    0.0		 ,   20.0/64.0  ,   20.0/64.0  ,	    0.0		  ,   0.0		,
	0.0       ,    0.0       ,   20.0/64.0	,   20.0/64.0  ,		0.0       ,   0.0       ,
	0.0       ,  -5.0/64.0   ,    0.0		,    0.0       ,    -5.0/64.0     ,   0.0       ,
  1.0/64.0    ,    0.0       ,    0.0		,    0.0       ,		0.0       ,   1.0/64.0
	},//new j_pos
// j_pos
{  1.0/128.0,			0.0	,		0.0		,		0.0		,		0.0		,		1.0/128.0,
	0.0		,   -5.0/128.0	,	 0.0		,		0.0		,   -5.0/128.0	,			0.0 ,
	0.0		, 0.0			,   20/128.0	,		52/128.0,		0.0		,			0.0 ,
	0.0		, 0.0			,   20/128.0	,		52/128.0,		0.0		,			0.0 ,
	0.0		,   -5.0/128.0	,		0.0		,		0.0		,	 -5.0/128.0	,			0.0 ,
  1.0/128.0	,   0.0 	   ,	0.0			,		0.0		,		0.0		,		1.0/128.0
	}, // k_pos
{ 0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   20.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   52.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   -5.0/64.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    1.0/64.0  ,    0.0       ,    0.0       ,   0.0
}, // l_pos
{   0.0       ,    0.0		,    0.0        ,    0.0       ,	  0.0         ,	1.0/64.0    ,
	0.0       ,    0.0       ,    0.0		,    0.0       ,	-5.0/64.0	  ,   0.0       ,
	0.0		  ,    0.0		 ,    0.0	    ,	20.0/64.0  ,		0.0		  ,   0.0		,
	0.0       ,    0.0       ,    52.0/64.0 ,	 0.0	   ,    	0.0       ,   0.0       ,
	0.0       ,    -5.0/64.0 ,    0.0		,    0.0       ,		0.0       ,   0.0       ,
	1.0/64.0  ,    0.0       ,    0.0		,    0.0       ,		0.0       ,   0.0       
	}, // m_pos
{  1.0/128.0,			0.0	,		0.0		,		0.0		,		0.0		,		1.0/128.0,
	0.0		,   -5.0/128.0	,	 0.0		,		0.0		,   -5.0/128.0	,			0.0 ,
	0.0		, 0.0			,   20/128.0	,		20/128.0,		0.0		,			0.0 ,
	0.0		, 0.0			,   52/128.0	,		52/128.0,		0.0		,			0.0 ,
	0.0		,   -5.0/128.0	,		0.0		,		0.0		,	 -5.0/128.0	,			0.0 ,
  1.0/128.0	,   0.0 	   ,	0.0			,		0.0		,		0.0		,		1.0/128.0
	}, // n_pos
{ 1.0/64.0    ,    0.0       ,    0.0		,    0.0       ,		0.0       ,   0.0       ,
	0.0       ,  -5.0/64.0   ,    0.0		,    0.0       ,		0.0       ,   0.0       ,
	0.0		  ,    0.0		 ,    20.0/64.0 ,    0.0	   ,	    0.0		  ,   0.0		,
	0.0       ,    0.0       ,    0.0		,    52.0/64.0 ,		0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0		,    0.0       ,    -5.0/64.0     ,   0.0       ,
	0.0       ,    0.0       ,    0.0		,    0.0       ,		0.0       ,   1.0/64.0
	}, // o_pos
};

static int TwoDSymmetricPattern_v4[15][SQR_FILTER]  =		// get one filter from another one, if symmetry properties are used (used e.g. in ExtendFilterCoefficients)
{
	{ 0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0,
		12, 13, 14, 15, 16, 17,        // a_pos  - not used actually
		0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0  },
	{ 0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0,
		12, 13, 14, 14, 13, 12,        // b_pos
		0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0  },
	{ 0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0,
		17, 16, 15, 14, 13, 12,        // c_pos
		0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0  },
	{ 0,  0, 12,  0,  0,  0,
		0,  0, 13,  0,  0,  0,
		0,  0, 14,  0,  0,  0,        // d_pos
		0,  0, 15,  0,  0,  0,
		0,  0, 16,  0,  0,  0,
		0,  0, 17,  0,  0,  0  },
	{ 0, 0,  0,  0,  0,  0,
	  0, 7,  0,  0,  0,  0,
	  0, 0, 14,  0,  0,  0,        // e_pos
	  0, 0,  0, 21,  0,  0,
	  0, 0,  0,  0, 28,  0,
	  0, 0,  0,  0,  0, 35  },
	{ 0,  0,  0,  0,  0,  0,
		0,  7,  0,  0,  7,  0,
		0,  0, 14, 14,  0,  0,        // f_pos
		0,  0, 21, 21,  0,  0,
		0, 28,  0,  0, 28,  0,
		35, 0,  0,  0,  0, 35  },
	{   0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  7,  0,
		0,  0,  0, 14,  0,  0,        // g_pos
		0,  0, 21,  0,   0,  0,
		0, 28,  0,  0,   0,  0,
		35, 0,  0,  0,   0,  0  },
	{ 0,  0, 12,  0,  0,  0,
		0,  0, 13,  0,  0,  0,
		0,  0, 14,  0,  0,  0,        // h_pos
		0,  0, 14,  0,  0,  0,
		0,  0, 13,  0,  0,  0,
		0,  0, 12,  0,  0,  0  },
	{   0,  0,  0,  0,  0, 30,
		0,  7,  0,  0, 25,  0,
		0,  0, 14, 20,  0,  0,        // i_pos
		0,  0, 14, 20,  0,  0,
		0,  7,  0,  0, 25,  0,
		0,  0,  0,  0,  0, 30  },
	{ 0, 0,  0,  0,  0,  0,
	  0, 7,  0,  0,  7,  0,
	  0, 0, 14, 14,  0,  0,        // j_pos
	  0, 0, 14, 14,  0,  0,
	  0, 7,  0,  0,  7,  0,
	  0, 0,  0,  0,  0, 0  },
	{30,  0,  0,  0,  0,  0,
	  0, 25,  0,  0,  7,  0,
	  0,  0, 20, 14,  0,  0,   			// k_pos
	  0,  0, 20, 14,  0,  0,
	  0, 25,  0,  0,  7,  0,
	 30,  0,  0,  0,  0,  0  },

	{ 0,  0, 17,  0,  0,  0,
		0,  0, 16,  0,  0,  0,
		0,  0, 15,  0,  0,  0,         // l_pos
		0,  0, 14,  0,  0,  0,
		0,  0, 13,  0,  0,  0,
		0,  0, 12,  0,  0,  0  },
	{	0,	0,  0,  0,  0,  35,
		0,  0,  0,  0, 28,  0,
		0,  0,  0, 21,  0,  0,        // m_pos
		0,  0, 14,  0,  0,  0,
		0,  7,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0 },

	{30,  0,  0,  0,  0, 30,
      0, 25,  0,  0, 25,  0,
	  0,  0, 20, 20,  0,  0,   			// n_pos
	  0,  0, 14, 14,  0,  0,
	  0,  7,  0,  0,  7,  0,
	  0,  0,  0,  0,  0,  0  },
	{35, 0,  0,   0,  0,  0,
	 0, 28,  0,   0,  0,  0,
	 0,  0, 21,   0,  0,  0,   			// o_pos
	 0,  0,  0,  14,  0,  0,
	 0,  0,  0,   0,  7,  0,
	 0,  0,  0,   0,  0,  0	}
};


int SymmetryPosition_v4[15] = {1,1,0,0,1,1,0,0,0,1,0,0,0,0,0};	// if 0, the position is copied from another one
int Is1DPosition_v4[15] =	  {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; // if 0, the position is a 2D one
int IsDiagonal1D_v4[15] =	  {0,0,0,0,1,3,2,0,3,3,3,0,2,3,1}; // if 1, the filter alighned NW-SE, 2 - NE-SW, 3 - diagonal cross

int POS_EQUATION_NUMBER_v4[15] =  {    6, 3, 6,
									6, 6, 6, 6,
									3, 6, 3, 6,
									6, 6, 6, 6};

int FILTER_NEW_SUB_POS_v4[15] = {       
		   a_pos, b_pos, a_pos,
	a_pos, e_pos, f_pos, e_pos,
	b_pos, f_pos, j_pos, f_pos,
	a_pos, e_pos, f_pos, e_pos };


#ifdef EAIF


static int TwoDEquationPattern_EAIF[15][SQR_FILTER]  =		// get equation number for symmetric filter coefficients
{
	{ 0,  1,  2,  3,  4,  5,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,        // a_pos
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0  
  ,0
  },
	{ 0,  1,  2,  2,  1,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,        // b_pos
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0
  ,0
  },
	{ 5,  4,  3,  2,  1,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,        // c_pos
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0
  ,0
  },
	{ 0,  1,  2,  3,  4,  5,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,        // d_pos
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0
  ,0
  },
	{ 0,  0,  0,  0,  0,  0,
	0,  0,  0,  1,  0,  0,
	0,  2,  3,  4,  5,  0,        // e_pos 
	0,  6,  7,  8,  9,  0,
	0,  0, 10, 11,  0,  0,
	0,  0,  0,  0,  0,  0
  ,0
  },
	{ 0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,
	0,  1,  2,  2,  1,  0,        // f_pos 
	0,  3,  4,  4,  3,  0,
	0,  0,  5,  5,  0,  0,
	0,  0,  0,  0,  0,  0
  ,0
  },
	{ 0,  0,  0,  0,  0,  0,
	0,  0,  1,  0,  0,  0,
	0,  5,  4,  3,  2,  0,        //g_pos 
	0,  9,  8,  7,  6,  0,
	0,  0, 11, 10,  0,  0,
	0,  0,  0,  0,  0,  0
  ,0
  },
	{ 0,  1,  2,  2,  1,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,        // h_pos
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0
  ,0
  },
	{ 0,  0,  0,  0,  0,  0,
	0,  0,  0,  1,  0,  0,
	0,  2,  3,  4,  5,  0,        // i_pos 
	0,  2,  3,  4,  5,  0,
	0,  0,  0,  1,  0,  0,
	0,  0,  0,  0,  0,  0
  ,0
  },
	{ 0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,
	0,  1,  2,  2,  1,  0,        // j_pos 
	0,  1,  2,  2,  1,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0
  ,0
  },
	{ 0,  0,  0,  0,  0,  0,
	0,  0,  1,  0,  0,  0,
	0,  5,  4,  3,  2,  0,        // k_pos 
	0,  5,  4,  3,  2,  0,
	0,  0,  1,  0,  0,  0,
	0,  0,  0,  0,  0,  0
  ,0
  },
	{ 5,  4,  3,  2,  1,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,         // l_pos
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0
  ,0
  },
	{ 0,  0,  0,  0,  0,  0,
	0,  0, 10, 11,  0,  0,
	0,  6,  7,  8,  9,  0,        // m_pos 
	0,  2,  3,  4,  5,  0,
	0,  0,  0,  1,  0,  0,
	0,  0,  0,  0,  0,  0
  ,0
  },
	{ 0,  0,  0,  0,  0,  0,
	0,  0,  5,  5,  0,  0,
	0,  3,  4,  4,  3,  0,        // n_pos 
	0,  1,  2,  2,  1,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0
  ,0
  },
	{ 0,  0,  0,  0,  0,  0,
	0,  0, 11, 10,  0,  0,
	0,  9,  8,  7,  6,  0,        // o_pos 
	0,  5,  4,  3,  2,  0,
	0,  0,  1,  0,  0,  0,
	0,  0,  0,  0,  0,  0
  ,0
  },
};

int SymmetryPosition_EAIF[15] = {1,1,0,1,1,1,0,1,1,1,0,0,0,0,0};

int POS_EQUATION_NUMBER_EAIF[15] =  {
                     6, 3, 6,
									6,12, 6,12,
									3, 6, 3, 6,
									6,12, 6,12};


int FILTER_NEW_SUB_POS_EAIF[15] = {       
         a_pos, b_pos, a_pos,
	d_pos, e_pos, f_pos, e_pos,
	h_pos, i_pos, j_pos, i_pos,
	d_pos, e_pos, f_pos, e_pos };

#endif

#endif

