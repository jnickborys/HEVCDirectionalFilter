#ifdef EDAIF2
#include "aif_common.h"

static int TwoDEquationPattern_DAIF[MAX_NUM_AIF][SQR_FILTER]  =		// get equation number for symmetric filter coefficients
{
	{ 0,  1,  2,  3,  4,  5,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,        // a_pos
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0  },

	{ 0,  1,  2,  3,  4,  5,
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

	{ 0,  0,  0,  0,  0,  6,
	0,  1,  0,  0,  7,  0,
	0,  0,  2,  8,  0,  0,        // f_pos - diagonal 12-tap filter
	0,  0,  9,  3,  0,  0,
	0,  10,  0,  0,  4,  0,
	11,  0,  0,  0,  0,  5  },

	{ 0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  1,  0,
	0,  0,  0,  2,  0,  0,        //g_pos - diagonal SW-NE filter
	0,  0,  3,  0,  0,  0,
	0,  4,  0,  0,  0,  0,
	5,  0,  0,  0,  0,  0  },

	{ 0,  1,  2,  3,  4,  5,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,        // h_pos
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0  },

	{ 0,  0,  0,  0,  0,  11,
	0,  1,  0,  0,  10,  0,
	0,  0,  2,  9,  0,  0,        // i_pos - diagonal 12-tap filter
	0,  0,  8,  3,  0,  0,
	0,  7,  0,  0,  4,  0,
	6,  0,  0,  0,  0,  5  },

	{ 0,  0,  0,  0,  0,  6,
	0,  1,  0,  0,  7,  0,
	0,  0,  2,  8,  0,  0,        // j_pos - diagonal 12-tap filter
	0,  0,  9,  3,  0,  0,
	0,  10,  0,  0,  4,  0,
	11,  0,  0,  0,  0,  5  },

	{ 5,  0,  0,  0,  0,  6,
	0,  4,  0,  0,  7,  0,
	0,  0,  3,  8,  0,  0,        // k_pos - diagonal 12-tap filter
	0,  0,  9,  2,  0,  0,
	0,  10,  0,  0,  1,  0,
	11,  0,  0,  0,  0,  0  },

	{ 5,  4,  3,  2,  1,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,         // l_pos
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0  },

	{ 0,  0,  0,  0,  0,  5,
	0,  0,  0,  0,  4,  0,
	0,  0,  0,  3,  0,  0,        //m_pos - diagonal SW-NE filter
	0,  0,  2,  0,  0,  0,
	0,  1,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0  },
	{ 5,  0,  0,  0,  0,  11,
	0,  4,  0,  0,  10,  0,
	0,  0,  3,  9,  0,  0,        // n_pos - diagonal 12-tap filter
	0,  0,  8,  2,  0,  0,        
	0, 7,  0,  0,  1,  0,        
	6, 0,  0,  0,  0,  0	},
	{ 5,  0,  0,  0,  0,  0,
	  0,  4,  0,  0,  0,  0,
	  0,  0,  3,  0,  0,  0,        // o_pos - 
	  0,  0,  0,  2,  0,  0,        
	  0,  0,  0,  0,  1,  0,        
	  0,  0,  0,  0,  0,  0	},
};

static int Calc2Filt_Indexes_DAIF[MAX_NUM_AIF][SQR_FILTER]  =		// get one filter from another one, if symmetry properties are used (used e.g. in ExtendFilterCoefficients)
{
	{ 12, 13, 14, 15, 16, 17,
    0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0,  // a_pos 
		0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0  },

	{ 12, 13, 14, 15, 16, 17,
    0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0,        // b_pos
		0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0  },

	{ 17, 16, 15, 14, 13, 12,
    0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0,        // c_pos
		0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0  },

 	{ 2,  8, 14, 20, 26, 32,
    0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0,        // d_pos
		0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0  },

	{ 0, 7,  14,  21,  28,  35,
	  0, 0,  0,  0,  0,  0,
	  0, 0,  0,  0,  0,  0,        // e_pos
	  0, 0,  0,  0,  0,  0,
	  0, 0,  0,  0,  0,  0,
	  0, 0,  0,  0,  0,  0},

	{0,  7,  14,  21,  28, 35,
    5, 10,  15, 20, 25,  30,
	  0,  0,  0,  0,  0,  0,   			// f_pos
	  0,  0,  0,  0,  0,  0,
	  0,  0,  0,  0,  0,  0,
	  0,  0,  0,  0,  0,  0  },

	{ 5,  10,  15,  20,  25,  30,
		0,  0,  0,  0,  0,  0,
		0,  0,  0, 0,  0,  0,        // g_pos
		0,  0, 0,  0,   0,  0,
		0, 0,  0,  0,   0,  0,
		0, 0,  0,  0,   0,  0  },

	{ 2,  8, 14,  20,  26, 32,
		0,  0, 0,  0,  0,  0,
		0,  0, 0,  0,  0,  0,        // h_pos
		0,  0, 0,  0,  0,  0,
		0,  0, 0,  0,  0,  0,
		0,  0, 0,  0,  0,  0  },

	{ 0,  7,  14,  21,  28, 35,
    30, 25, 20,  15, 10,  5,
	  0,  0,  0,  0,   0,   0,   			// i_pos
	  0,  0,  0,  0,   0,   0,
	  0,  0,  0,  0,   0,   0,
	  0,  0,  0,  0,   0,   0  },

	{0,  7,  14,  21,  28, 35,
    5, 10,  15, 20, 25,  30,
	  0,  0,  0,  0,  0,  0,   			// j_pos
	  0,  0,  0,  0,  0,  0,
	  0,  0,  0,  0,  0,  0,
	  0,  0,  0,  0,  0,  0  },

	{35,  28, 21, 14,  7, 0,
    5, 10, 15, 20, 25, 30,
	  0,  0,  0,  0,  0,  0,   			// k_pos
	  0,  0,  0,  0,  0,  0,
	  0,  0,  0,  0,  0,  0,
	  0,  0,  0,  0,  0,  0  },

	{ 32,  26, 20,  14,  8,  2,
		0,  0, 0,  0,  0,  0,
		0,  0, 0,  0,  0,  0,        // l_pos
		0,  0, 0,  0,  0,  0,
		0,  0, 0,  0,  0,  0,
		0,  0, 0,  0,  0,  0  },

	{	30,	25,  20,  15,  10,  5,
		0,  0,  0,  0, 0,  0,
		0,  0,  0,  0,  0,  0,        // m_pos
		0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0,
		0,  0,  0,  0,  0,  0 },

	{ 35,  28,  21,  14,  7, 0,
    30, 25,  20,  15, 10,  5,
	  0,  0, 0, 0,  0,  0,   			// n_pos
	  0,  0, 0, 0,  0,  0,
	  0,  0,  0,  0,  0, 0,
	  0,  0,  0,  0,  0,  0},

  { 35,  28,  21,  14,  7, 0,
	  0, 0,  0,  0,  0,  0,
	  0, 0, 0,  0,  0,  0,        // o_pos
	  0, 0,  0, 0,  0,  0,
	  0, 0,  0,  0, 0,  0,
	  0, 0,  0,  0,  0, 0  },
/*
     0  1  2   3 4  5
    6  7  8   9 10 11
    12 13 14 15 16 17
    18 19 20 21 22 23
    24 25 26 27 28 29
    30 31 32 33 34 35*/
};

static double STANDARD_2D_FILTER_DAIF[MAX_NUM_AIF][SQR_FILTER] = 
{
{ 0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	3.0/128.0  ,   -15.0/128.0  ,   111.0/128.0  ,   37.0/128.0  ,   -10.0/128.0  ,   2.0/128.0  ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0
  ,0.0
}, // a_pos
{ 0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	3.0/128.0  ,   -17.0/128.0  ,   78.0/128.0  ,   78.0/128.0  ,   -17.0/128.0  ,   3.0/128.0  ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0
  ,0.0
}, // b_pos
{ 0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	2.0/128.0  ,   -10.0/128.0  ,   37.0/128.0  ,   111.0/128.0  ,   -15.0/128.0  ,   3.0/128.0  ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0       ,    0.0       ,    0.0       ,   0.0
  ,0.0
}, // c_pos
{ 0.0       ,    0.0       ,    3.0/128.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   -15.0/128.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   111.0/128.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   37.0/128.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   -10.0/128.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    2.0/128.0  ,    0.0       ,    0.0       ,   0.0
  ,0.0
}, // d_pos
{ 3.0/128.0    ,    0.0       ,    0.0		,    0.0       ,		0.0       ,   0.0       ,
	0.0       ,  -15.0/128.0   ,    0.0		,    0.0       ,		0.0       ,   0.0       ,
	0.0		  ,    0.0		 ,   111.0/128.0  ,   0.0		   ,	    0.0		  ,   0.0		,
	0.0       ,    0.0       ,    0.0		,    37.0/128.0 ,		0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0		,    0.0       ,    -10.0/128.0     ,   0.0       ,
	0.0       ,    0.0       ,    0.0		,    0.0       ,		0.0       ,   2.0/128.0
  ,0.0
	}, // e_pos
	
{ 3.0/256.0,			0.0	,		0.0		,		0.0		,		0.0		,		3.0/256.0,
	0.0		,   -15.0/256.0	,	 0.0		,		0.0		,   -15.0/256.0	,			0.0 ,
	0.0		, 0.0			,   111/256.0	,		111/256.0,		0.0		,			0.0 ,
	0.0		, 0.0			,   37/256.0	,		37/256.0,		0.0		,			0.0 ,
	0.0		,   -10.0/256.0	,		0.0		,		0.0		,	 -10.0/256.0	,			0.0 ,
  2.0/256.0	,   0.0 	   ,	0.0			,		0.0		,		0.0		,		2.0/256.0
  ,0.0
}, // f_pos
{   0.0       ,    0.0		,    0.0        ,    0.0       ,	  0.0         ,	3.0/128.0    ,
	0.0       ,    0.0       ,    0.0		,    0.0       ,	-15.0/128.0	  ,   0.0       ,
	0.0		  ,    0.0		 ,    0.0	    ,	111.0/128.0  ,		0.0		  ,   0.0		,
	0.0       ,    0.0       ,   37.0/128.0  ,	 0.0	   ,    	0.0       ,   0.0       ,
	0.0       ,    -10.0/128.0 ,    0.0		,    0.0       ,		0.0       ,   0.0       ,
	2.0/128.0  ,    0.0       ,    0.0		,    0.0       ,		0.0       ,   0.0       
  ,0.0
	}, // g_pos
{ 0.0       ,    0.0       ,    3.0/128.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   -17.0/128.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   78.0/128.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   78.0/128.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   -17.0/128.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    3.0/128.0  ,    0.0       ,    0.0       ,   0.0
  ,0.0
}, // h_pos
{	3.0/256   ,    0.0       ,    0.0	  	,    0.0       ,		0.0       ,   2.0/256.0 ,
	0.0       ,  -15.0/256.0 ,    0.0		,    0.0       ,	-10.0/256.0 ,   0.0       ,
	0.0		    ,    0.0		   ,   111.0/256.0,   37.0/256.0 ,	    0.0		  ,   0.0    		,
	0.0       ,    0.0       ,   111.0/256.0,   37.0/256.0 ,		0.0       ,   0.0       ,
	0.0       ,  -15.0/256.0 ,    0.0		,    0.0       ,   -10.0/256.0,   0.0       ,
	3.0/256	  ,    0.0       ,    0.0		  ,    0.0       ,		0.0       ,   2.0/256.0
  ,0.0
	},//new i_pos
{  3.0/256.0,			0.0	,		0.0		,		0.0		,		0.0		,		3.0/256.0,
	0.0		,   -17.0/256.0	,	 0.0		,		0.0		,   -17.0/256.0	,			0.0 ,
	0.0		, 0.0			,   78/256.0	,		78/256.0,		0.0		,			0.0 ,
	0.0		, 0.0			,   78/256.0	,		78/256.0,		0.0		,			0.0 ,
	0.0		,   -17.0/256.0	,		0.0		,		0.0		,	 -17.0/256.0	,			0.0 ,
  3.0/256.0	,   0.0 	   ,	0.0			,		0.0		,		0.0		,		3.0/256.0
  ,0.0
	},// j_pos
{  2.0/256.0,			0.0	,		0.0		,		0.0		,		0.0		,		3.0/256.0,
	0.0		,   -10.0/256.0	,	 0.0		,		0.0		,   -15.0/256.0	,			0.0 ,
	0.0		, 0.0			,   37/256.0	,		111/256.0,		0.0		,			0.0 ,
	0.0		, 0.0			,   37/256.0	,		111/256.0,		0.0		,			0.0 ,
	0.0		,   -10.0/256.0	,		0.0		,		0.0		,	 -15.0/256.0	,			0.0 ,
  2.0/256.0	,   0.0 	   ,	0.0			,		0.0		,		0.0		,		3.0/256.0
  ,0.0
	}, // k_pos
{ 0.0       ,    0.0       ,    2.0/128.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   -10.0/128.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   37.0/128.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   111.0/128.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,   -15.0/128.0  ,    0.0       ,    0.0       ,   0.0       ,
	0.0       ,    0.0       ,    3.0/128.0  ,    0.0       ,    0.0       ,   0.0
  ,0.0
}, // l_pos
{   0.0       ,    0.0		,    0.0        ,    0.0       ,	  0.0         ,	2.0/128.0    ,
	0.0       ,    0.0       ,    0.0		,    0.0       ,	-10.0/128.0	  ,   0.0       ,
	0.0		  ,    0.0		 ,    0.0	    ,	37.0/128.0  ,		0.0		  ,   0.0		,
	0.0       ,    0.0       ,    111.0/128.0 ,	 0.0	   ,    	0.0       ,   0.0       ,
	0.0       ,    -15.0/128.0 ,    0.0		,    0.0       ,		0.0       ,   0.0       ,
	3.0/128.0  ,    0.0       ,    0.0		,    0.0       ,		0.0       ,   0.0       
  ,0.0
	}, // m_pos
{  2.0/256.0,			0.0	,		0.0		,		0.0		,		0.0		,		2.0/256.0,
	0.0		,   -10.0/256.0	,	 0.0		,		0.0		,   -10.0/256.0	,			0.0 ,
	0.0		, 0.0			,   37/256.0	,		37/256.0,		0.0		,			0.0 ,
	0.0		, 0.0			,   111/256.0	,		111/256.0,		0.0		,			0.0 ,
	0.0		,   -15.0/256.0	,		0.0		,		0.0		,	 -15.0/256.0	,			0.0 ,
  3.0/256.0	,   0.0 	   ,	0.0			,		0.0		,		0.0		,		3.0/256.0
  ,0.0
	}, // n_pos
{ 2.0/128.0    ,    0.0       ,    0.0		,    0.0       ,		0.0       ,   0.0       ,
	0.0       ,  -10.0/128.0   ,    0.0		,    0.0       ,		0.0       ,   0.0       ,
	0.0		  ,    0.0		 ,    37.0/128.0 ,    0.0	   ,	    0.0		  ,   0.0		,
	0.0       ,    0.0       ,    0.0		,    111.0/128.0 ,		0.0       ,   0.0       ,
	0.0       ,    0.0       ,    0.0		,    0.0       ,    -15.0/128.0     ,   0.0       ,
	0.0       ,    0.0       ,    0.0		,    0.0       ,		0.0       ,   3.0/128.0
  ,0.0
	}, // o_pos
};

static int Filt_Indexes_offset_DAIF[MAX_NUM_AIF] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; // how many times we copu data from CalcFilters.

static int SymmetryPosition_DAIF[MAX_NUM_AIF] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};	// if 0, the position is copied from another one
static int Is1DPosition_DAIF[MAX_NUM_AIF] =	  {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; // if 0, the position is a 2D one
static int IsDiagonal1D_DAIF[MAX_NUM_AIF] =	  {0,0,0,0,1,3,2,0,3,3,3,0,2,3,1}; // if 1, the filter alighned NW-SE, 2 - NE-SW, 3 - diagonal cross
static int nBitsIntRepresent_DIAF[MAX_NUM_AIF] =	  {7,7,7,7,7,8,7,7,8,8,8,7,7,8,7}; // if 1, the filter alighned NW-SE, 2 - NE-SW, 3 - diagonal cross

static int POS_EQUATION_NUMBER_DAIF[MAX_NUM_AIF] =  {    
                      6, 6,   6,
									6,  6, 12,  6,
									6, 12, 12, 12,
									6,  6, 12,  6,};
static int FILTER_NEW_SUB_POS_DAIF[MAX_NUM_AIF] = {       
		   a_pos, b_pos, c_pos,
	d_pos, e_pos, f_pos, g_pos,
	h_pos, i_pos, j_pos, k_pos,
	l_pos, m_pos, n_pos, o_pos};

#endif
