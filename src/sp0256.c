// license:BSD-3-Clause
// copyright-holders:Joseph Zbiciak,Tim Lindner
// modified for WORM by Martin Howse

/**********************************************************************

    SP0256 Narrator Speech Processor emulation

**********************************************************************
                            _____   _____
                   Vss   1 |*    \_/     | 28  OSC 2
                _RESET   2 |             | 27  OSC 1
           ROM DISABLE   3 |             | 26  ROM CLOCK
                    C1   4 |             | 25  _SBY RESET
                    C2   5 |             | 24  DIGITAL OUT
                    C3   6 |             | 23  Vdi
                   Vdd   7 |    SP0256   | 22  TEST
                   SBY   8 |             | 21  SER IN
                  _LRQ   9 |             | 20  _ALD
                    A8  10 |             | 19  SE
                    A7  11 |             | 18  A1
               SER OUT  12 |             | 17  A2
                    A6  13 |             | 16  A3
                    A5  14 |_____________| 15  A4

**********************************************************************/

/*
   GI SP0256 Narrator Speech Processor

   By Joe Zbiciak. Ported to MESS by tim lindner.

*/


#include "sp0256.h"
#include "audio.h"
#include "english2phoneme/TTS.h"
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "sp0romstest.h" // this has all roms as: m_rom12, m_rom19, m_romAL2 // latter is with vocabs and TTS
#include "sp0256vocab.h"
#include "resources.h"
extern float _selx, _sely, _selz, _mode, _speed;
//extern u8 TTS;
static const unsigned char *m_romm;
static u8 modus=0;

typedef unsigned char UINT8;
typedef signed char INT8;
typedef u16 UINT16;
typedef int16_t INT16;
typedef uint32_t UINT32;
typedef int32_t INT32;

void bitrevbuff(UINT8 *buffer, unsigned int start, unsigned int length);

struct lpc12_t
{
  INT16     rpt, cnt;       /* Repeat counter, Period down-counter.         */
  UINT32  per, per_orig, rng;       /* Period, Amplitude, Random Number Generator   */
	INT16     amp;
	INT16   f_coef[6];      /* F0 through F5.                               */
	INT16   b_coef[6];      /* B0 through B5.                               */
	INT16   z_data[6][2];   /* Time-delay data for the filter stages.       */
	UINT8   r[16];          /* The encoded register set.                    */
	INT16     interp;
};

u8           m_silent;          /* Flag: SP0256 is silent.                      */


static struct lpc12_t m_filt;            /* 12-pole filter                               */
static uint16_t            m_lrq;             /* Load ReQuest.  == 0 if we can accept a load  */
static int32_t            m_ald;             /* Address LoaD.  < 0 if no command pending.    */
static int32_t            m_pc;              /* Microcontroller's PC value.                  */
static uint32_t            m_stack;           /* Microcontroller's PC stack.                  */
static u8            m_fifo_sel;        /* True when executing from FIFO.               */
static u8            m_halted;          /* True when CPU is halted.                     */
//static u8 really=0;
static UINT32         m_mode;            /* Mode register.                               */
static UINT32         m_page;            /* Page set by SETPAGE                          */

static UINT32         m_fifo_head;       /* FIFO head pointer (where new data goes).     */
static UINT32         m_fifo_tail;       /* FIFO tail pointer (where data comes from).   */
static UINT32         m_fifo_bitp;       /* FIFO bit-pointer (for partial decles).       */
static UINT16         m_fifo[64];        /* The 64-decle FIFO.                           */


#define CLOCK_DIVIDER (7*6*8)
#define HIGH_QUALITY

#define PER_PAUSE    (64)               /* Equiv timing period for pauses.  */
#define PER_NOISE    (64)               /* Equiv timing period for noise.   */

#define FIFO_ADDR    (0x1800 << 3)      /* SP0256 address of SPB260 speech FIFO. = 49152 decimal  */


/* ======================================================================== */
/*  qtbl  -- Coefficient Quantization Table.  This comes from a             */
/*              SP0250 data sheet, and should be correct for SP0256.        */
/* ======================================================================== */
const INT16 qtbl[128]     =
{
	0,      9,      17,     25,     33,     41,     49,     57,
	65,     73,     81,     89,     97,     105,    113,    121,
	129,    137,    145,    153,    161,    169,    177,    185,
	193,    201,    209,    217,    225,    233,    241,    249,
	257,    265,    273,    281,    289,    297,    301,    305,
	309,    313,    317,    321,    325,    329,    333,    337,
	341,    345,    349,    353,    357,    361,    365,    369,
	373,    377,    381,    385,    389,    393,    397,    401,
	405,    409,    413,    417,    421,    425,    427,    429,
	431,    433,    435,    437,    439,    441,    443,    445,
	447,    449,    451,    453,    455,    457,    459,    461,
	463,    465,    467,    469,    471,    473,    475,    477,
	479,    481,    482,    483,    484,    485,    486,    487,
	488,    489,    490,    491,    492,    493,    494,    495,
	496,    497,    498,    499,    500,    501,    502,    503,
	504,    505,    506,    507,    508,    509,    510,    511
};

static void sp0256_iinit()
{

	/* -------------------------------------------------------------------- */
	/*  Configure our internal variables.                                   */
	/* -------------------------------------------------------------------- */
	m_filt.rng = 1;

	/* -------------------------------------------------------------------- */
	/*  Set up the microsequencer's initial state.                          */
	/* -------------------------------------------------------------------- */

	m_halted   = 1; // was 1
	m_filt.rpt = -1;
	m_lrq      = 0x8000;
 	m_page     = 0x1000 << 3; //32768 =0x8000
	m_silent   = 1;

	/* -------------------------------------------------------------------- */
	/*  Setup the ROM.                                                      */
	/* -------------------------------------------------------------------- */
	// the rom is not supposed to be reversed first; according to Joe Zbiciak.
	// see http://forums.bannister.org/ubbthreads.php?ubb=showflat&Number=72385#Post72385
	// TODO: because of this, check if the bitrev functions are even used anywhere else
	//	bitrevbuff(m_rom, 0, 0xffff);
	sp0256_newsay();
}


//-------------------------------------------------
//  device_reset - device-specific reset
//-------------------------------------------------

static void reset()
{
	// reset FIFO and SP0256
	m_fifo_head = m_fifo_tail = m_fifo_bitp = 0;

	memset(&m_filt, 0, sizeof(m_filt));
	m_halted   = 1;
	m_filt.rpt = -1;
	m_filt.rng = 1;
	m_lrq      = 0x8000;
	m_ald      = 0x0000;
	m_pc       = 0x0000;
	m_stack    = 0x0000;
	m_fifo_sel = 0;
	m_mode     = 0;
	m_page     = 0x1000 << 3;
	m_silent   = 1;
	//	m_sby_line = 0;

	  m_lrq = 0;
}



/* ======================================================================== */
/*  LIMIT            -- Limiter function for digital sample output.         */
/* ======================================================================== */
static inline INT16 limit(INT16 s)
{
#ifdef HIGH_QUALITY /* Higher quality than the original, but who cares? */
	if (s >  8191) return  8191;
	if (s < -8192) return -8192;
#else
	if (s >  127) return  127;
	if (s < -128) return -128;
#endif
	return s;
}

/* ======================================================================== */
/*  LPC12_UPDATE     -- Update the 12-pole filter, outputting samples.      */
/* ======================================================================== */

static inline void lpc12_update(struct lpc12_t *f, INT16* out)
{
	u8 j;
	INT16 samp;
	u8 do_int;
	int16_t val;
	//	static int count=128;
	/* -------------------------------------------------------------------- */
	/*  Iterate up to the desired number of samples.  We actually may       */
	/*  break out early if our repeat count expires.                        */
	/* -------------------------------------------------------------------- */
		/* ---------------------------------------------------------------- */
		/*  Generate a series of periodic impulses, or random noise.        */
		/* ---------------------------------------------------------------- */
		do_int = 0;
		//		samp   = 0;
		//		f->per_orig=0;
		//		f->amp=512;
		//		f->per_orig=128;
		if (f->per_orig)
		{
		  //		  val=_selx*1029.0f;
		  //		  MAXED(val,1023);
		  if (modus!=1)	f->per=f->per_orig*logspeed[val];
		  		  else f->per=64.0f;//*logspeed[val];
		  //		  f->per=f->per_orig;
		  //		  		  f->per=f->per_orig;
		  //		  f->per=count++;
		  //		  if (count>255) count=128;
			if (f->cnt <= 0)
			{
				f->cnt += f->per;
				samp    = f->amp;
				f->rpt--;
				do_int  = f->interp;

				for (j = 0; j < 6; j++)
				  f->z_data[j][1] = f->z_data[j][0] = 0;

			} else
			{
				samp = 0;
				f->cnt--;
			}

		} else
		{
			u8 bit;

			if (--f->cnt <= 0)
			{
				do_int = f->interp;
				f->cnt = PER_NOISE;
				f->rpt--;
				for (j = 0; j < 6; j++)
				  f->z_data[j][0] = f->z_data[j][1] = 0;
			}

			bit = f->rng & 1;
			f->rng = (f->rng >> 1) ^ (bit ? 0x4001 : 0);
			//bit=rand()%2;
			if (bit) { samp =  f->amp; }
			else     { samp = -f->amp; }
			//			samp=0;
		}

		/* ---------------------------------------------------------------- */
		/*  If we need to, process the interpolation registers.             */
		/* ---------------------------------------------------------------- */
		if (do_int)
		{
			f->r[0] += f->r[14];
			f->r[1] += f->r[15];

			f->amp   = (f->r[0] & 0x1F) << (((f->r[0] & 0xE0) >> 5) + 0);
			//			f->amp=160;
			f->per_orig   = f->r[1];

			do_int   = 0;
		}

		/* ---------------------------------------------------------------- */
		/*  Stop if we expire our repeat counter and return the actual      */
		/*  number of samples we did.                                       */
		/* ---------------------------------------------------------------- */
				if (f->rpt <= 0) {
		  		  *out=0;
		  		  return; 
				}
		/* ---------------------------------------------------------------- */
		/*  Each 2nd order stage looks like one of these.  The App. Manual  */
		/*  gives the first form, the patent gives the second form.         */
		/*  They're equivalent except for time delay.  I implement the      */
		/*  first form.   (Note: 1/Z == 1 unit of time delay.)              */
		/*                                                                  */
		/*          ---->(+)-------->(+)----------+------->                 */
		/*                ^           ^           |                         */
		/*                |           |           |                         */
		/*                |           |           |                         */
		/*               [B]        [2*F]         |                         */
		/*                ^           ^           |                         */
		/*                |           |           |                         */
		/*                |           |           |                         */
		/*                +---[1/Z]<--+---[1/Z]<--+                         */
		/*                                                                  */
		/*                                                                  */
		/*                +---[2*F]<---+                                    */
		/*                |            |                                    */
		/*                |            |                                    */
		/*                v            |                                    */
		/*          ---->(+)-->[1/Z]-->+-->[1/Z]---+------>                 */
		/*                ^                        |                        */
		/*                |                        |                        */
		/*                |                        |                        */
		/*                +-----------[B]<---------+                        */
		/*                                                                  */
		/* ---------------------------------------------------------------- */

				/*				for (j = 0; j < 6; j++)
		{
		  samp += (((int32_t)f->b_coef[j] * (int32_t)f->z_data[j][1]) >> 9);
		  samp += (((int32_t)f->f_coef[j] * (int32_t)f->z_data[j][0]) >> 8);

			f->z_data[j][1] = f->z_data[j][0];
			f->z_data[j][0] = samp;
		}
				*/
						//						*out= limit(samp);
				samp=samp<<8;
				if (samp>32767) samp=32767;
				else if (samp<-32767) samp=-32767;
				*out=samp;
				
		//		return 1;
}

static u8 stage_map[6] = { 0, 1, 2, 3, 4, 5 };

/* ======================================================================== */
/*  LPC12_REGDEC -- Decode the register set in the filter bank.             */
/* ======================================================================== */
static inline void lpc12_regdec(struct lpc12_t *f)
{
	u8 i;

	/* -------------------------------------------------------------------- */
	/*  Decode the Amplitude and Period registers.  Force the 'cnt' to 0    */
	/*  to get an initial impulse.  We compensate elsewhere by setting      */
	/*  the repeat count to "repeat + 1".                                   */
	/* -------------------------------------------------------------------- */
	f->amp = (f->r[0] & 0x1F) << (((f->r[0] & 0xE0) >> 5) + 0);
	//	f->amp=160;			
	//	if (modus==3 || modus==4) f->cnt = 0; // why just 3 and 4 - non single phon modes????
	f->cnt = 0;
	f->per_orig = f->r[1];

	/* -------------------------------------------------------------------- */
	/*  Decode the filter coefficients from the quant table.                */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 6; i++)
	{
#define IQ(x) (((x) & 0x80) ? qtbl[0x7F & -(x)] : -qtbl[(x)])

		f->b_coef[stage_map[i]] = IQ(f->r[2 + 2*i]);
		f->f_coef[stage_map[i]] = IQ(f->r[3 + 2*i]);
	}

	/* -------------------------------------------------------------------- */
	/*  Set the Interp flag based on whether we have interpolation parms    */
	/* -------------------------------------------------------------------- */
	f->interp = f->r[14] || f->r[15];

	return;
}

/* ======================================================================== */
/*  SP0256_DATAFMT   -- Data format table for the SP0256's microsequencer   */
/*                                                                          */
/*  len     4 bits      Length of field to extract                          */
/*  lshift  4 bits      Left-shift amount on field                          */
/*  param   4 bits      Parameter number being updated                      */
/*  delta   1 bit       This is a delta-update.  (Implies sign-extend)      */
/*  field   1 bit       This is a field replace.                            */
/*  clr5    1 bit       Clear F5, B5.                                       */
/*  clrall  1 bit       Clear all before doing this update                  */
/* ======================================================================== */

#define CR(l,s,p,d,f,c5,ca)         \
		(                           \
			(((l)  & 15) <<  0) |   \
			(((s)  & 15) <<  4) |   \
			(((p)  & 15) <<  8) |   \
			(((d)  &  1) << 12) |   \
			(((f)  &  1) << 13) |   \
			(((c5) &  1) << 14) |   \
			(((ca) &  1) << 15)     \
		)

#define CR_DELTA  CR(0,0,0,1,0,0,0)
#define CR_FIELD  CR(0,0,0,0,1,0,0)
#define CR_CLR5   CR(0,0,0,0,0,1,0)
#define CR_CLRA   CR(0,0,0,0,0,0,1)
#define CR_LEN(x) ((x) & 15)
#define CR_SHF(x) (((x) >> 4) & 15)
#define CR_PRM(x) (((x) >> 8) & 15)

enum { AM = 0, PR, B0, F0, B1, F1, B2, F2, B3, F3, B4, F4, B5, F5, IA, IP };

const UINT16 sp0256_datafmt[]     =
{
	/* -------------------------------------------------------------------- */
	/*  OPCODE 1111: PAUSE                                                  */
	/* -------------------------------------------------------------------- */
	/*    0 */  CR( 0,  0,  0,  0,  0,  0,  1),     /*  Clear all   */

	/* -------------------------------------------------------------------- */
	/*  Opcode 0001: LOADALL                                                */
	/* -------------------------------------------------------------------- */
				/* All modes                */
	/*    1 */  CR( 8,  0,  AM, 0,  0,  0,  1),     /*  Amplitude   */
	/*    2 */  CR( 8,  0,  PR, 0,  0,  0,  0),     /*  Period      */
	/*    3 */  CR( 8,  0,  B0, 0,  0,  0,  0),     /*  B0          */
	/*    4 */  CR( 8,  0,  F0, 0,  0,  0,  0),     /*  F0          */
	/*    5 */  CR( 8,  0,  B1, 0,  0,  0,  0),     /*  B1          */
	/*    6 */  CR( 8,  0,  F1, 0,  0,  0,  0),     /*  F1          */
	/*    7 */  CR( 8,  0,  B2, 0,  0,  0,  0),     /*  B2          */
	/*    8 */  CR( 8,  0,  F2, 0,  0,  0,  0),     /*  F2          */
	/*    9 */  CR( 8,  0,  B3, 0,  0,  0,  0),     /*  B3          */
	/*   10 */  CR( 8,  0,  F3, 0,  0,  0,  0),     /*  F3          */
	/*   11 */  CR( 8,  0,  B4, 0,  0,  0,  0),     /*  B4          */
	/*   12 */  CR( 8,  0,  F4, 0,  0,  0,  0),     /*  F4          */
	/*   13 */  CR( 8,  0,  B5, 0,  0,  0,  0),     /*  B5          */
	/*   14 */  CR( 8,  0,  F5, 0,  0,  0,  0),     /*  F5          */
				/* Mode 01 and 11 only      */
	/*   15 */  CR( 8,  0,  IA, 0,  0,  0,  0),     /*  Amp Interp  */
	/*   16 */  CR( 8,  0,  IP, 0,  0,  0,  0),     /*  Pit Interp  */

	/* -------------------------------------------------------------------- */
	/*  Opcode 0100: LOAD_4                                                 */
	/* -------------------------------------------------------------------- */
				/* Mode 00 and 01           */
	/*   17 */  CR( 6,  2,  AM, 0,  0,  0,  1),     /*  Amplitude   */
	/*   18 */  CR( 8,  0,  PR, 0,  0,  0,  0),     /*  Period      */
	/*   19 */  CR( 4,  3,  B3, 0,  0,  0,  0),     /*  B3 (S=0)    */
	/*   20 */  CR( 6,  2,  F3, 0,  0,  0,  0),     /*  F3          */
	/*   21 */  CR( 7,  1,  B4, 0,  0,  0,  0),     /*  B4          */
	/*   22 */  CR( 6,  2,  F4, 0,  0,  0,  0),     /*  F4          */
				/* Mode 01 only             */
	/*   23 */  CR( 8,  0,  B5, 0,  0,  0,  0),     /*  B5          */
	/*   24 */  CR( 8,  0,  F5, 0,  0,  0,  0),     /*  F5          */

				/* Mode 10 and 11           */
	/*   25 */  CR( 6,  2,  AM, 0,  0,  0,  1),     /*  Amplitude   */
	/*   26 */  CR( 8,  0,  PR, 0,  0,  0,  0),     /*  Period      */
	/*   27 */  CR( 6,  1,  B3, 0,  0,  0,  0),     /*  B3 (S=0)    */
	/*   28 */  CR( 7,  1,  F3, 0,  0,  0,  0),     /*  F3          */
	/*   29 */  CR( 8,  0,  B4, 0,  0,  0,  0),     /*  B4          */
	/*   30 */  CR( 8,  0,  F4, 0,  0,  0,  0),     /*  F4          */
				/* Mode 11 only             */
	/*   31 */  CR( 8,  0,  B5, 0,  0,  0,  0),     /*  B5          */
	/*   32 */  CR( 8,  0,  F5, 0,  0,  0,  0),     /*  F5          */

	/* -------------------------------------------------------------------- */
	/*  Opcode 0110: SETMSB_6                                               */
	/* -------------------------------------------------------------------- */
				/* Mode 00 only             */
	/*   33 */  CR( 0,  0,  0,  0,  0,  1,  0),     /*  Clear 5     */
				/* Mode 00 and 01           */
	/*   34 */  CR( 6,  2,  AM, 0,  0,  0,  0),     /*  Amplitude   */
	/*   35 */  CR( 6,  2,  F3, 0,  1,  0,  0),     /*  F3 (5 MSBs) */
	/*   36 */  CR( 6,  2,  F4, 0,  1,  0,  0),     /*  F4 (5 MSBs) */
				/* Mode 01 only             */
	/*   37 */  CR( 8,  0,  F5, 0,  1,  0,  0),     /*  F5 (5 MSBs) */

				/* Mode 10 only             */
	/*   38 */  CR( 0,  0,  0,  0,  0,  1,  0),     /*  Clear 5     */
				/* Mode 10 and 11           */
	/*   39 */  CR( 6,  2,  AM, 0,  0,  0,  0),     /*  Amplitude   */
	/*   40 */  CR( 7,  1,  F3, 0,  1,  0,  0),     /*  F3 (6 MSBs) */
	/*   41 */  CR( 8,  0,  F4, 0,  1,  0,  0),     /*  F4 (6 MSBs) */
				/* Mode 11 only             */
	/*   42 */  CR( 8,  0,  F5, 0,  1,  0,  0),     /*  F5 (6 MSBs) */

	/*   43 */  0,  /* unused */
	/*   44 */  0,  /* unused */

	/* -------------------------------------------------------------------- */
	/*  Opcode 1001: DELTA_9                                                */
	/* -------------------------------------------------------------------- */
				/* Mode 00 and 01           */
	/*   45 */  CR( 4,  2,  AM, 1,  0,  0,  0),     /*  Amplitude   */
	/*   46 */  CR( 5,  0,  PR, 1,  0,  0,  0),     /*  Period      */
	/*   47 */  CR( 3,  4,  B0, 1,  0,  0,  0),     /*  B0 4 MSBs   */
	/*   48 */  CR( 3,  3,  F0, 1,  0,  0,  0),     /*  F0 5 MSBs   */
	/*   49 */  CR( 3,  4,  B1, 1,  0,  0,  0),     /*  B1 4 MSBs   */
	/*   50 */  CR( 3,  3,  F1, 1,  0,  0,  0),     /*  F1 5 MSBs   */
	/*   51 */  CR( 3,  4,  B2, 1,  0,  0,  0),     /*  B2 4 MSBs   */
	/*   52 */  CR( 3,  3,  F2, 1,  0,  0,  0),     /*  F2 5 MSBs   */
	/*   53 */  CR( 3,  3,  B3, 1,  0,  0,  0),     /*  B3 5 MSBs   */
	/*   54 */  CR( 4,  2,  F3, 1,  0,  0,  0),     /*  F3 6 MSBs   */
	/*   55 */  CR( 4,  1,  B4, 1,  0,  0,  0),     /*  B4 7 MSBs   */
	/*   56 */  CR( 4,  2,  F4, 1,  0,  0,  0),     /*  F4 6 MSBs   */
				/* Mode 01 only             */
	/*   57 */  CR( 5,  0,  B5, 1,  0,  0,  0),     /*  B5 8 MSBs   */
	/*   58 */  CR( 5,  0,  F5, 1,  0,  0,  0),     /*  F5 8 MSBs   */

				/* Mode 10 and 11           */
	/*   59 */  CR( 4,  2,  AM, 1,  0,  0,  0),     /*  Amplitude   */
	/*   60 */  CR( 5,  0,  PR, 1,  0,  0,  0),     /*  Period      */
	/*   61 */  CR( 4,  1,  B0, 1,  0,  0,  0),     /*  B0 7 MSBs   */
	/*   62 */  CR( 4,  2,  F0, 1,  0,  0,  0),     /*  F0 6 MSBs   */
	/*   63 */  CR( 4,  1,  B1, 1,  0,  0,  0),     /*  B1 7 MSBs   */
	/*   64 */  CR( 4,  2,  F1, 1,  0,  0,  0),     /*  F1 6 MSBs   */
	/*   65 */  CR( 4,  1,  B2, 1,  0,  0,  0),     /*  B2 7 MSBs   */
	/*   66 */  CR( 4,  2,  F2, 1,  0,  0,  0),     /*  F2 6 MSBs   */
	/*   67 */  CR( 4,  1,  B3, 1,  0,  0,  0),     /*  B3 7 MSBs   */
	/*   68 */  CR( 5,  1,  F3, 1,  0,  0,  0),     /*  F3 7 MSBs   */
	/*   69 */  CR( 5,  0,  B4, 1,  0,  0,  0),     /*  B4 8 MSBs   */
	/*   70 */  CR( 5,  0,  F4, 1,  0,  0,  0),     /*  F4 8 MSBs   */
				/* Mode 11 only             */
	/*   71 */  CR( 5,  0,  B5, 1,  0,  0,  0),     /*  B5 8 MSBs   */
	/*   72 */  CR( 5,  0,  F5, 1,  0,  0,  0),     /*  F5 8 MSBs   */

	/* -------------------------------------------------------------------- */
	/*  Opcode 1010: SETMSB_A                                               */
	/* -------------------------------------------------------------------- */
				/* Mode 00 only             */
	/*   73 */  CR( 0,  0,  0,  0,  0,  1,  0),     /*  Clear 5     */
				/* Mode 00 and 01           */
	/*   74 */  CR( 6,  2,  AM, 0,  0,  0,  0),     /*  Amplitude   */
	/*   75 */  CR( 5,  3,  F0, 0,  1,  0,  0),     /*  F0 (5 MSBs) */
	/*   76 */  CR( 5,  3,  F1, 0,  1,  0,  0),     /*  F1 (5 MSBs) */
	/*   77 */  CR( 5,  3,  F2, 0,  1,  0,  0),     /*  F2 (5 MSBs) */

				/* Mode 10 only             */
	/*   78 */  CR( 0,  0,  0,  0,  0,  1,  0),     /*  Clear 5     */
				/* Mode 10 and 11           */
	/*   79 */  CR( 6,  2,  AM, 0,  0,  0,  0),     /*  Amplitude   */
	/*   80 */  CR( 6,  2,  F0, 0,  1,  0,  0),     /*  F0 (6 MSBs) */
	/*   81 */  CR( 6,  2,  F1, 0,  1,  0,  0),     /*  F1 (6 MSBs) */
	/*   82 */  CR( 6,  2,  F2, 0,  1,  0,  0),     /*  F2 (6 MSBs) */

	/* -------------------------------------------------------------------- */
	/*  Opcode 0010: LOAD_2  Mode 00 and 10                                 */
	/*  Opcode 1100: LOAD_C  Mode 00 and 10                                 */
	/* -------------------------------------------------------------------- */
				/* LOAD_2, LOAD_C  Mode 00  */
	/*   83 */  CR( 6,  2,  AM, 0,  0,  0,  1),     /*  Amplitude   */
	/*   84 */  CR( 8,  0,  PR, 0,  0,  0,  0),     /*  Period      */
	/*   85 */  CR( 3,  4,  B0, 0,  0,  0,  0),     /*  B0 (S=0)    */
	/*   86 */  CR( 5,  3,  F0, 0,  0,  0,  0),     /*  F0          */
	/*   87 */  CR( 3,  4,  B1, 0,  0,  0,  0),     /*  B1 (S=0)    */
	/*   88 */  CR( 5,  3,  F1, 0,  0,  0,  0),     /*  F1          */
	/*   89 */  CR( 3,  4,  B2, 0,  0,  0,  0),     /*  B2 (S=0)    */
	/*   90 */  CR( 5,  3,  F2, 0,  0,  0,  0),     /*  F2          */
	/*   91 */  CR( 4,  3,  B3, 0,  0,  0,  0),     /*  B3 (S=0)    */
	/*   92 */  CR( 6,  2,  F3, 0,  0,  0,  0),     /*  F3          */
	/*   93 */  CR( 7,  1,  B4, 0,  0,  0,  0),     /*  B4          */
	/*   94 */  CR( 6,  2,  F4, 0,  0,  0,  0),     /*  F4          */
				/* LOAD_2 only              */
	/*   95 */  CR( 5,  0,  IA, 0,  0,  0,  0),     /*  Ampl. Intr. */
	/*   96 */  CR( 5,  0,  IP, 0,  0,  0,  0),     /*  Per. Intr.  */

				/* LOAD_2, LOAD_C  Mode 10  */
	/*   97 */  CR( 6,  2,  AM, 0,  0,  0,  1),     /*  Amplitude   */
	/*   98 */  CR( 8,  0,  PR, 0,  0,  0,  0),     /*  Period      */
	/*   99 */  CR( 6,  1,  B0, 0,  0,  0,  0),     /*  B0 (S=0)    */
	/*  100 */  CR( 6,  2,  F0, 0,  0,  0,  0),     /*  F0          */
	/*  101 */  CR( 6,  1,  B1, 0,  0,  0,  0),     /*  B1 (S=0)    */
	/*  102 */  CR( 6,  2,  F1, 0,  0,  0,  0),     /*  F1          */
	/*  103 */  CR( 6,  1,  B2, 0,  0,  0,  0),     /*  B2 (S=0)    */
	/*  104 */  CR( 6,  2,  F2, 0,  0,  0,  0),     /*  F2          */
	/*  105 */  CR( 6,  1,  B3, 0,  0,  0,  0),     /*  B3 (S=0)    */
	/*  106 */  CR( 7,  1,  F3, 0,  0,  0,  0),     /*  F3          */
	/*  107 */  CR( 8,  0,  B4, 0,  0,  0,  0),     /*  B4          */
	/*  108 */  CR( 8,  0,  F4, 0,  0,  0,  0),     /*  F4          */
				/* LOAD_2 only              */
	/*  109 */  CR( 5,  0,  IA, 0,  0,  0,  0),     /*  Ampl. Intr. */
	/*  110 */  CR( 5,  0,  IP, 0,  0,  0,  0),     /*  Per. Intr.  */

	/* -------------------------------------------------------------------- */
	/*  OPCODE 1101: DELTA_D                                                */
	/* -------------------------------------------------------------------- */
				/* Mode 00 and 01           */
	/*  111 */  CR( 4,  2,  AM, 1,  0,  0,  0),     /*  Amplitude   */
	/*  112 */  CR( 5,  0,  PR, 1,  0,  0,  0),     /*  Period      */
	/*  113 */  CR( 3,  3,  B3, 1,  0,  0,  0),     /*  B3 5 MSBs   */
	/*  114 */  CR( 4,  2,  F3, 1,  0,  0,  0),     /*  F3 6 MSBs   */
	/*  115 */  CR( 4,  1,  B4, 1,  0,  0,  0),     /*  B4 7 MSBs   */
	/*  116 */  CR( 4,  2,  F4, 1,  0,  0,  0),     /*  F4 6 MSBs   */
				/* Mode 01 only             */
	/*  117 */  CR( 5,  0,  B5, 1,  0,  0,  0),     /*  B5 8 MSBs   */
	/*  118 */  CR( 5,  0,  F5, 1,  0,  0,  0),     /*  F5 8 MSBs   */

				/* Mode 10 and 11           */
	/*  119 */  CR( 4,  2,  AM, 1,  0,  0,  0),     /*  Amplitude   */
	/*  120 */  CR( 5,  0,  PR, 1,  0,  0,  0),     /*  Period      */
	/*  121 */  CR( 4,  1,  B3, 1,  0,  0,  0),     /*  B3 7 MSBs   */
	/*  122 */  CR( 5,  1,  F3, 1,  0,  0,  0),     /*  F3 7 MSBs   */
	/*  123 */  CR( 5,  0,  B4, 1,  0,  0,  0),     /*  B4 8 MSBs   */
	/*  124 */  CR( 5,  0,  F4, 1,  0,  0,  0),     /*  F4 8 MSBs   */
				/* Mode 11 only             */
	/*  125 */  CR( 5,  0,  B5, 1,  0,  0,  0),     /*  B5 8 MSBs   */
	/*  126 */  CR( 5,  0,  F5, 1,  0,  0,  0),     /*  F5 8 MSBs   */

	/* -------------------------------------------------------------------- */
	/*  OPCODE 1110: LOAD_E                                                 */
	/* -------------------------------------------------------------------- */
	/*  127 */  CR( 6,  2,  AM, 0,  0,  0,  0),     /*  Amplitude   */
	/*  128 */  CR( 8,  0,  PR, 0,  0,  0,  0),     /*  Period      */

	/* -------------------------------------------------------------------- */
	/*  Opcode 0010: LOAD_2  Mode 01 and 11                                 */
	/*  Opcode 1100: LOAD_C  Mode 01 and 11                                 */
	/* -------------------------------------------------------------------- */
				/* LOAD_2, LOAD_C  Mode 01  */
	/*  129 */  CR( 6,  2,  AM, 0,  0,  0,  1),     /*  Amplitude   */
	/*  130 */  CR( 8,  0,  PR, 0,  0,  0,  0),     /*  Period      */
	/*  131 */  CR( 3,  4,  B0, 0,  0,  0,  0),     /*  B0 (S=0)    */
	/*  132 */  CR( 5,  3,  F0, 0,  0,  0,  0),     /*  F0          */
	/*  133 */  CR( 3,  4,  B1, 0,  0,  0,  0),     /*  B1 (S=0)    */
	/*  134 */  CR( 5,  3,  F1, 0,  0,  0,  0),     /*  F1          */
	/*  135 */  CR( 3,  4,  B2, 0,  0,  0,  0),     /*  B2 (S=0)    */
	/*  136 */  CR( 5,  3,  F2, 0,  0,  0,  0),     /*  F2          */
	/*  137 */  CR( 4,  3,  B3, 0,  0,  0,  0),     /*  B3 (S=0)    */
	/*  138 */  CR( 6,  2,  F3, 0,  0,  0,  0),     /*  F3          */
	/*  139 */  CR( 7,  1,  B4, 0,  0,  0,  0),     /*  B4          */
	/*  140 */  CR( 6,  2,  F4, 0,  0,  0,  0),     /*  F4          */
	/*  141 */  CR( 8,  0,  B5, 0,  0,  0,  0),     /*  B5          */
	/*  142 */  CR( 8,  0,  F5, 0,  0,  0,  0),     /*  F5          */
				/* LOAD_2 only              */
	/*  143 */  CR( 5,  0,  IA, 0,  0,  0,  0),     /*  Ampl. Intr. */
	/*  144 */  CR( 5,  0,  IP, 0,  0,  0,  0),     /*  Per. Intr.  */

				/* LOAD_2, LOAD_C  Mode 11  */
	/*  145 */  CR( 6,  2,  AM, 0,  0,  0,  1),     /*  Amplitude   */
	/*  146 */  CR( 8,  0,  PR, 0,  0,  0,  0),     /*  Period      */
	/*  147 */  CR( 6,  1,  B0, 0,  0,  0,  0),     /*  B0 (S=0)    */
	/*  148 */  CR( 6,  2,  F0, 0,  0,  0,  0),     /*  F0          */
	/*  149 */  CR( 6,  1,  B1, 0,  0,  0,  0),     /*  B1 (S=0)    */
	/*  150 */  CR( 6,  2,  F1, 0,  0,  0,  0),     /*  F1          */
	/*  151 */  CR( 6,  1,  B2, 0,  0,  0,  0),     /*  B2 (S=0)    */
	/*  152 */  CR( 6,  2,  F2, 0,  0,  0,  0),     /*  F2          */
	/*  153 */  CR( 6,  1,  B3, 0,  0,  0,  0),     /*  B3 (S=0)    */
	/*  154 */  CR( 7,  1,  F3, 0,  0,  0,  0),     /*  F3          */
	/*  155 */  CR( 8,  0,  B4, 0,  0,  0,  0),     /*  B4          */
	/*  156 */  CR( 8,  0,  F4, 0,  0,  0,  0),     /*  F4          */
	/*  157 */  CR( 8,  0,  B5, 0,  0,  0,  0),     /*  B5          */
	/*  158 */  CR( 8,  0,  F5, 0,  0,  0,  0),     /*  F5          */
				/* LOAD_2 only              */
	/*  159 */  CR( 5,  0,  IA, 0,  0,  0,  0),     /*  Ampl. Intr. */
	/*  160 */  CR( 5,  0,  IP, 0,  0,  0,  0),     /*  Per. Intr.  */

	/* -------------------------------------------------------------------- */
	/*  Opcode 0011: SETMSB_3                                               */
	/*  Opcode 0101: SETMSB_5                                               */
	/* -------------------------------------------------------------------- */
				/* Mode 00 only             */
	/*  161 */  CR( 0,  0,  0,  0,  0,  1,  0),     /*  Clear 5     */
				/* Mode 00 and 01           */
	/*  162 */  CR( 6,  2,  AM, 0,  0,  0,  0),     /*  Amplitude   */
	/*  163 */  CR( 8,  0,  PR, 0,  0,  0,  0),     /*  Period      */
	/*  164 */  CR( 5,  3,  F0, 0,  1,  0,  0),     /*  F0 (5 MSBs) */
	/*  165 */  CR( 5,  3,  F1, 0,  1,  0,  0),     /*  F1 (5 MSBs) */
	/*  166 */  CR( 5,  3,  F2, 0,  1,  0,  0),     /*  F2 (5 MSBs) */
				/* SETMSB_3 only            */
	/*  167 */  CR( 5,  0,  IA, 0,  0,  0,  0),     /*  Ampl. Intr. */
	/*  168 */  CR( 5,  0,  IP, 0,  0,  0,  0),     /*  Per. Intr.  */

				/* Mode 10 only             */
	/*  169 */  CR( 0,  0,  0,  0,  0,  1,  0),     /*  Clear 5     */
				/* Mode 10 and 11           */
	/*  170 */  CR( 6,  2,  AM, 0,  0,  0,  0),     /*  Amplitude   */
	/*  171 */  CR( 8,  0,  PR, 0,  0,  0,  0),     /*  Period      */
	/*  172 */  CR( 6,  2,  F0, 0,  1,  0,  0),     /*  F0 (6 MSBs) */
	/*  173 */  CR( 6,  2,  F1, 0,  1,  0,  0),     /*  F1 (6 MSBs) */
	/*  174 */  CR( 6,  2,  F2, 0,  1,  0,  0),     /*  F2 (6 MSBs) */
				/* SETMSB_3 only            */
	/*  175 */  CR( 5,  0,  IA, 0,  0,  0,  0),     /*  Ampl. Intr. */
	/*  176 */  CR( 5,  0,  IP, 0,  0,  0,  0),     /*  Per. Intr.  */
};

const INT16 sp0256_df_idx[16 * 8]    =
{
	/*  OPCODE 0000 */      -1, -1,     -1, -1,     -1, -1,     -1, -1,
	/*  OPCODE 1000 */      -1, -1,     -1, -1,     -1, -1,     -1, -1,
	/*  OPCODE 0100 */      17, 22,     17, 24,     25, 30,     25, 32,
	/*  OPCODE 1100 */      83, 94,     129,142,    97, 108,    145,158,
	/*  OPCODE 0010 */      83, 96,     129,144,    97, 110,    145,160,
	/*  OPCODE 1010 */      73, 77,     74, 77,     78, 82,     79, 82,
	/*  OPCODE 0110 */      33, 36,     34, 37,     38, 41,     39, 42,
	/*  OPCODE 1110 */      127,128,    127,128,    127,128,    127,128,
	/*  OPCODE 0001 */      1,  14,     1,  16,     1,  14,     1,  16,
	/*  OPCODE 1001 */      45, 56,     45, 58,     59, 70,     59, 72,
	/*  OPCODE 0101 */      161,166,    162,166,    169,174,    170,174,
	/*  OPCODE 1101 */      111,116,    111,118,    119,124,    119,126,
	/*  OPCODE 0011 */      161,168,    162,168,    169,176,    170,176,
	/*  OPCODE 1011 */      -1, -1,     -1, -1,     -1, -1,     -1, -1,
	/*  OPCODE 0111 */      -1, -1,     -1, -1,     -1, -1,     -1, -1,
	/*  OPCODE 1111 */      0,  0,      0,  0,      0,  0,      0,  0
};

/* ======================================================================== */
/*  BITREV32       -- Bit-reverse a 32-bit number.                            */
/* ======================================================================== */
static inline UINT32 bitrev32(UINT32 val)
{
	val = ((val & 0xFFFF0000) >> 16) | ((val & 0x0000FFFF) << 16);
	val = ((val & 0xFF00FF00) >>  8) | ((val & 0x00FF00FF) <<  8);
	val = ((val & 0xF0F0F0F0) >>  4) | ((val & 0x0F0F0F0F) <<  4);
	val = ((val & 0xCCCCCCCC) >>  2) | ((val & 0x33333333) <<  2);
	val = ((val & 0xAAAAAAAA) >>  1) | ((val & 0x55555555) <<  1);

	return val;
}

/* ======================================================================== */
/*  BITREV8       -- Bit-reverse a 8-bit number.                            */
/* ======================================================================== */
static inline UINT8 bitrev8(UINT8 val)
{
	val = ((val & 0xF0) >>  4) | ((val & 0x0F) <<  4);
	val = ((val & 0xCC) >>  2) | ((val & 0x33) <<  2);
	val = ((val & 0xAA) >>  1) | ((val & 0x55) <<  1);

	return val;
}

/* ======================================================================== */
/*  BITREVBUFF       -- Bit-reverse a buffer.                               */
/* ======================================================================== */
void bitrevbuff(UINT8 *buffer, unsigned int start, unsigned int length)
{
	for (unsigned int i = start; i < length; i++ )
		buffer[i] = bitrev8(buffer[i]);
}

/* ======================================================================== */
/*  SP0256_GETB  -- Get up to 8 bits at the current PC.                     */
/* ======================================================================== */
UINT32 getb( int len )
{
  UINT32 data, minus;
	u16 d0, d1;

	//	fprintf(stderr,"m_pc %d\n",m_pc>>3);

	/* -------------------------------------------------------------------- */
	/*  Fetch data from the FIFO or from the MASK                           */
	/* -------------------------------------------------------------------- */
	if (m_fifo_sel)
	{
		d0 = m_fifo[(m_fifo_tail    ) & 63];
		d1 = m_fifo[(m_fifo_tail + 1) & 63];

		data = ((d1 << 10) | d0) >> m_fifo_bitp;


		/* ---------------------------------------------------------------- */
		/*  Note the PC doesn't advance when we execute from FIFO.          */
		/*  Just the FIFO's bit-pointer advances.   (That's not REALLY      */
		/*  what happens, but that's roughly how it behaves.)               */
		/* ---------------------------------------------------------------- */
		m_fifo_bitp += len;
		if (m_fifo_bitp >= 10)
		{
			m_fifo_tail++;
			m_fifo_bitp -= 10;
		}
	} else
	{
		/* ---------------------------------------------------------------- */
		/*  Figure out which ROMs are being fetched into, and grab two      */
		/*  adjacent bytes.  The byte we're interested in is extracted      */
		/*  from the appropriate bit-boundary between them.                 */
		/* ---------------------------------------------------------------- */
	  int32_t idx0 = (m_pc    ) >> 3, d0; //???
	  int32_t idx1 = (m_pc + 8) >> 3, d1;
		
	  //	int firstadd=(idx0 & 0xffff);
	  //	int secondadd=(idx1 & 0xffff);

	  /*	  if (firstadd<0 || secondadd>0x800) data=0; // better check on this - NOT in case of m_rom19
				else
				{*/

	  data=0; 

	  /*	  if (idx0>=0x1000 && idx0<0x1800) {// this works but not for further roms
	      minus=0x1000;
	      data=1;
	    }
	    else data=0;*/

	  if (idx0<0x1800 && idx0>=0x1000) {
	  data=1; //
	  minus=0x1000; // default
	}

	  //else if (idx0>=0x1800 && idx0<0x4000){
	    //  fprintf(stderr, "fifo?????? 0x%X\n", idx0);
	  //	    data=0;
	  //	  }
	  else if (idx0>=0x4000 && idx0<0x8000) {
		    m_romm=m_rom003; // 003 has phonemes as AL2 and some phrases but not so many WHY?
		    minus=0x4000;
		    data=1;
	  }
	  else if (idx0>=0x8000 && idx0<0xC000) {
	    m_romm=m_rom004; 
	    minus=0x8000;
	    data=1;
	    }



	  if (data!=0){
	  int32_t firstadd=(idx0 & 0xffff)-minus;
	  int32_t secondadd=(idx1 & 0xffff)-minus;
	  

		d0 = m_romm[firstadd];
		d1 = m_romm[secondadd]; // was 0xffff

		data = ((d1 << 8) | d0) >> (m_pc & 7);
		//				}

		m_pc += len;
		  

	/* -------------------------------------------------------------------- */
	/*  Mask data to the requested length.                                  */
	/* -------------------------------------------------------------------- */
	data &= ((1 << len) - 1);
	  }
	}

	return data;
}

/* ======================================================================== */
/*  SP0256_MICRO -- Emulate the microsequencer in the SP0256.  Executes     */
/*                  instructions either until the repeat count != 0 or      */
/*                  the sequencer gets halted by a RTS to 0.                */
/* ======================================================================== */

//const float repeatre_map[32]={0.1f, 0.119587327404f, 0.143011288757f, 0.171023378111f, 0.20452228712f, 0.244582737113f, 0.292489958607f, 0.349780924424f, 0.418293659289f, 0.500226207846f, 0.598207152939f, 0.715379946542f, 0.855503758857f, 1.02307408106f, 1.22346695091f, 1.46311142827f, 1.74969585401f, 2.09241450952f, 2.50226259016f, 2.99238895621f, 3.57851797828f, 4.27945401091f, 5.11768467915f, 6.12010233278f, 7.31886681418f, 8.75243721937f, 10.4668057534f, 12.5169732651f, 14.9687137997f, 17.9006847798f, 21.4069505153f, 25.6f}; // how did we generate this?

void micro()
{
	UINT8  immed4;
	UINT8  opcode;
	UINT16 cr;
	u8 ctrl_xfer;
	int16_t repeat;
	u8 i, idx0, idx1;
	int16_t val;
	/* -------------------------------------------------------------------- */
	/*  Only execute instructions while the filter is not busy.             */
	/* -------------------------------------------------------------------- */
	while (m_filt.rpt <= 0)
	{
		/* ---------------------------------------------------------------- */
		/*  If the CPU is halted, see if we have a new command pending      */
		/*  in the Address LoaD buffer.                                     */
		/* ---------------------------------------------------------------- */
		if (m_halted && !m_lrq)
		{
		 	  m_pc       = m_ald | (0x1000 << 3); // OR with 0x8000 this adds 0x1000 which we subtract later when shifts back
			  // m_pc = m_ald;
			m_fifo_sel = 0;
			m_halted   = 0;
			m_lrq      = 0x8000;
			m_ald      = 0;
			for (i = 0; i < 16; i++)
				m_filt.r[i] = 0;
			//			m_drq_cb(1);
		}

		/* ---------------------------------------------------------------- */
		/*  If we're still halted, do nothing.                              */
		/* ---------------------------------------------------------------- */
		if (m_halted)
		{
			m_filt.rpt = 1;
			m_lrq      = 0x8000;
			m_ald      = 0;
			for (i = 0; i < 16; i++)
				m_filt.r[i] = 0;

			//			SET_SBY(1)
			//			really=1;
			return;
		}

		/* ---------------------------------------------------------------- */
		/*  Fetch the first 8 bits of the opcode, which are always in the   */
		/*  same approximate format -- immed4 followed by opcode.           */
		/* ---------------------------------------------------------------- */
		immed4 = getb(4);
		opcode = getb(4);
		repeat = 0;
		ctrl_xfer = 0;

		//		printf("$%.4X.%.1X: OPCODE %d%d%d%d.%d%d\n",
		//		(m_pc >> 3) - 1, m_pc & 7,
		//		!!(opcode & 1), !!(opcode & 2),
		//		!!(opcode & 4), !!(opcode & 8),
		//		!!(m_mode&4), !!(m_mode&2));

		/* ---------------------------------------------------------------- */
		/*  Handle the special cases for specific opcodes.                  */
		/* ---------------------------------------------------------------- */
		switch (opcode)
		{
			/* ------------------------------------------------------------ */
			/*  OPCODE 0000:  RTS / SETPAGE                                 */
			/* ------------------------------------------------------------ */
			case 0x0:
			{
				/* -------------------------------------------------------- */
				/*  If immed4 != 0, then this is a SETPAGE instruction.     */
				/* -------------------------------------------------------- */
				if (immed4)     /* SETPAGE */
				{
					m_page = bitrev32(immed4) >> 13;
				} else
				/* -------------------------------------------------------- */
				/*  Otherwise, this is an RTS / HLT.                        */
				/* -------------------------------------------------------- */
				{
					UINT32 btrg;

					/* ---------------------------------------------------- */
					/*  Figure out our branch target.                       */
					/* ---------------------------------------------------- */
					btrg = m_stack;

					m_stack = 0;

					/* ---------------------------------------------------- */
					/*  If the branch target is zero, this is a HLT.        */
					/*  Otherwise, it's an RTS, so set the PC.              */
					/* ---------------------------------------------------- */
					if (!btrg)
					{
						m_halted   = 1;
						m_pc       = 0;
						ctrl_xfer  = 1;
					} else
					{
						m_pc      = btrg;
						ctrl_xfer = 1;
					}
				}

				break;
			}

			/* ------------------------------------------------------------ */
			/*  OPCODE 0111:  JMP          Jump to 12-bit/16-bit Abs Addr   */
			/*  OPCODE 1011:  JSR          Jump to Subroutine               */
			/* ------------------------------------------------------------ */
			case 0xE:
			case 0xD:
			{
				uint32_t btrg;

				/* -------------------------------------------------------- */
				/*  Figure out our branch target.                           */
				/* -------------------------------------------------------- */
				btrg = m_page                     |
						(bitrev32(immed4)  >> 17) |
						(bitrev32(getb(8)) >> 21);
				ctrl_xfer = 1;

				/* -------------------------------------------------------- */
				/*  If this is a JSR, push our return address on the        */
				/*  stack.  Make sure it's byte aligned.                    */
				/* -------------------------------------------------------- */
				if (opcode == 0xD)
					m_stack = (m_pc + 7) & ~7;

				/* -------------------------------------------------------- */
				/*  Jump to the new location!                               */
				/* -------------------------------------------------------- */
				m_pc = btrg;
				break;
			}

			/* ------------------------------------------------------------ */
			/*  OPCODE 1000:  SETMODE      Set the Mode and Repeat MSBs     */
			/* ------------------------------------------------------------ */
			case 0x1:
			{
				m_mode = ((immed4 & 8) >> 2) | (immed4 & 4) | ((immed4 & 3) << 4);
				break;
			}

			/* ------------------------------------------------------------ */
			/*  OPCODE 0001:  LOADALL      Load All Parameters              */
			/*  OPCODE 0010:  LOAD_2       Load Per, Ampl, Coefs, Interp.   */
			/*  OPCODE 0011:  SETMSB_3     Load Pitch, Ampl, MSBs, & Intrp  */
			/*  OPCODE 0100:  LOAD_4       Load Pitch, Ampl, Coeffs         */
			/*  OPCODE 0101:  SETMSB_5     Load Pitch, Ampl, and Coeff MSBs */
			/*  OPCODE 0110:  SETMSB_6     Load Ampl, and Coeff MSBs.       */
			/*  OPCODE 1001:  DELTA_9      Delta update Ampl, Pitch, Coeffs */
			/*  OPCODE 1010:  SETMSB_A     Load Ampl and MSBs of 3 Coeffs   */
			/*  OPCODE 1100:  LOAD_C       Load Pitch, Ampl, Coeffs         */
			/*  OPCODE 1101:  DELTA_D      Delta update Ampl, Pitch, Coeffs */
			/*  OPCODE 1110:  LOAD_E       Load Pitch, Amplitude            */
			/*  OPCODE 1111:  PAUSE        Silent pause                     */
			/* ------------------------------------------------------------ */
			default:
			{			  
			  repeat = immed4 | (m_mode & 0x30);
			  if (modus!=4) {// TTS mode
			    val = (_sely*131.0f); // only if we are not in TTS mode
			    MAXED(val,127);
			    if (modus==1) repeat=64.0f*logpitch[val]; // make much longer for singing
			    else if (modus==16) repeat=((float)(repeat)*8.0f*logpitch[val]);
			    else repeat=((float)(repeat)*logpitch[val]);
			  }
			  if (repeat<1) repeat=1;
			  break;
			}
		}
		if (opcode != 1) m_mode &= 0xF;

		/* ---------------------------------------------------------------- */
		/*  If this was a control transfer, handle setting "fifo_sel"       */
		/*  and all that ugliness.                                          */
		/* ---------------------------------------------------------------- */
		if (ctrl_xfer)
		{

			/* ------------------------------------------------------------ */
			/*  Set our "FIFO Selected" flag based on whether we're going   */
			/*  to the FIFO's address.                                      */
			/* ------------------------------------------------------------ */
			m_fifo_sel = m_pc == FIFO_ADDR;


			/* ------------------------------------------------------------ */
			/*  Control transfers to the FIFO cause it to discard the       */
			/*  partial decle that's at the front of the FIFO.              */
			/* ------------------------------------------------------------ */
			if (m_fifo_sel && m_fifo_bitp)
			{

				/* Discard partially-read decle. */
				if (m_fifo_tail < m_fifo_head) m_fifo_tail++;
				m_fifo_bitp = 0;
			}


			continue;
		}

		/* ---------------------------------------------------------------- */
		/*  Otherwise, if we have a repeat count, then go grab the data     */
		/*  block and feed it to the filter.                                */
		/* ---------------------------------------------------------------- */
		if (!repeat) continue;

		m_filt.rpt = repeat + 1;

		i = (opcode << 3) | (m_mode & 6);
		idx0 = sp0256_df_idx[i++];
		idx1 = sp0256_df_idx[i  ];

		//		assert(idx0 >= 0 && idx1 >= 0 && idx1 >= idx0);

		/* ---------------------------------------------------------------- */
		/*  Step through control words in the description for data block.   */
		/* ---------------------------------------------------------------- */
		for (i = idx0; i <= idx1; i++)
		{
			int16_t len, shf, delta, field, prm, clra, clr5;
			INT8 value;

			/* ------------------------------------------------------------ */
			/*  Get the control word and pull out some important fields.    */
			/* ------------------------------------------------------------ */
			cr = sp0256_datafmt[i];

			len   = CR_LEN(cr);
			shf   = CR_SHF(cr);
			prm   = CR_PRM(cr);
			clra  = cr & CR_CLRA;
			clr5  = cr & CR_CLR5;
			delta = cr & CR_DELTA;
			field = cr & CR_FIELD;
			value = 0;

			/* ------------------------------------------------------------ */
			/*  Clear any registers that were requested to be cleared.      */
			/* ------------------------------------------------------------ */
			if (clra)
			{
				for (u8 j = 0; j < 16; j++)
					m_filt.r[j] = 0;

				m_silent = 1;
			}

			if (clr5)
				m_filt.r[B5] = m_filt.r[F5] = 0;

			/* ------------------------------------------------------------ */
			/*  If this entry has a bitfield with it, grab the bitfield.    */
			/* ------------------------------------------------------------ */
			if (len)
			{
				value = getb(len);
			}
			else
			{
				continue;
			}

			/* ------------------------------------------------------------ */
			/*  Sign extend if this is a delta update.                      */
			/* ------------------------------------------------------------ */
			if (delta)  /* Sign extend */
			{
				if (value & (1 << (len - 1))) value |= -1 << len;
			}

			/* ------------------------------------------------------------ */
			/*  Shift the value to the appropriate precision.               */
			/* ------------------------------------------------------------ */
			if (shf)
				value <<= shf;


			m_silent = 0;

			/* ------------------------------------------------------------ */
			/*  If this is a field-replace, insert the field.               */
			/* ------------------------------------------------------------ */
			if (field)
			{

				m_filt.r[prm] &= ~(~0 << shf); /* Clear the old bits.     */
				m_filt.r[prm] |= value;        /* Merge in the new bits.  */


				continue;
			}

			/* ------------------------------------------------------------ */
			/*  If this is a delta update, add to the appropriate field.    */
			/* ------------------------------------------------------------ */
			if (delta)
			{

				m_filt.r[prm] += value;

				continue;
			}

			/* ------------------------------------------------------------ */
			/*  Otherwise, just write the new value.                        */
			/* ------------------------------------------------------------ */
			m_filt.r[prm] = value;
		}

		/* ---------------------------------------------------------------- */
		/*  Special case:  Set PAUSE's equivalent period.                   */
		/* ---------------------------------------------------------------- */
		if (opcode == 0xF)
		{
			m_silent = 1;
			m_filt.r[1] = PER_PAUSE;
		}

		/* ---------------------------------------------------------------- */
		/*  Now that we've updated the registers, go decode them.           */
		/* ---------------------------------------------------------------- */
		lpc12_regdec(&m_filt);

		/* ---------------------------------------------------------------- */
		/*  Break out since we now have a repeat count.                     */
		/* ---------------------------------------------------------------- */
		break;
	}
}


int16_t sp0256_get_sample1219(void){
  static int16_t output; 
  modus=3;
  micro();
  lpc12_update(&m_filt, &output);
   
   if (m_halted==1 && m_filt.rpt <= 0)     {
     sp0256_newsay1219();
   }
   return output;
}

const unsigned char remap19[49]   ={64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 115, 116, 117, 118, 119, 120, 121, 122, 123, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28}; // there are 49

void sp0256_newsay1219(void){
  u8 dada, indexy;
  m_lrq=0; m_halted=1; m_filt.rpt=0;

  u8 selector=_selz*89.0f; // total is 36+49=85
  MAXED(selector, 85); //0-85
  selector=85-selector; // inverted
  
  //  selector=34;

  if (selector<37) { // so top is 36
    m_page=0x1000<<3;
   m_romm=m_rom12;
   dada=6+selector; // they are 6->42
    }
    else {    
      m_romm=m_rom19;
      indexy=selector-37;
      dada=remap19[indexy]; // checked
      
      if (indexy>19) m_page=0x8000<<3;
      else m_page=0x1000<<3;
    }
      m_ald = ((dada) << 4); // or do as index <<3 and store this index 		
      m_lrq = 0; //from 8 bit write
    }


int16_t sp0256_get_sample(void){
  static int16_t output; 
  modus=16;
  micro();
  lpc12_update(&m_filt, &output);
  m_lrq=0x8000;
  
  if (m_halted==1)     {
  //  if (really==1){
  //    really=0;
	sp0256_newsay();
   }
   return output;
 }

int16_t sp0256_get_sample_sing(void){
  static int16_t output; 
  modus=1; // singing
  micro();
  lpc12_update(&m_filt, &output);
  
  if (m_halted==1 && m_filt.rpt <= 0)     {
    sp0256_newsay();
   }
  
   return output;
 }


 int16_t sp0256_get_sampleTTS(void){
  static int16_t output; 
  modus=4;

   micro();
   lpc12_update(&m_filt, &output);

   
   if (m_halted==1 && m_filt.rpt <= 0)     {
          sp0256_newsayTTS();
   }

   return output;
 }


 int16_t sp0256_get_samplevocabbankone(void){
  static int16_t output; 
  modus=3;
   micro();
   lpc12_update(&m_filt, &output);

   if (m_halted==1 && m_filt.rpt <= 0)     {
          sp0256_newsayvocabbankone(0);
   }

   
return output;
 }

 int16_t sp0256_get_samplevocabbanktwo(void){
  static int16_t output;   
  modus=3; // was zero

  micro();
  lpc12_update(&m_filt, &output);
   
   if (m_halted==1 && m_filt.rpt <= 0)     {
          sp0256_newsayvocabbanktwo(0);
   }

   return output;
 }

extern char TTSinarray[17];
static u8 TTSoutarray[256];
static u8 TTSindex=0;
static u8 TTSlength=0;

inline void sp0256_newsay(void){
   u8 dada;
   dada=_selz*67.0f; 
   MAXED(dada,63);
   dada=63-dada;
   //   dada=31;
   m_lrq=0; m_halted=1; m_filt.rpt=0;
   m_page     = 0x1000 << 3; //32768 =0x8000
   m_romm=m_romAL2;
   m_ald = ((dada&63) << 4);  		
   m_lrq = 0; //from 8 bit write
 }

void sp0256_newsayTTS(void){// called at end of phoneme
   u8 dada;
   //   m_halted=1;
   m_lrq=0; m_halted=1; m_filt.rpt=0;
   m_page     = 0x1000 << 3; //32768 =0x8000
   m_romm=m_romAL2;
   
   dada=TTSoutarray[TTSindex]; 
   TTSindex++;
   if (TTSindex>=TTSlength) {
     TTSindex=0;
     TTSlength= text2speechfor256(16,TTSinarray,TTSoutarray);
   }

   m_ald = ((dada&63) << 4); 
   m_lrq = 0; //from 8 bit write
 }

// new code to retrigger TTS

void sp0256_retriggerTTS(void){// called on an outside trigger
   u8 dada;
   //   m_halted=1;
   m_lrq=0; m_halted=1; m_filt.rpt=0;
   m_page     = 0x1000 << 3; //32768 =0x8000
   m_romm=m_romAL2;
   
   //     TTSindex=0;
     TTSlength= text2speechfor256(16,TTSinarray,TTSoutarray);
     dada=TTSoutarray[0]; 
     TTSindex=1;
     m_ald = ((dada&63) << 4); 
     m_lrq = 0; //from 8 bit write
 }

void sp0256_newsayvocabbankonea(){ // calls a reset back to start of phoneme vocab
  sp0256_newsayvocabbankone(1);
}

void sp0256_newsayvocabbanktwoa(){
  sp0256_newsayvocabbanktwo(1);
}


void sp0256_newsayvocabbankone(u8 reset){// called at end of phoneme
   u8 dada;
   m_lrq=0; m_halted=1; m_filt.rpt=0;
   static u8 vocabindex=0, whichone=0;
   //   m_halted=1;

   m_page     = 0x1000 << 3; //32768 =0x8000
   m_romm=m_romAL2;   
   if (reset==1 || *(vocab_sp0256_bankone[whichone]+vocabindex)==255 ){ // safer so reset==1  never evaluates past end of array
     vocabindex=0;
     whichone=_selz*92.0f; // split vocab into banks - on this one we need 0-150 values 
     MAXED(whichone,87);
     whichone=87-whichone;
   }   

   dada=*(vocab_sp0256_bankone[whichone]+vocabindex);  
   vocabindex++;

   m_ald = ((dada&63) << 4); 		
   m_lrq = 0; //from 8 bit write
 }


void sp0256_newsayvocabbanktwo(u8 reset){// called at end of phoneme
   u8 dada;
   m_lrq=0;m_halted=1;      m_filt.rpt=0;
   static u8 vocabindex=0, whichone=0;
   //   m_halted=1;

   m_page     = 0x1000 << 3; //32768 =0x8000
   m_romm=m_romAL2;   
   if (reset==1 || *(vocab_sp0256_banktwo[whichone]+vocabindex)==255){
     vocabindex=0;
     whichone=_selz*94.0f; // this is a different bank so we have 94
     MAXED(whichone,89);
     whichone=89-whichone;
   }
   dada=*(vocab_sp0256_banktwo[whichone]+vocabindex);  // in this case end switch
   vocabindex++;
   m_ald = ((dada&63) << 4);  		
   //   m_lrq = 0; //from 8 bit write
 }


 void sp0256_init(void){
   sp0256_iinit();
   reset();
 }
