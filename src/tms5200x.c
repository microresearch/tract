// license:BSD-3-Clause
// copyright-holders:Frank Palazzolo, Aaron Giles, Jonathan Gevaryahu, Raphael Nabet, Couriersud, Michael Zapf

// modified for WORM by Martin Howse


#include "audio.h"
#include "stdio.h"
#include "math.h"
#include "time.h"

//extern uint32_t triggermode;
extern float _selx, _sely, _selz;
extern float exy[64][240];
static u8 modus=0;
extern u8 _intmode;

const u8 koffset[10]={0, 32, 64, 80, 96, 112, 128, 144, 152, 160};

typedef u8 UINT8;
typedef uint16_t UINT16;
typedef int16_t INT16;
typedef uint32_t UINT32;
typedef int32_t INT32;

#include "tms5110r.inc"
#include "english2phoneme/TTS.h"
#include "resources.h"

//** 5100-speak and spell
//** 5200- early echo II, disks of TRON????, TI99/4 
//** 5220- later echo II, BBC MICRO

// +2 for spanish

#define TMS_VOCAB 72 // TOTAL
#define TMS_VOCAB_TOP 71 // 0-
#define TMS_VOCAB_F 75.0f 
#define ALLPHON_BANK 0 // CHECKed - but moved alphon!
//#define EXTRAX 55 // this is the extra vocab - replaced +3

// group first 5100:

// refreshing

// then 5200
#include "LPC/roms/vocab_testroms.h" // 5200=wordlist_alphons[127] // and // wordlist_arcade_one[74] // and // wordlist_arcade_two[67]
#include "LPC/roms/vocaid.h"
// list of which_coeff:

// 5100: &T0280B_0281A_coeff
// 5200: &T0285_2501E_coeff
// 5220: &tms5220_coeff

typedef struct TMS_vocab__ {
  const uint8_t **wordlist;
  const struct tms5100_coeffs *m_coeff;
  uint16_t extent;
  float extentplus;
} TMS_vocab;

// alt coeffs for spk and spell

// T0280D_0281D_coeff - later spk and spell
// T0280F_2801A_coeff - speak and math.spk and read??? - english spk and spell CD2801

const TMS_vocab vocab_alphons={wordlist_alphons, &T0285_2501E_coeff, 126, 130.0f};//ti99 alphons
const TMS_vocab vocab_vocaid={wordlist_spellvocaid, &T0280F_2801A_coeff, 123, 127.0f}; ///*

static u8 last5100=1;
static u8 last5200=58;

const TMS_vocab *allTMSvocabs[TMS_VOCAB]={&vocab_alphons, &vocab_vocaid};

//&vocab_tntellET

//LIST: &0vocab_2303, &1vocab_2304, &2vocab_2321, &3vocab_2322, &4vocab_2350_1, &5vocab_2350_2, &6vocab_2352_1, &7vocab_2352_2, &8vocab_spkspellone_1, &9vocab_spkspellone_2, &10vocab_spkspelltwo_1, &11vocab_spkspelltwo_2, &12vocab_D000_1, &13vocab_D000_2, &14vocab_D001_1, &15vocab_D001_2, &16vocab_D002_1, &17vocab_D002_2, &18vocab_D003, &19vocab_D004_1, &20vocab_D004_2, &21vocab_D005_1, &22vocab_D005_2, &23vocab_D006_1, &24vocab_D006_2, &25vocab_D007_1, &26vocab_D007_2, &27vocab_D008, &28vocab_D009_1, &29vocab_D009_2, &30vocab_D010_1, &31vocab_D010_2, &32vocab_D011, &33vocab_D012_1, &34vocab_D012_2, &35vocab_D013, &36vocab_D014, &37vocab_D015_1, &38vocab_D015_2, &39vocab_D016, &40vocab_D017, &41vocab_D018, &42vocab_D019, &43vocab_D020_21, &44vocab_D022, &45vocab_D023, &46vocab_D024_25, &47vocab_D026, &48vocab_D027_34, &49vocab_echofemale_1, &50vocab_echofemale_2, &51vocab_mpf, &52vocab_alphons, &53vocab_arcade_one, &54vocab_arcade_two, &55vocab_ti99_one_1, &56vocab_ti99_one_2, &57vocab_ti99_two_1, &58vocab_ti99_two_2, &59vocab_acorn_1, &60vocab_acorn_2, &61vocab_male_1, &62vocab_male_2};//, &63vocab_large_male_one_1, &64vocab_large_male_one_2, &65vocab_large_male_two_1, &66vocab_large_male_two_2}; // 67

//};//


//const TMS_vocab *allTMSvocabs[46]={&vocab_2303, &vocab_alphons};

static const uint8_t* ptrAddr; static uint8_t ptrBit;
uint8_t byte_rev[256];

extern char TTSinarray[17];
static u8 TTSoutarray[256];
static u8 TTSindex=0;
static u8 TTSlength=0;


#define PERFECT_INTERPOLATION_HACK 0

static INT16 clip_analog(INT16 cliptemp);

/* *****optional defines***** */

/* Hacky improvements which don't match patent: */
/* Interpolation shift logic:
 * One of the following two lines should be used, and the other commented
 * The second line is more accurate mathematically but not accurate to the patent
 */
#define INTERP_SHIFT >> m_coeff->interp_coeff[m_IP]
//#define INTERP_SHIFT / (1<<m_coeff->interp_coeff[m_IP])

/* Other hacks */
/* HACK?: if defined, outputs the low 4 bits of the lattice filter to the i/o
 * or clip logic, even though the real hardware doesn't do this, partially verified by decap */
#undef ALLOW_4_LSB
//#define ALLOW_4_LSB 1

/* *****configuration of chip connection stuff***** */
/* must be defined; if 0, output the waveform as if it was tapped on the speaker pin as usual, if 1, output the waveform as if it was tapped on the i/o pin (volume is much lower in the latter case) */
#define FORCE_DIGITAL 0

/* must be defined; if 1, normal speech (one A cycle, one B cycle per interpolation step); if 0; speak as if SPKSLOW was used (two A cycles, one B cycle per interpolation step) */
#define FORCE_SUBC_RELOAD 1 // appears to be unused


/* *****debugging defines***** */
#undef VERBOSE
// above is general, somewhat obsolete, catch all for debugs which don't fit elsewhere
#define DEBUG_DUMP_INPUT_DATA 1
// above dumps the data input to the tms52xx to stdout, useful for making logged data dumps for real hardware tests
#undef DEBUG_FIFO
// above debugs fifo stuff: writes, reads and flag updates
#undef DEBUG_PARSE_FRAME_DUMP
// above dumps each frame to stderr: be sure to select one of the options below if you define it!
#undef DEBUG_PARSE_FRAME_DUMP_BIN
// dumps each speech frame as binary
#undef DEBUG_PARSE_FRAME_DUMP_HEX
// dumps each speech frame as hex
#undef DEBUG_FRAME_ERRORS
// above dumps info if a frame ran out of data
#undef DEBUG_COMMAND_DUMP
// above dumps all non-speech-data command writes
#undef DEBUG_PIN_READS
// above spams the errorlog with i/o ready messages whenever the ready or irq pin is read
#undef DEBUG_GENERATION
// above dumps debug information related to the sample generation loop, i.e. whether interpolation is inhibited or not, and what the current and target values for each frame are.
#undef DEBUG_GENERATION_VERBOSE
// above dumps MUCH MORE debug information related to the sample generation loop, namely the excitation, energy, pitch, k*, and output values for EVERY SINGLE SAMPLE during a frame.
#undef DEBUG_LATTICE
// above dumps the lattice filter state data each sample.
#undef DEBUG_CLIP
// above dumps info to stderr whenever the analog clip hardware is (or would be) clipping the signal.
#undef DEBUG_IO_READY
// above debugs the io ready callback
#undef DEBUG_RS_WS
// above debugs the tms5220_data_r and data_w access methods which actually respect rs and ws

#define MAX_SAMPLE_CHUNK    512

/* Variants */

#define TMS5220_IS_5220C    (4)
#define TMS5220_IS_5200     (5)
#define TMS5220_IS_5220     (6)
#define TMS5220_IS_CD2501ECD (7)

#define TMS5220_IS_CD2501E  TMS5220_IS_5200

#define TMS5220_HAS_RATE_CONTROL ((m_variant == TMS5220_IS_5220C) || (m_variant == TMS5220_IS_CD2501ECD))
#define TMS5220_IS_52xx ((m_variant == TMS5220_IS_5220C) || (m_variant == TMS5220_IS_5200) || (m_variant == TMS5220_IS_5220) || (m_variant == TMS5220_IS_CD2501ECD))


static const UINT8 reload_table[4] = { 0, 2, 4, 6 }; //sample count reload for 5220c and cd2501ecd only; 5200 and 5220 always reload with 0; keep in mind this is loaded on IP=0 PC=12 subcycle=1 so it immediately will increment after one sample, effectively being 1,3,5,7 as in the comments above.

	int m_variant=0;

	/* coefficient tables */
const struct tms5100_coeffs *m_coeff; // this is fine with MAX_K (max number of coeffs)

	/* these contain global status bits */
	UINT8 m_speaking_now;     /* True only if actual speech is being generated right now. Is set when a speak vsm command happens OR when speak external happens and buffer low becomes nontrue; Is cleared when speech halts after the last stop frame or the last frame after talk status is otherwise cleared.*/
	UINT8 m_speak_external;   /* If 1, DDIS is 1, i.e. Speak External command in progress, writes go to FIFO. */
UINT8 m_talk_status;      /* If 1, TS status bit is 1, i.e. speak or speak external is in progress and we have not encountered a stop frame yet; */

	/* these contain global status bits */
	UINT8 m_previous_TALK_STATUS;      /* this is the OLD value of TALK_STATUS (i.e. previous value of m_SPEN|m_TALKD), needed for generating interrupts on a falling TALK_STATUS edge */
	UINT8 m_SPEN;             /* set on speak(or speak external and BL falling edge) command, cleared on stop command, reset command, or buffer out */
	UINT8 m_DDIS;             /* If 1, DDIS is 1, i.e. Speak External command in progress, writes go to FIFO. */
	UINT8 m_TALK;             /* set on SPEN & RESETL4(pc12->pc0 transition), cleared on stop command or reset command */
#define TALK_STATUS (m_SPEN|m_TALKD)
	UINT8 m_TALKD;            /* TALK(TCON) value, latched every RESETL4 */
	UINT8 m_buffer_low;       /* If 1, FIFO has less than 8 bytes in it */
	UINT8 m_buffer_empty;     /* If 1, FIFO is empty */
	UINT8 m_irq_pin;          /* state of the IRQ pin (output) */
	UINT8 m_ready_pin;        /* state of the READY pin (output) */

	/* these contain data describing the current and previous voice frames */
#define OLD_FRAME_SILENCE_FLAG m_OLDE // 1 if E=0, 0 otherwise.
#define OLD_FRAME_UNVOICED_FLAG m_OLDP // 1 if P=0 (unvoiced), 0 if voiced
	UINT8 m_OLDE;
	UINT8 m_OLDP;

#define NEW_FRAME_STOP_FLAG (m_new_frame_energy_idx == 0xF) // 1 if this is a stop (Energy = 0xF) frame
#define NEW_FRAME_SILENCE_FLAG (m_new_frame_energy_idx == 0) // ditto as above
#define NEW_FRAME_UNVOICED_FLAG (m_new_frame_pitch_idx == 0) // ditto as above
	UINT8 m_new_frame_energy_idx;
	UINT8 m_new_frame_pitch_idx;
	UINT8 m_new_frame_k_idx[10];

	INT16 m_target_energy;
	INT16 m_target_pitch;
	INT16 m_target_k[10];


	/* these are all used to contain the current state of the sound generation */
#ifndef PERFECT_INTERPOLATION_HACK
	INT16 m_current_energy;
	INT16 m_current_pitch;
	INT16 m_current_k[10];
#else
	UINT8 m_old_frame_energy_idx;
	UINT8 m_old_frame_pitch_idx;
	UINT8 m_old_frame_k_idx[10];
	UINT8 m_old_zpar;
	UINT8 m_old_uv_zpar;

	INT32 m_current_energy;
	INT32 m_current_pitch;
	INT32 m_current_k[10];
#endif

	UINT16 m_previous_energy; /* needed for lattice filter to match patent */

	UINT8 m_subcycle;         /* contains the current subcycle for a given PC: 0 is A' (only used on SPKSLOW mode on 51xx), 1 is A, 2 is B */
UINT8 m_subc_reload=1;      /* contains 1 for normal speech, 0 when SPKSLOW is active */ // wasn't set before
	UINT8 m_PC;               /* current parameter counter (what param is being interpolated), ranges from 0 to 12 */
	/* NOTE: the interpolation periods counts 1,2,3,4,5,6,7,0 for divide by 8,8,8,4,4,2,2,1 */
	UINT8 m_IP;               /* the current interpolation period */
	UINT8 m_inhibit;          /* If 1, interpolation is inhibited until the DIV1 period */
	UINT8 m_uv_zpar;          /* If 1, zero k5 thru k10 coefficients */
	UINT8 m_zpar;             /* If 1, zero ALL parameters. */
	UINT8 m_pitch_zero;       /* circuit 412; pitch is forced to zero under certain circumstances */
	UINT8 m_c_variant_rate=0;    /* only relevant for tms5220C's multi frame rate feature; is the actual 4 bit value written on a 0x2* or 0x0* command */
	UINT16 m_pitch_count;     /* pitch counter; provides chirp rom address */

	INT32 m_u[11];
	INT32 m_x[10];

	UINT16 m_RNG;             /* the random noise generator configuration is: 1 + x + x^3 + x^4 + x^13 TODO: no it isn't */
	INT16 m_excitation_data;

	/* The TMS52xx has two different ways of providing output data: the
	   analog speaker pin (which was usually used) and the Digital I/O pin.
	   The internal DAC used to feed the analog pin is only 8 bits, and has the
	   funny clipping/clamping logic, while the digital pin gives full 10 bit
	   resolution of the output data.
	   set/reset this other than the FORCE_DIGITAL define
	 */
	UINT8 m_digital_select=1;

	/* io_ready: page 3 of the datasheet specifies that READY will be asserted until
	 * data is available or processed by the system.
	 */
	UINT8 m_io_ready;

	/* flag for "true" timing involving rs/ws */
	UINT8 m_true_timing;


// helpers

int extract_bits(int count) // extract from rom image/array
{

  uint8_t value;
  UINT16 data;
  INT8 num_bits=count;

  	data = byte_rev[*ptrAddr]<<8;
	if (ptrBit+num_bits > 8)
	{
	    data |= byte_rev[*(ptrAddr+1)];
	}
	data <<= ptrBit;
	value = data >> (16-num_bits);
	ptrBit += num_bits;
	if (ptrBit >= 8)
	{
		ptrBit -= 8;
		ptrAddr++;
		if (ptrBit==0) {
		}
	}
	return value;
}


/**********************************************************************************************

     clip_analog -- clips the 14 bit return value from the lattice filter to its final 10 bit value (-512 to 511), and upshifts/range extends this to 16 bits

***********************************************************************************************/

static INT16 clip_analog(INT16 cliptemp)
{
	/* clipping, just like the patent shows:
	 * the top 10 bits of this result are visible on the digital output IO pin.
	 * next, if the top 3 bits of the 14 bit result are all the same, the lowest of those 3 bits plus the next 7 bits are the signed analog output, otherwise the low bits are all forced to match the inverse of the topmost bit, i.e.:
	 * 1x xxxx xxxx xxxx -> 0b10000000
	 * 11 1bcd efgh xxxx -> 0b1bcdefgh
	 * 00 0bcd efgh xxxx -> 0b0bcdefgh
	 * 0x xxxx xxxx xxxx -> 0b01111111
	 */
	if (cliptemp > 2047) cliptemp = 2047;
	else if (cliptemp < -2048) cliptemp = -2048;
	/* at this point the analog output is tapped */
#ifdef ALLOW_4_LSB
	// input:  ssss snnn nnnn nnnn
	// N taps:       ^^^ ^         = 0x0780
	// output: snnn nnnn nnnn NNNN
	return (cliptemp << 4)|((cliptemp&0x780)>>7); // upshift and range adjust
#else
	cliptemp &= ~0xF;
	// input:  ssss snnn nnnn 0000
	// N taps:       ^^^ ^^^^      = 0x07F0
	// P taps:       ^             = 0x0400
	// output: snnn nnnn NNNN NNNP
	return (cliptemp << 4)|((cliptemp&0x7F0)>>3)|((cliptemp&0x400)>>10); // upshift and range adjust
#endif
}


/**********************************************************************************************

     matrix_multiply -- does the proper multiply and shift
     a is the k coefficient and is clamped to 10 bits (9 bits plus a sign)
     b is the running result and is clamped to 14 bits.
     output is 14 bits, but note the result LSB bit is always 1.
     Because the low 4 bits of the result are trimmed off before
     output, this makes almost no difference in the computation.

**********************************************************************************************/
static INT32 matrix_multiply(INT32 a, INT32 b)
{
	INT32 result;
	while (a>511) { a-=1024; }
	while (a<-512) { a+=1024; }
	while (b>16383) { b-=32768; }
	while (b<-16384) { b+=32768; }
	result = ((a*b)>>9)|1;//&(~1);
	return result;
}

/**********************************************************************************************

     lattice_filter -- executes one 'full run' of the lattice filter on a specific byte of
     excitation data, and specific values of all the current k constants,  and returns the
     resulting sample.

***********************************************************************************************/

INT32 lattice_filter()
{
	// Lattice filter here
	// Aug/05/07: redone as unrolled loop, for clarity - LN
	/* Originally Copied verbatim from table I in US patent 4,209,804, now updated to be in same order as the actual chip does it, not that it matters.
	  notation equivalencies from table:
	  Yn(i) == m_u[n-1]
	  Kn = m_current_k[n-1]
	  bn = m_x[n-1]
	 */
	/*
	    int ep = matrix_multiply(m_previous_energy, (m_excitation_data<<6));  //Y(11)
	     m_u[10] = ep;
	    for (int i = 0; i < 10; i++)
	    {
	        int ii = 10-i; // for m = 10, this would be 11 - i, and since i is from 1 to 10, then ii ranges from 10 to 1
	        //int jj = ii+1; // this variable, even on the fortran version, is never used. it probably was intended to be used on the two lines below the next one to save some redundant additions on each.
	        ep = ep - (((m_current_k[ii-1] * m_x[ii-1])>>9)|1); // subtract reflection from lower stage 'top of lattice'
	         m_u[ii-1] = ep;
	        m_x[ii] = m_x[ii-1] + (((m_current_k[ii-1] * ep)>>9)|1); // add reflection from upper stage 'bottom of lattice'
	    }
	m_x[0] = ep; // feed the last section of the top of the lattice directly to the bottom of the lattice
	*/

		m_u[10] = matrix_multiply(m_previous_energy, (m_excitation_data<<6));  //Y(11)
		m_u[9] = m_u[10] - matrix_multiply(m_current_k[9], m_x[9]);
		m_u[8] = m_u[9] - matrix_multiply(m_current_k[8], m_x[8]);
		m_u[7] = m_u[8] - matrix_multiply(m_current_k[7], m_x[7]);
		m_u[6] = m_u[7] - matrix_multiply(m_current_k[6], m_x[6]);
		m_u[5] = m_u[6] - matrix_multiply(m_current_k[5], m_x[5]);
		m_u[4] = m_u[5] - matrix_multiply(m_current_k[4], m_x[4]);
		m_u[3] = m_u[4] - matrix_multiply(m_current_k[3], m_x[3]);
		m_u[2] = m_u[3] - matrix_multiply(m_current_k[2], m_x[2]);
		m_u[1] = m_u[2] - matrix_multiply(m_current_k[1], m_x[1]);
		m_u[0] = m_u[1] - matrix_multiply(m_current_k[0], m_x[0]);
		m_x[9] = m_x[8] + matrix_multiply(m_current_k[8], m_u[8]);
		m_x[8] = m_x[7] + matrix_multiply(m_current_k[7], m_u[7]);
		m_x[7] = m_x[6] + matrix_multiply(m_current_k[6], m_u[6]);
		m_x[6] = m_x[5] + matrix_multiply(m_current_k[5], m_u[5]);
		m_x[5] = m_x[4] + matrix_multiply(m_current_k[4], m_u[4]);
		m_x[4] = m_x[3] + matrix_multiply(m_current_k[3], m_u[3]);
		m_x[3] = m_x[2] + matrix_multiply(m_current_k[2], m_u[2]);
		m_x[2] = m_x[1] + matrix_multiply(m_current_k[1], m_u[1]);
		m_x[1] = m_x[0] + matrix_multiply(m_current_k[0], m_u[0]);
		m_x[0] = m_u[0];
		m_previous_energy = m_current_energy;
		return m_u[0];
}


// parse frame

void parse_frame()
{
	int indx, i, rep_flag;

	// We actually don't care how many bits are left in the fifo here; the frame subpart will be processed normally, and any bits extracted 'past the end' of the fifo will be read as zeroes; the fifo being emptied will set the /BE latch which will halt speech exactly as if a stop frame had been encountered (instead of whatever partial frame was read); the same exact circuitry is used for both on the real chip, see us patent 4335277 sheet 16, gates 232a (decode stop frame) and 232b (decode /BE plus DDIS (decode disable) which is active during speak external).

	/* if the chip is a tms5220C, and the rate mode is set to that each frame (0x04 bit set)
	has a 2 bit rate preceding it, grab two bits here and store them as the rate; */
	if ((TMS5220_HAS_RATE_CONTROL) && (m_c_variant_rate & 0x04))
	{
		indx = extract_bits(2);
		m_IP = reload_table[indx];
	}
	else // non-5220C and 5220C in fixed rate mode
	m_IP = reload_table[m_c_variant_rate&0x3];

	// attempt to extract the energy index
	m_new_frame_energy_idx = extract_bits(m_coeff->energy_bits);
	// if the energy index is 0 or 15, we're done
	if ((m_new_frame_energy_idx == 0) || (m_new_frame_energy_idx == 15))
		return;


	// attempt to extract the repeat flag
	rep_flag = extract_bits(1);

	// attempt to extract the pitch
	m_new_frame_pitch_idx = extract_bits(m_coeff->pitch_bits);
	// if this is a repeat frame, just do nothing, it will reuse the old coefficients
	if (rep_flag)
		return;

	// extract first 4 K coefficients
	for (i = 0; i < 4; i++)
	{
		m_new_frame_k_idx[i] = extract_bits(m_coeff->kbits[i]);

	}

	// if the pitch index was zero, we only need 4 K's...
	if (m_new_frame_pitch_idx == 0)
	{
		return;
	}

	// If we got here, we need the remaining 6 K's
	for (i = 4; i < m_coeff->num_k; i++)
	{
		m_new_frame_k_idx[i] = extract_bits(m_coeff->kbits[i]);
	}
	return;
}


////////////////bend

void parse_frame_bend_5200()
{
	int indx, i, rep_flag;

	if ((TMS5220_HAS_RATE_CONTROL) && (m_c_variant_rate & 0x04))
	{
		indx = extract_bits(2);
		m_IP = reload_table[indx];
	}
	else // non-5220C and 5220C in fixed rate mode
	m_IP = reload_table[m_c_variant_rate&0x3];

	// attempt to extract the energy index
	m_new_frame_energy_idx = extract_bits(m_coeff->energy_bits);

	// if the energy index is 0 or 15, we're done
	if ((m_new_frame_energy_idx == 0) || (m_new_frame_energy_idx == 15))
		return;

	// mult our energy bits if not 15 or 0
	m_new_frame_energy_idx*=2.0f*(1.0f-exy[_intmode][0]);
	if (m_new_frame_energy_idx>14) m_new_frame_energy_idx=14;

	// attempt to extract the repeat flag
	rep_flag = extract_bits(1);

	// attempt to extract the pitch
	m_new_frame_pitch_idx = extract_bits(m_coeff->pitch_bits);

	// if this is a repeat frame, just do nothing, it will reuse the old coefficients
	if (rep_flag)
		return;
	m_new_frame_pitch_idx*=2.0f*(1.0f-exy[_intmode][1]);
	if (m_coeff->pitch_bits==5 && m_new_frame_pitch_idx>31) m_new_frame_pitch_idx=31;
	if  (m_new_frame_pitch_idx>63) m_new_frame_pitch_idx=63; // fixed bug as was energy_idx???
	
	// extract first 4 K coefficients - deal with length of these individually
	//	{ 5, 5, 4, 4, 4, 4, 4, 3, 3, 3 },
		m_new_frame_k_idx[0] = extract_bits(m_coeff->kbits[0]);
		m_new_frame_k_idx[0]*=2.0f*exy[_intmode][2];
		if (m_new_frame_k_idx[0]>31) m_new_frame_k_idx[0]=31;

		m_new_frame_k_idx[1] = extract_bits(m_coeff->kbits[1]);
		m_new_frame_k_idx[1]*=2.0f*exy[_intmode][3];
		if (m_new_frame_k_idx[1]>31) m_new_frame_k_idx[1]=31;

		m_new_frame_k_idx[2] = extract_bits(m_coeff->kbits[2]);
		m_new_frame_k_idx[2]*=2.0f*exy[_intmode][4];
		if (m_new_frame_k_idx[2]>15) m_new_frame_k_idx[2]=15;

		m_new_frame_k_idx[3] = extract_bits(m_coeff->kbits[3]);
		m_new_frame_k_idx[3]*=2.0f*exy[_intmode][5];
		if (m_new_frame_k_idx[3]>15) m_new_frame_k_idx[3]=15;

	// if the pitch index was zero, we only need 4 K's...
	if (m_new_frame_pitch_idx == 0)
	{
		return;
	}

	// If we got here, we need the remaining 6 K's

		m_new_frame_k_idx[4] = extract_bits(m_coeff->kbits[4]);
		m_new_frame_k_idx[4]*=2.0f*exy[_intmode][6];
		if (m_new_frame_k_idx[4]>15) m_new_frame_k_idx[4]=15;

		m_new_frame_k_idx[5] = extract_bits(m_coeff->kbits[5]);
		m_new_frame_k_idx[5]*=2.0f*exy[_intmode][7];
		if (m_new_frame_k_idx[5]>15) m_new_frame_k_idx[5]=15;

		m_new_frame_k_idx[6] = extract_bits(m_coeff->kbits[6]);
		m_new_frame_k_idx[6]*=2.0f*exy[_intmode][8];
		if (m_new_frame_k_idx[6]>15) m_new_frame_k_idx[6]=15;

		m_new_frame_k_idx[7] = extract_bits(m_coeff->kbits[7]);
		m_new_frame_k_idx[7]*=2.0f*exy[_intmode][9];
		if (m_new_frame_k_idx[7]>7) m_new_frame_k_idx[7]=7;

		m_new_frame_k_idx[8] = extract_bits(m_coeff->kbits[8]);
		m_new_frame_k_idx[8]*=2.0f*exy[_intmode][10];
		if (m_new_frame_k_idx[8]>7) m_new_frame_k_idx[8]=7;

		m_new_frame_k_idx[9] = extract_bits(m_coeff->kbits[9]);
		m_new_frame_k_idx[9]*=2.0f*exy[_intmode][11];
		if (m_new_frame_k_idx[9]>7) m_new_frame_k_idx[9]=7;
	return;
}



void parse_frame_raw_5100() //- for our 3 sets of coeffs - exy is 0-10
{
  // 5100. energy bits=4=16 , pitch bits=5=32 , numk= 10 kbits= 	{ 5, 5, 4, 4, 4, 4, 4, 3, 3, 3 },
  // 5220. energy bits=4 , pitch bits=6=64 , numk= 10 kbits= 	{ 5, 5, 4, 4, 4, 4, 4, 3, 3, 3 },

  // exy is 0-1.0f 

	int indx, i, rep_flag;

	// attempt to extract the energy index
	m_new_frame_energy_idx = 13-(exy[_intmode][0]*13.0f); // exy[0]
	m_new_frame_energy_idx+=1;
	// if the energy index is 0 or 15, we're done
	//if (m_new_frame_energy_idx == 0)
	//		return;

	// attempt to extract the pitch
	//	m_new_frame_pitch_idx = extract_bits(m_coeff->pitch_bits); // exy[1] make sure it can be 0
		m_new_frame_pitch_idx = 1+( _selz*30.0f); // NO exp! - never 0 = unvoiced!!!
	//	m_new_frame_pitch_idx = 0;
	//	MAXED(m_new_frame_pitch_idx, 63);

	// extract first 4 K coefficients
	/*	for (i = 0; i < 4; i++)
	{
	  m_new_frame_k_idx[i] = extract_bits(m_coeff->kbits[i]); // exy[2]-exy[5]
	}
	*/
	m_new_frame_k_idx[0] =   exy[_intmode][1]*31.0f; //5
	m_new_frame_k_idx[1] =   exy[_intmode][2]*31.0f;
	m_new_frame_k_idx[2] =   exy[_intmode][3]*15.0f; //4
	m_new_frame_k_idx[3] =   exy[_intmode][4]*15.0f; //4

	// if the pitch index was zero, we only need 4 K's...
	if (m_new_frame_pitch_idx == 0)
	{
		/* and the rest of the coefficients are zeroed, but that's done in the generator code */
		return;
	}

	// If we got here, we need the remaining 6 K's
	/*	for (i = 4; i < m_coeff->num_k; i++)
	{
	  m_new_frame_k_idx[i] = extract_bits(m_coeff->kbits[i]); // exy[_intmode][6]-exy[_intmode][11]
	  }*/

	m_new_frame_k_idx[4] =   exy[_intmode][5]*15.0f; 
	m_new_frame_k_idx[5] =   exy[_intmode][6]*15.0f;
	m_new_frame_k_idx[6] =   exy[_intmode][7]*15.0f;
	m_new_frame_k_idx[7] =   exy[_intmode][8]*7.0f; //3
	m_new_frame_k_idx[8] =   exy[_intmode][9]*7.0f;
	m_new_frame_k_idx[9] =   exy[_intmode][10]*7.0f;
}


void parse_frame_raw_5200() //  - for our 3 sets of coeffs - exy is 0-10
{
  // 5100. energy bits=4=16 , pitch bits=5=32 , numk= 10 kbits= 	{ 5, 5, 4, 4, 4, 4, 4, 3, 3, 3 },
  // 5220. energy bits=4 , pitch bits=6=64 , numk= 10 kbits= 	{ 5, 5, 4, 4, 4, 4, 4, 3, 3, 3 },

  // exy is 0-1.0f 

	int indx, i, rep_flag;

	// attempt to extract the energy index
	m_new_frame_energy_idx = 13-(exy[_intmode][0]*13.0f); // exy[_intmode][0]
	m_new_frame_energy_idx+=1;
	// if the energy index is 0 or 15, we're done
	//	if ((m_new_frame_energy_idx == 0) || (m_new_frame_energy_idx == 15))
	//		return;

	// attempt to extract the pitch
	//	m_new_frame_pitch_idx = extract_bits(m_coeff->pitch_bits); // exy[_intmode][1] make sure it can be 0
	m_new_frame_pitch_idx = 1+( _selz*62.0f);
	//	if (m_new_frame_pitch_idx>=60) m_new_frame_pitch_idx=0;
	//	MAXED(m_new_frame_pitch_idx, 63);

	// extract first 4 K coefficients
	m_new_frame_k_idx[0] =   exy[_intmode][1]*31.0f;
	m_new_frame_k_idx[1] =   exy[_intmode][2]*31.0f;
	m_new_frame_k_idx[2] =   exy[_intmode][3]*15.0f;
	m_new_frame_k_idx[3] =   exy[_intmode][4]*15.0f;

	// if the pitch index was zero, we only need 4 K's...
	if (m_new_frame_pitch_idx == 0)
	{
		return;
	}
	
	m_new_frame_k_idx[4] =   exy[_intmode][5]*15.0f;
	m_new_frame_k_idx[5] =   exy[_intmode][6]*15.0f;
	m_new_frame_k_idx[6] =   exy[_intmode][7]*15.0f;
	m_new_frame_k_idx[7] =   exy[_intmode][8]*7.0f;
	m_new_frame_k_idx[8] =   exy[_intmode][9]*7.0f;
	m_new_frame_k_idx[9] =   exy[_intmode][10]*7.0f;
}

int16_t process(u8 *ending)
{
  static u8 counter;
	int i, bitout, zpar;
	INT32 this_sample;
	int16_t sample, val;
		if ((m_IP == 0) && (m_PC == 0) && (m_subcycle < 2))
		{
			m_OLDE = (m_new_frame_energy_idx == 0);
			m_OLDP = (m_new_frame_pitch_idx == 0);
		}
		if ((m_IP == 0) && (m_PC == 12) && (m_subcycle == 1))
		{
			m_IP = reload_table[m_c_variant_rate&0x3];

#ifdef PERFECT_INTERPOLATION_HACK
			/* remember previous frame energy, pitch, and coefficients */
			m_old_frame_energy_idx = m_new_frame_energy_idx;
			m_old_frame_pitch_idx = m_new_frame_pitch_idx;
			for (i = 0; i < m_coeff->num_k; i++)
				m_old_frame_k_idx[i] = m_new_frame_k_idx[i];
#endif

			/* if the talk status was clear last frame, halt speech now. */
			if (m_talk_status == 0)
			{
				//				m_speaking_now = 1; // finally halt speech
				*ending=1;
				return 0;
			}


			/* Parse a new frame into the new_target_energy, new_target_pitch and new_target_k[] */
			parse_frame();
			/* if the new frame is a stop frame, set an interrupt and set talk status to 0 */
			if (NEW_FRAME_STOP_FLAG == 1)
				{
					m_talk_status = m_speak_external = 0;
					//					set_interrupt_state(1);
					//					update_fifo_status_and_ints();
				}

			/* in all cases where interpolation would be inhibited, set the inhibit flag; otherwise clear it.
			   Interpolation inhibit cases:
			 * Old frame was voiced, new is unvoiced
			 * Old frame was silence/zero energy, new has nonzero energy
			 * Old frame was unvoiced, new is voiced (note this is the case on the patent but may not be correct on the real final chip)
			 */
			if ( ((OLD_FRAME_UNVOICED_FLAG == 0) && (NEW_FRAME_UNVOICED_FLAG == 1))
				|| ((OLD_FRAME_UNVOICED_FLAG == 1) && (NEW_FRAME_UNVOICED_FLAG == 0)) /* this line needs further investigation, starwars tie fighters may sound better without it */
				|| ((OLD_FRAME_SILENCE_FLAG == 1) && (NEW_FRAME_SILENCE_FLAG == 0)) )
				m_inhibit = 1;
			else // normal frame, normal interpolation
				m_inhibit = 0;

			/* load new frame targets from tables, using parsed indices */
			m_target_energy = m_coeff->energytable[m_new_frame_energy_idx];
			m_target_pitch = m_coeff->pitchtable[m_new_frame_pitch_idx];
			zpar = NEW_FRAME_UNVOICED_FLAG; // find out if parameters k5-k10 should be zeroed
			for (i = 0; i < 4; i++)
				m_target_k[i] = m_coeff->ktable[i][m_new_frame_k_idx[i]];
			for (i = 4; i < m_coeff->num_k; i++)
				m_target_k[i] = (m_coeff->ktable[i][m_new_frame_k_idx[i]] * (1-zpar));

			/* if TS is now 0, ramp the energy down to 0. Is this really correct to hardware? */
			if (m_talk_status == 0)
			{
			  m_target_energy = 0;
			  //			  for (i = 0; i < 4; i++) m_target_k[i] = 0;
			  //			  for (i = 4; i < m_coeff->num_k; i++)	m_target_k[i] = 0;
			}
		}
		else // Not a new frame, just interpolate the existing frame.
		{
			int inhibit_state = ((m_inhibit==1)&&(m_IP != 0)); // disable inhibit when reaching the last interp period, but don't overwrite the m_inhibit value
#ifdef PERFECT_INTERPOLATION_HACK
			int samples_per_frame = m_subc_reload?175:266; // either (13 A cycles + 12 B cycles) * 7 interps for normal SPEAK/SPKEXT, or (13*2 A cycles + 12 B cycles) * 7 interps for SPKSLOW
			//int samples_per_frame = m_subc_reload?200:304; // either (13 A cycles + 12 B cycles) * 8 interps for normal SPEAK/SPKEXT, or (13*2 A cycles + 12 B cycles) * 8 interps for SPKSLOW
			// try changing this up/down
			
			int current_sample = (m_subcycle - m_subc_reload)+(m_PC*(3-m_subc_reload))+((m_subc_reload?25:38)*((m_IP-1)&7));

			zpar = OLD_FRAME_UNVOICED_FLAG;
			//fprintf(stderr, "CS: %03d", current_sample);
			// reset the current energy, pitch, etc to what it was at frame start
			m_current_energy = m_coeff->energytable[m_old_frame_energy_idx];
			m_current_pitch = m_coeff->pitchtable[m_old_frame_pitch_idx];
			for (i = 0; i < 4; i++)
				m_current_k[i] = m_coeff->ktable[i][m_old_frame_k_idx[i]];
			for (i = 4; i < m_coeff->num_k; i++)
				m_current_k[i] = (m_coeff->ktable[i][m_old_frame_k_idx[i]] * (1-zpar));
			// now adjust each value to be exactly correct for each of the samples per frame
			if (m_IP != 0) // if we're still interpolating...
			{
				m_current_energy += (((m_target_energy - m_current_energy)*(1-inhibit_state))*current_sample)/samples_per_frame;
				m_current_pitch += (((m_target_pitch - m_current_pitch)*(1-inhibit_state))*current_sample)/samples_per_frame;
				for (i = 0; i < m_coeff->num_k; i++)
					m_current_k[i] += (((m_target_k[i] - m_current_k[i])*(1-inhibit_state))*current_sample)/samples_per_frame;
			}
			else // we're done, play this frame for 1/8 frame.
			{
				m_current_energy = m_target_energy;
				m_current_pitch = m_target_pitch;
				for (i = 0; i < m_coeff->num_k; i++)
					m_current_k[i] = m_target_k[i];
			}
#else
			//Updates to parameters only happen on subcycle '2' (B cycle) of PCs.
			if (m_subcycle == 2)
			{
				switch(m_PC)
				{
					case 0: /* PC = 0, B cycle, write updated energy */
					m_current_energy += (((m_target_energy - m_current_energy)*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 1: /* PC = 1, B cycle, write updated pitch */
					m_current_pitch += (((m_target_pitch - m_current_pitch)*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 9: case 10: case 11:
					/* PC = 2 through 11, B cycle, write updated K1 through K10 */
					m_current_k[m_PC-2] += (((m_target_k[m_PC-2] - m_current_k[m_PC-2])*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 12: /* PC = 12, do nothing */
					break;
				}
			}
#endif
		}

		// calculate the output
		if (OLD_FRAME_UNVOICED_FLAG == 1)
		{
			// generate unvoiced samples here
			if (m_RNG & 1)
				m_excitation_data = ~0x3F; /* according to the patent it is (either + or -) half of the maximum value in the chirp table, so either 01000000(0x40) or 11000000(0xC0)*/
			else
			  m_excitation_data = 0x40; // these values can be start of bend
		}
		else /* (OLD_FRAME_UNVOICED_FLAG == 0) */
		{
			// generate voiced samples here
			/* US patent 4331836 Figure 14B shows, and logic would hold, that a pitch based chirp
			 * function has a chirp/peak and then a long chain of zeroes.
			 * The last entry of the chirp rom is at address 0b110011 (51d), the 52nd sample,
			 * and if the address reaches that point the ADDRESS incrementer is
			 * disabled, forcing all samples beyond 51d to be == 51d
			 */
			if (m_pitch_count >= 51)
				m_excitation_data = (INT8)m_coeff->chirptable[51];
			else /*m_pitch_count < 51*/
				m_excitation_data = (INT8)m_coeff->chirptable[m_pitch_count];
		}

		// Update LFSR *20* times every sample (once per T cycle), like patent shows
	for (i=0; i<20; i++)
	{
		bitout = ((m_RNG >> 12) & 1) ^
				((m_RNG >>  3) & 1) ^
				((m_RNG >>  2) & 1) ^
				((m_RNG >>  0) & 1);
		m_RNG <<= 1;
		m_RNG |= bitout;
	}
	/////////////////////////////

	//			this_sample = lattice_filter(); /* execute lattice filter */
			this_sample=m_excitation_data<<4;

		/* next, force result to 14 bits (since its possible that the addition at the final (k1) stage of the lattice overflowed) */
		while (this_sample > 16383) this_sample -= 32768;
		while (this_sample < -16384) this_sample += 32768;
		if (m_digital_select == 0) // analog SPK pin output is only 8 bits, with clipping
		  //			buffer[buf_count] = clip_analog(this_sample);
		  sample= clip_analog(this_sample);
		else // digital I/O pin output is 12 bits
		{
#ifdef ALLOW_4_LSB
			// input:  ssss ssss ssss ssss ssnn nnnn nnnn nnnn
			// N taps:                       ^                 = 0x2000;
			// output: ssss ssss ssss ssss snnn nnnn nnnn nnnN
		  //			buffer[buf_count] = (this_sample<<1)|((this_sample&0x2000)>>13);
		  sample=(this_sample<<1)|((this_sample&0x2000)>>13);
#else
			this_sample &= ~0xF;
			// input:  ssss ssss ssss ssss ssnn nnnn nnnn 0000
			// N taps:                       ^^ ^^^            = 0x3E00;
			// output: ssss ssss ssss ssss snnn nnnn nnnN NNNN
			//			buffer[buf_count] = (this_sample<<1)|((this_sample&0x3E00)>>9);
			sample=(this_sample<<1)|((this_sample&0x3E00)>>9);
#endif
		}
		// Update all counts

		m_subcycle++;
		if ((m_subcycle == 2) && (m_PC == 12))
		{
			/* Circuit 412 in the patent acts a reset, resetting the pitch counter to 0
			 * if INHIBIT was true during the most recent frame transition.
			 * The exact time this occurs is betwen IP=7, PC=12 sub=0, T=t12
			 * and m_IP = 0, PC=0 sub=0, T=t12, a period of exactly 20 cycles,
			 * which overlaps the time OLDE and OLDP are updated at IP=7 PC=12 T17
			 * (and hence INHIBIT itself 2 t-cycles later). We do it here because it is
			 * convenient and should make no difference in output.
			 */
			if ((m_IP == 7)&&(m_inhibit==1)) m_pitch_count = 0;
			m_subcycle = m_subc_reload;
			m_PC = 0;
			m_IP++;
			m_IP&=0x7;
		}
		else if (m_subcycle == 3)
		{
			m_subcycle = m_subc_reload;
			if (modus&4){
			if (counter++>=(16-((1.0f-_selx)*17.0f))){
			  m_PC++;
			counter=0;
			}
			}
			else if (modus&8){
			if (counter++>=(16-((1.0f-_sely)*17.0f))){
			  m_PC++;
			counter=0;
			}
			}
			else m_PC++;
		}
		m_pitch_count++;
		val=_selz*131.0f;
		MAXED(val,127);
		val=127-val;
		if (modus&1){
		  if (m_pitch_count >= (m_current_pitch*logpitch[val])) m_pitch_count = 0; //  - pitch bend
		  //		  		  if (m_pitch_count >= (m_current_pitch)) m_pitch_count = 0; 
		}
		else if (modus&2)
		  {
		    if (m_pitch_count >= (64.0f*logpitch[val])) m_pitch_count = 0; //  - pitch bend
		  }
		else
		   if (m_pitch_count >= m_current_pitch) m_pitch_count = 0; //
		  
		m_pitch_count &= 0x1FF;
		if (m_digital_select==0) return sample>>1; // 0 is lowbit
		else
		  return sample<<1; // tested and leave as is 
}

int16_t process_k_tabled5100(u8 *ending) 
{
  static u8 counter;
	int i, bitout, zpar;
	INT32 this_sample;
	int16_t sample, val;
		if ((m_IP == 0) && (m_PC == 0) && (m_subcycle < 2))
		{
			m_OLDE = (m_new_frame_energy_idx == 0);
			m_OLDP = (m_new_frame_pitch_idx == 0);
		}
		if ((m_IP == 0) && (m_PC == 12) && (m_subcycle == 1))
		{
			m_IP = reload_table[m_c_variant_rate&0x3];

			if (m_talk_status == 0)
			{
				*ending=1;
				return 0;
			}
			parse_frame(); 
			if (NEW_FRAME_STOP_FLAG == 1)
				{
					m_talk_status = m_speak_external = 0;
				}

			if ( ((OLD_FRAME_UNVOICED_FLAG == 0) && (NEW_FRAME_UNVOICED_FLAG == 1))
				|| ((OLD_FRAME_UNVOICED_FLAG == 1) && (NEW_FRAME_UNVOICED_FLAG == 0)) /* this line needs further investigation, starwars tie fighters may sound better without it */
				|| ((OLD_FRAME_SILENCE_FLAG == 1) && (NEW_FRAME_SILENCE_FLAG == 0)) )
				m_inhibit = 1;
			else // normal frame, normal interpolation
				m_inhibit = 0;

			/* load new frame targets from tables, using parsed indices */
			m_target_energy = m_coeff->energytable[m_new_frame_energy_idx];
			m_target_pitch = m_coeff->pitchtable[m_new_frame_pitch_idx];
			zpar = NEW_FRAME_UNVOICED_FLAG; // find out if parameters k5-k10 should be zeroed
			// 168 in total
			for (i = 0; i < 4; i++)
				m_target_k[i] = m_coeff->ktable[i][m_new_frame_k_idx[i]]*(2.0f*(exy[_intmode][koffset[i]+m_new_frame_k_idx[i]]+0.1f));
			for (i = 4; i < m_coeff->num_k; i++)
			  m_target_k[i] = (m_coeff->ktable[i][m_new_frame_k_idx[i]] * (1-zpar)*(2.0f*(exy[_intmode][koffset[i]+m_new_frame_k_idx[i]]+0.1f))); // ktable size is 168

			if (m_talk_status == 0)
			{
				m_target_energy = 0;
			}
		}
		else // Not a new frame, just interpolate the existing frame.
		{
			int inhibit_state = ((m_inhibit==1)&&(m_IP != 0)); // disable inhibit when reaching the last interp period, but don't overwrite the m_inhibit value
			//Updates to parameters only happen on subcycle '2' (B cycle) of PCs.
			if (m_subcycle == 2)
			{
				switch(m_PC)
				{
					case 0: /* PC = 0, B cycle, write updated energy */
					m_current_energy += (((m_target_energy - m_current_energy)*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 1: /* PC = 1, B cycle, write updated pitch */
					m_current_pitch += (((m_target_pitch - m_current_pitch)*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 9: case 10: case 11:
					/* PC = 2 through 11, B cycle, write updated K1 through K10 */
					m_current_k[m_PC-2] += (((m_target_k[m_PC-2] - m_current_k[m_PC-2])*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 12: /* PC = 12, do nothing */
					break;
				}
			}
		}

		// calculate the output
		if (OLD_FRAME_UNVOICED_FLAG == 1)
		{
			// generate unvoiced samples here
			if (m_RNG & 1)
				m_excitation_data = ~0x3F; /* according to the patent it is (either + or -) half of the maximum value in the chirp table, so either 01000000(0x40) or 11000000(0xC0)*/
			else
				m_excitation_data = 0x40;
		}
		else /* (OLD_FRAME_UNVOICED_FLAG == 0) */
		{
			if (m_pitch_count >= 51)
			  m_excitation_data = (INT8)m_coeff->chirptable[51];//*(2.0f*(exy[_intmode][51]+0.1f));
			else /*m_pitch_count < 51*/
			  m_excitation_data = (INT8)m_coeff->chirptable[m_pitch_count];//*(2.0f*(exy[_intmode][m_pitch_count]+0.1f));
		}

		// Update LFSR *20* times every sample (once per T cycle), like patent shows
	for (i=0; i<20; i++)
	{
		bitout = ((m_RNG >> 12) & 1) ^
				((m_RNG >>  3) & 1) ^
				((m_RNG >>  2) & 1) ^
				((m_RNG >>  0) & 1);
		m_RNG <<= 1;
		m_RNG |= bitout;
	}
	/////////////////////////////

		this_sample = lattice_filter(); /* execute lattice filter */

		while (this_sample > 16383) this_sample -= 32768;
		while (this_sample < -16384) this_sample += 32768;
		if (m_digital_select == 0) // analog SPK pin output is only 8 bits, with clipping
		  //			buffer[buf_count] = clip_analog(this_sample);
		  sample= clip_analog(this_sample);
		else // digital I/O pin output is 12 bits
		{
#ifdef ALLOW_4_LSB
		  sample=(this_sample<<1)|((this_sample&0x2000)>>13);
#else
			this_sample &= ~0xF;
			sample=(this_sample<<1)|((this_sample&0x3E00)>>9);
#endif
		}

		m_subcycle++;
		if ((m_subcycle == 2) && (m_PC == 12))
		{
			if ((m_IP == 7)&&(m_inhibit==1)) m_pitch_count = 0;
			m_subcycle = m_subc_reload;
			m_PC = 0;
			m_IP++;
			m_IP&=0x7;
		}
		else if (m_subcycle == 3)
		{
			m_subcycle = m_subc_reload;
			m_PC++;
		}
		m_pitch_count++;
		if (m_pitch_count >= m_current_pitch) m_pitch_count = 0; 
		m_pitch_count &= 0x1FF;
	return sample;
}


int16_t process_pitchk_tabled5100(u8 *ending) // for 5100 we have 32+168 in exy= 200 5200 is 232
{
  static u8 counter;
  u8 offsetage;
  int i, bitout, zpar;
	INT32 this_sample;
	int16_t sample, val;
		if ((m_IP == 0) && (m_PC == 0) && (m_subcycle < 2))
		{
			m_OLDE = (m_new_frame_energy_idx == 0);
			m_OLDP = (m_new_frame_pitch_idx == 0);
		}
		if ((m_IP == 0) && (m_PC == 12) && (m_subcycle == 1))
		{
			m_IP = reload_table[m_c_variant_rate&0x3];

			if (m_talk_status == 0)
			{
				*ending=1;
				return 0;
			}
			parse_frame(); 
			if (NEW_FRAME_STOP_FLAG == 1)
				{
					m_talk_status = m_speak_external = 0;
				}

			if ( ((OLD_FRAME_UNVOICED_FLAG == 0) && (NEW_FRAME_UNVOICED_FLAG == 1))
				|| ((OLD_FRAME_UNVOICED_FLAG == 1) && (NEW_FRAME_UNVOICED_FLAG == 0)) /* this line needs further investigation, starwars tie fighters may sound better without it */
				|| ((OLD_FRAME_SILENCE_FLAG == 1) && (NEW_FRAME_SILENCE_FLAG == 0)) )
				m_inhibit = 1;
			else // normal frame, normal interpolation
				m_inhibit = 0;

			/* load new frame targets from tables, using parsed indices */
			m_target_energy = m_coeff->energytable[m_new_frame_energy_idx];
			u8 val=exy[_intmode][m_new_frame_pitch_idx]*131.0f;
			MAXED(val,127);
			m_target_pitch = m_coeff->pitchtable[m_new_frame_pitch_idx]*logpitch[val];

			offsetage=1<<m_coeff->pitch_bits; // 32 in case of 5100 or 64 for others
			
			
			for (i = 0; i < 4; i++) {
			  //			  val=exy[_intmode][offsetage+(koffset[i]+m_new_frame_k_idx[i])]*131.0f;
			  //			MAXED(val,127);
			  m_target_k[i] = m_coeff->ktable[i][m_new_frame_k_idx[i]]*(2.0f*exy[_intmode][offsetage+(koffset[i]+m_new_frame_k_idx[i])]+0.1f);
			  //			m_target_k[i] = m_coeff->ktable[i][m_new_frame_k_idx[i]]*logpitch[val];
			}
			for (i = 4; i < m_coeff->num_k; i++){
			  //			  val=exy[_intmode][offsetage+(koffset[i]+m_new_frame_k_idx[i])]*131.0f;
			  //			  MAXED(val,127);
			  //			  m_target_k[i] = m_coeff->ktable[i][m_new_frame_k_idx[i]]*(1-zpar) * logpitch[val];


			  m_target_k[i] = (m_coeff->ktable[i][m_new_frame_k_idx[i]] * (1-zpar)*(2.0f*exy[_intmode][offsetage+(koffset[i]+m_new_frame_k_idx[i])]+0.1f));
			}

			if (m_talk_status == 0)
			{
				m_target_energy = 0;
			}
		}
		else // Not a new frame, just interpolate the existing frame.
		{
			int inhibit_state = ((m_inhibit==1)&&(m_IP != 0)); // disable inhibit when reaching the last interp period, but don't overwrite the m_inhibit value
			//Updates to parameters only happen on subcycle '2' (B cycle) of PCs.
			if (m_subcycle == 2)
			{
				switch(m_PC)
				{
					case 0: /* PC = 0, B cycle, write updated energy */
					m_current_energy += (((m_target_energy - m_current_energy)*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 1: /* PC = 1, B cycle, write updated pitch */
					m_current_pitch += (((m_target_pitch - m_current_pitch)*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 9: case 10: case 11:
					/* PC = 2 through 11, B cycle, write updated K1 through K10 */
					m_current_k[m_PC-2] += (((m_target_k[m_PC-2] - m_current_k[m_PC-2])*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 12: /* PC = 12, do nothing */
					break;
				}
			}
		}

		// calculate the output
		if (OLD_FRAME_UNVOICED_FLAG == 1)
		{
			// generate unvoiced samples here
			if (m_RNG & 1)
				m_excitation_data = ~0x3F; /* according to the patent it is (either + or -) half of the maximum value in the chirp table, so either 01000000(0x40) or 11000000(0xC0)*/
			else
				m_excitation_data = 0x40;
		}
		else /* (OLD_FRAME_UNVOICED_FLAG == 0) */
		{
			if (m_pitch_count >= 51)
			  m_excitation_data = (INT8)m_coeff->chirptable[51];//*(2.0f*(exy[_intmode][51]+0.1f));
			else /*m_pitch_count < 51*/
			  m_excitation_data = (INT8)m_coeff->chirptable[m_pitch_count];//*(2.0f*(exy[_intmode][m_pitch_count]+0.1f));
		}

		// Update LFSR *20* times every sample (once per T cycle), like patent shows
	for (i=0; i<20; i++)
	{
		bitout = ((m_RNG >> 12) & 1) ^
				((m_RNG >>  3) & 1) ^
				((m_RNG >>  2) & 1) ^
				((m_RNG >>  0) & 1);
		m_RNG <<= 1;
		m_RNG |= bitout;
	}
	/////////////////////////////

		this_sample = lattice_filter(); /* execute lattice filter */

		while (this_sample > 16383) this_sample -= 32768;
		while (this_sample < -16384) this_sample += 32768;
		if (m_digital_select == 0) // analog SPK pin output is only 8 bits, with clipping
		  //			buffer[buf_count] = clip_analog(this_sample);
		  sample= clip_analog(this_sample);
		else // digital I/O pin output is 12 bits
		{
#ifdef ALLOW_4_LSB
		  sample=(this_sample<<1)|((this_sample&0x2000)>>13);
#else
			this_sample &= ~0xF;
			sample=(this_sample<<1)|((this_sample&0x3E00)>>9);
#endif
		}

		m_subcycle++;
		if ((m_subcycle == 2) && (m_PC == 12))
		{
			if ((m_IP == 7)&&(m_inhibit==1)) m_pitch_count = 0;
			m_subcycle = m_subc_reload;
			m_PC = 0;
			m_IP++;
			m_IP&=0x7;
		}
		else if (m_subcycle == 3)
		{
			m_subcycle = m_subc_reload;
			m_PC++;
		}
		m_pitch_count++;
		if (m_pitch_count >= m_current_pitch) m_pitch_count = 0; // TEST - pitch bend
		m_pitch_count &= 0x1FF;
	return sample;
}

int16_t process_pitch_tabled5100(u8 *ending)  // also stripped down
{
  static u8 counter;
	int i, bitout, zpar;
	INT32 this_sample;
	int16_t sample, val;
		if ((m_IP == 0) && (m_PC == 0) && (m_subcycle < 2))
		{
			m_OLDE = (m_new_frame_energy_idx == 0);
			m_OLDP = (m_new_frame_pitch_idx == 0);
		}
		if ((m_IP == 0) && (m_PC == 12) && (m_subcycle == 1))
		{
			m_IP = reload_table[m_c_variant_rate&0x3];

			if (m_talk_status == 0)
			{
				*ending=1;
				return 0;
			}
			parse_frame(); 
			if (NEW_FRAME_STOP_FLAG == 1)
				{
					m_talk_status = m_speak_external = 0;
				}

			if ( ((OLD_FRAME_UNVOICED_FLAG == 0) && (NEW_FRAME_UNVOICED_FLAG == 1))
				|| ((OLD_FRAME_UNVOICED_FLAG == 1) && (NEW_FRAME_UNVOICED_FLAG == 0)) /* this line needs further investigation, starwars tie fighters may sound better without it */
				|| ((OLD_FRAME_SILENCE_FLAG == 1) && (NEW_FRAME_SILENCE_FLAG == 0)) )
				m_inhibit = 1;
			else // normal frame, normal interpolation
				m_inhibit = 0;

			/* load new frame targets from tables, using parsed indices */
			m_target_energy = m_coeff->energytable[m_new_frame_energy_idx];
			u8 val=exy[_intmode][m_new_frame_pitch_idx]*131.0f;
			MAXED(val,127);
			m_target_pitch = m_coeff->pitchtable[m_new_frame_pitch_idx]*logpitch[val];

			zpar = NEW_FRAME_UNVOICED_FLAG; // find out if parameters k5-k10 should be zeroed
			for (i = 0; i < 4; i++)
				m_target_k[i] = m_coeff->ktable[i][m_new_frame_k_idx[i]];
			for (i = 4; i < m_coeff->num_k; i++)
				m_target_k[i] = (m_coeff->ktable[i][m_new_frame_k_idx[i]] * (1-zpar));

			if (m_talk_status == 0)
			{
				m_target_energy = 0;
			}
		}
		else // Not a new frame, just interpolate the existing frame.
		{
			int inhibit_state = ((m_inhibit==1)&&(m_IP != 0)); // disable inhibit when reaching the last interp period, but don't overwrite the m_inhibit value
			//Updates to parameters only happen on subcycle '2' (B cycle) of PCs.
			if (m_subcycle == 2)
			{
				switch(m_PC)
				{
					case 0: /* PC = 0, B cycle, write updated energy */
					m_current_energy += (((m_target_energy - m_current_energy)*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 1: /* PC = 1, B cycle, write updated pitch */
					m_current_pitch += (((m_target_pitch - m_current_pitch)*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 9: case 10: case 11:
					/* PC = 2 through 11, B cycle, write updated K1 through K10 */
					m_current_k[m_PC-2] += (((m_target_k[m_PC-2] - m_current_k[m_PC-2])*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 12: /* PC = 12, do nothing */
					break;
				}
			}
		}

		// calculate the output
		if (OLD_FRAME_UNVOICED_FLAG == 1)
		{
			// generate unvoiced samples here
			if (m_RNG & 1)
				m_excitation_data = ~0x3F; /* according to the patent it is (either + or -) half of the maximum value in the chirp table, so either 01000000(0x40) or 11000000(0xC0)*/
			else
				m_excitation_data = 0x40;
		}
		else /* (OLD_FRAME_UNVOICED_FLAG == 0) */
		{
			if (m_pitch_count >= 51)
				m_excitation_data = (INT8)m_coeff->chirptable[51];
			else /*m_pitch_count < 51*/
				m_excitation_data = (INT8)m_coeff->chirptable[m_pitch_count];
		}

		// Update LFSR *20* times every sample (once per T cycle), like patent shows
	for (i=0; i<20; i++)
	{
		bitout = ((m_RNG >> 12) & 1) ^
				((m_RNG >>  3) & 1) ^
				((m_RNG >>  2) & 1) ^
				((m_RNG >>  0) & 1);
		m_RNG <<= 1;
		m_RNG |= bitout;
	}
	/////////////////////////////

		this_sample = lattice_filter(); /* execute lattice filter */

		while (this_sample > 16383) this_sample -= 32768;
		while (this_sample < -16384) this_sample += 32768;
		if (m_digital_select == 0) // analog SPK pin output is only 8 bits, with clipping
		  //			buffer[buf_count] = clip_analog(this_sample);
		  sample= clip_analog(this_sample);
		else // digital I/O pin output is 12 bits
		{
#ifdef ALLOW_4_LSB
		  sample=(this_sample<<1)|((this_sample&0x2000)>>13);
#else
			this_sample &= ~0xF;
			sample=(this_sample<<1)|((this_sample&0x3E00)>>9);
#endif
		}

		m_subcycle++;
		if ((m_subcycle == 2) && (m_PC == 12))
		{
			if ((m_IP == 7)&&(m_inhibit==1)) m_pitch_count = 0;
			m_subcycle = m_subc_reload;
			m_PC = 0;
			m_IP++;
			m_IP&=0x7;
		}
		else if (m_subcycle == 3)
		{
			m_subcycle = m_subc_reload;
			m_PC++;
		}
		m_pitch_count++;
		if (m_pitch_count >= m_current_pitch) m_pitch_count = 0; 
		m_pitch_count &= 0x1FF;
	return sample;
}


int16_t processbend5100(u8 *ending)
{
  static u8 counter;
	int i, bitout, zpar;
	INT32 this_sample;
	int16_t sample, val;
	/* loop until the buffer is full or we've stopped speaking */
	//	if (m_speaking_now)
	//	{
		/* if it is the appropriate time to update the old energy/pitch idxes,
		 * i.e. when IP=7, PC=12, T=17, subcycle=2, do so. Since IP=7 PC=12 T=17
		 * is JUST BEFORE the transition to IP=0 PC=0 T=0 sybcycle=(0 or 1),
		 * which happens 4 T-cycles later), we change on the latter.*/
		if ((m_IP == 0) && (m_PC == 0) && (m_subcycle < 2))
		{
			m_OLDE = (m_new_frame_energy_idx == 0);
			m_OLDP = (m_new_frame_pitch_idx == 0);
		}

		/* if we're ready for a new frame to be applied, i.e. when IP=0, PC=12, Sub=1
		 * (In reality, the frame was really loaded incrementally during the entire IP=0
		 * PC=x time period, but it doesn't affect anything until IP=0 PC=12 happens)
		 */
		if ((m_IP == 0) && (m_PC == 12) && (m_subcycle == 1))
		{
			// HACK for regression testing, be sure to comment out before release!
			//m_RNG = 0x1234;
			// end HACK

			/* appropriately override the interp count if needed; this will be incremented after the frame parse! */
			m_IP = reload_table[m_c_variant_rate&0x3];

#ifdef PERFECT_INTERPOLATION_HACK
			/* remember previous frame energy, pitch, and coefficients */
			m_old_frame_energy_idx = m_new_frame_energy_idx;
			m_old_frame_pitch_idx = m_new_frame_pitch_idx;
			for (i = 0; i < m_coeff->num_k; i++)
				m_old_frame_k_idx[i] = m_new_frame_k_idx[i];
#endif

			/* if the talk status was clear last frame, halt speech now. */
			if (m_talk_status == 0)
			{
				//				m_speaking_now = 1; // finally halt speech
				*ending=1;
				// keep speaking - RESET TOD!
				//				goto empty;
				return 0;
			}


			/* Parse a new frame into the new_target_energy, new_target_pitch and new_target_k[] */
			//			if (m_coeff->pitch_bits==5)			parse_frame_bend_5200();
			parse_frame_bend_5200();
			/* if the new frame is a stop frame, set an interrupt and set talk status to 0 */
			if (NEW_FRAME_STOP_FLAG == 1)
				{
					m_talk_status = m_speak_external = 0;
					//					set_interrupt_state(1);
					//					update_fifo_status_and_ints();
				}

			/* in all cases where interpolation would be inhibited, set the inhibit flag; otherwise clear it.
			   Interpolation inhibit cases:
			 * Old frame was voiced, new is unvoiced
			 * Old frame was silence/zero energy, new has nonzero energy
			 * Old frame was unvoiced, new is voiced (note this is the case on the patent but may not be correct on the real final chip)
			 */
			if ( ((OLD_FRAME_UNVOICED_FLAG == 0) && (NEW_FRAME_UNVOICED_FLAG == 1))
				|| ((OLD_FRAME_UNVOICED_FLAG == 1) && (NEW_FRAME_UNVOICED_FLAG == 0)) /* this line needs further investigation, starwars tie fighters may sound better without it */
				|| ((OLD_FRAME_SILENCE_FLAG == 1) && (NEW_FRAME_SILENCE_FLAG == 0)) )
				m_inhibit = 1;
			else // normal frame, normal interpolation
				m_inhibit = 0;

			/* load new frame targets from tables, using parsed indices */
			m_target_energy = m_coeff->energytable[m_new_frame_energy_idx];
			m_target_pitch = m_coeff->pitchtable[m_new_frame_pitch_idx];
			zpar = NEW_FRAME_UNVOICED_FLAG; // find out if parameters k5-k10 should be zeroed
			for (i = 0; i < 4; i++)
				m_target_k[i] = m_coeff->ktable[i][m_new_frame_k_idx[i]];
			for (i = 4; i < m_coeff->num_k; i++)
				m_target_k[i] = (m_coeff->ktable[i][m_new_frame_k_idx[i]] * (1-zpar));

			/* if TS is now 0, ramp the energy down to 0. Is this really correct to hardware? */
			if (m_talk_status == 0)
			{
				m_target_energy = 0;
			}
		}
		else // Not a new frame, just interpolate the existing frame.
		{
			int inhibit_state = ((m_inhibit==1)&&(m_IP != 0)); // disable inhibit when reaching the last interp period, but don't overwrite the m_inhibit value
#ifdef PERFECT_INTERPOLATION_HACK
			int samples_per_frame = m_subc_reload?175:266; // either (13 A cycles + 12 B cycles) * 7 interps for normal SPEAK/SPKEXT, or (13*2 A cycles + 12 B cycles) * 7 interps for SPKSLOW
			//int samples_per_frame = m_subc_reload?200:304; // either (13 A cycles + 12 B cycles) * 8 interps for normal SPEAK/SPKEXT, or (13*2 A cycles + 12 B cycles) * 8 interps for SPKSLOW
			// try changing this up/down
			
			int current_sample = (m_subcycle - m_subc_reload)+(m_PC*(3-m_subc_reload))+((m_subc_reload?25:38)*((m_IP-1)&7));

			zpar = OLD_FRAME_UNVOICED_FLAG;
			//fprintf(stderr, "CS: %03d", current_sample);
			// reset the current energy, pitch, etc to what it was at frame start
			m_current_energy = m_coeff->energytable[m_old_frame_energy_idx];
			m_current_pitch = m_coeff->pitchtable[m_old_frame_pitch_idx];
			for (i = 0; i < 4; i++)
				m_current_k[i] = m_coeff->ktable[i][m_old_frame_k_idx[i]];
			for (i = 4; i < m_coeff->num_k; i++)
				m_current_k[i] = (m_coeff->ktable[i][m_old_frame_k_idx[i]] * (1-zpar));
			// now adjust each value to be exactly correct for each of the samples per frame
			if (m_IP != 0) // if we're still interpolating...
			{
				m_current_energy += (((m_target_energy - m_current_energy)*(1-inhibit_state))*current_sample)/samples_per_frame;
				m_current_pitch += (((m_target_pitch - m_current_pitch)*(1-inhibit_state))*current_sample)/samples_per_frame;
				for (i = 0; i < m_coeff->num_k; i++)
					m_current_k[i] += (((m_target_k[i] - m_current_k[i])*(1-inhibit_state))*current_sample)/samples_per_frame;
			}
			else // we're done, play this frame for 1/8 frame.
			{
				m_current_energy = m_target_energy;
				m_current_pitch = m_target_pitch;
				for (i = 0; i < m_coeff->num_k; i++)
					m_current_k[i] = m_target_k[i];
			}
#else
			//Updates to parameters only happen on subcycle '2' (B cycle) of PCs.
			if (m_subcycle == 2)
			{
				switch(m_PC)
				{
					case 0: /* PC = 0, B cycle, write updated energy */
					m_current_energy += (((m_target_energy - m_current_energy)*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 1: /* PC = 1, B cycle, write updated pitch */
					m_current_pitch += (((m_target_pitch - m_current_pitch)*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 9: case 10: case 11:
					/* PC = 2 through 11, B cycle, write updated K1 through K10 */
					m_current_k[m_PC-2] += (((m_target_k[m_PC-2] - m_current_k[m_PC-2])*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 12: /* PC = 12, do nothing */
					break;
				}
			}
#endif
		}

		// calculate the output
		if (OLD_FRAME_UNVOICED_FLAG == 1)
		{
			// generate unvoiced samples here
			if (m_RNG & 1)
				m_excitation_data = ~0x3F; /* according to the patent it is (either + or -) half of the maximum value in the chirp table, so either 01000000(0x40) or 11000000(0xC0)*/
			else
				m_excitation_data = 0x40;
		}
		else /* (OLD_FRAME_UNVOICED_FLAG == 0) */
		{
			// generate voiced samples here
			/* US patent 4331836 Figure 14B shows, and logic would hold, that a pitch based chirp
			 * function has a chirp/peak and then a long chain of zeroes.
			 * The last entry of the chirp rom is at address 0b110011 (51d), the 52nd sample,
			 * and if the address reaches that point the ADDRESS incrementer is
			 * disabled, forcing all samples beyond 51d to be == 51d
			 */
			if (m_pitch_count >= 51)
				m_excitation_data = (INT8)m_coeff->chirptable[51];
			else /*m_pitch_count < 51*/
				m_excitation_data = (INT8)m_coeff->chirptable[m_pitch_count];
		}

		// Update LFSR *20* times every sample (once per T cycle), like patent shows
	for (i=0; i<20; i++)
	{
		bitout = ((m_RNG >> 12) & 1) ^
				((m_RNG >>  3) & 1) ^
				((m_RNG >>  2) & 1) ^
				((m_RNG >>  0) & 1);
		m_RNG <<= 1;
		m_RNG |= bitout;
	}
	/////////////////////////////

		this_sample = lattice_filter(); /* execute lattice filter */
		/* next, force result to 14 bits (since its possible that the addition at the final (k1) stage of the lattice overflowed) */
		while (this_sample > 16383) this_sample -= 32768;
		while (this_sample < -16384) this_sample += 32768;
		if (m_digital_select == 0) // analog SPK pin output is only 8 bits, with clipping
		  //			buffer[buf_count] = clip_analog(this_sample);
		  sample= clip_analog(this_sample);
		else // digital I/O pin output is 12 bits
		{
#ifdef ALLOW_4_LSB
			// input:  ssss ssss ssss ssss ssnn nnnn nnnn nnnn
			// N taps:                       ^                 = 0x2000;
			// output: ssss ssss ssss ssss snnn nnnn nnnn nnnN
		  //			buffer[buf_count] = (this_sample<<1)|((this_sample&0x2000)>>13);
		  sample=(this_sample<<1)|((this_sample&0x2000)>>13);
#else
			this_sample &= ~0xF;
			// input:  ssss ssss ssss ssss ssnn nnnn nnnn 0000
			// N taps:                       ^^ ^^^            = 0x3E00;
			// output: ssss ssss ssss ssss snnn nnnn nnnN NNNN
			//			buffer[buf_count] = (this_sample<<1)|((this_sample&0x3E00)>>9);
			sample=(this_sample<<1)|((this_sample&0x3E00)>>9);
#endif
		}
		// Update all counts

		m_subcycle++;
		if ((m_subcycle == 2) && (m_PC == 12))
		{
			/* Circuit 412 in the patent acts a reset, resetting the pitch counter to 0
			 * if INHIBIT was true during the most recent frame transition.
			 * The exact time this occurs is betwen IP=7, PC=12 sub=0, T=t12
			 * and m_IP = 0, PC=0 sub=0, T=t12, a period of exactly 20 cycles,
			 * which overlaps the time OLDE and OLDP are updated at IP=7 PC=12 T17
			 * (and hence INHIBIT itself 2 t-cycles later). We do it here because it is
			 * convenient and should make no difference in output.
			 */
			if ((m_IP == 7)&&(m_inhibit==1)) m_pitch_count = 0;
			m_subcycle = m_subc_reload;
			m_PC = 0;
			m_IP++;
			m_IP&=0x7;
		}
		else if (m_subcycle == 3)
		{
			m_subcycle = m_subc_reload;
			m_PC++;
		}
		m_pitch_count++;
		if (m_pitch_count >= m_current_pitch) m_pitch_count = 0; 
		m_pitch_count &= 0x1FF;
	return sample;
}

int16_t process5100raw()
{
  static u8 counter;
	int i, bitout, zpar;
	INT32 this_sample;
	int16_t sample, val;
	/* loop until the buffer is full or we've stopped speaking */
	//	if (m_speaking_now)
	//	{

		/* if we're ready for a new frame to be applied, i.e. when IP=0, PC=12, Sub=1
		 * (In reality, the frame was really loaded incrementally during the entire IP=0
		 * PC=x time period, but it doesn't affect anything until IP=0 PC=12 happens)
		 */
		if ((m_IP == 0) && (m_PC == 12) && (m_subcycle == 1))
		{
			// HACK for regression testing, be sure to comment out before release!
			//m_RNG = 0x1234;
			// end HACK

			/* appropriately override the interp count if needed; this will be incremented after the frame parse! */
			m_IP = reload_table[m_c_variant_rate&0x3];

#ifdef PERFECT_INTERPOLATION_HACK
			/* remember previous frame energy, pitch, and coefficients */
			m_old_frame_energy_idx = m_new_frame_energy_idx;
			m_old_frame_pitch_idx = m_new_frame_pitch_idx;
			for (i = 0; i < m_coeff->num_k; i++)
				m_old_frame_k_idx[i] = m_new_frame_k_idx[i];
#endif


			/* Parse a new frame into the new_target_energy, new_target_pitch and new_target_k[] */
			if (m_coeff->pitch_bits==5)			parse_frame_raw_5100();
			else parse_frame_raw_5200();

			/* in all cases where interpolation would be inhibited, set the inhibit flag; otherwise clear it.
			   Interpolation inhibit cases:
			 * Old frame was voiced, new is unvoiced
			 * Old frame was silence/zero energy, new has nonzero energy
			 * Old frame was unvoiced, new is voiced (note this is the case on the patent but may not be correct on the real final chip)
			 */

			/* load new frame targets from tables, using parsed indices */
			m_target_energy = m_coeff->energytable[m_new_frame_energy_idx];
			m_target_pitch = m_coeff->pitchtable[m_new_frame_pitch_idx];
			//zpar = NEW_FRAME_UNVOICED_FLAG; // find out if parameters k5-k10 should be zeroed
			for (i = 0; i < 4; i++)
				m_target_k[i] = m_coeff->ktable[i][m_new_frame_k_idx[i]];
			for (i = 4; i < m_coeff->num_k; i++)
				m_target_k[i] = (m_coeff->ktable[i][m_new_frame_k_idx[i]]);

		}
		else // Not a new frame, just interpolate the existing frame.
		{
		  //		  m_inhibit = 0; thus
		  //		  int inhibit_state = ((m_inhibit==1)&&(m_IP != 0)); // disable inhibit when reaching the last interp period, but don't overwrite the m_inhibit value
		  int inhibit_state = 0;
		  //#ifdef PERFECT_INTERPOLATION_HACK

#ifdef PERFECT_INTERPOLATION_HACK
			int samples_per_frame = m_subc_reload?175:266; // either (13 A cycles + 12 B cycles) * 7 interps for normal SPEAK/SPKEXT, or (13*2 A cycles + 12 B cycles) * 7 interps for SPKSLOW
			//int samples_per_frame = m_subc_reload?200:304; // either (13 A cycles + 12 B cycles) * 8 interps for normal SPEAK/SPKEXT, or (13*2 A cycles + 12 B cycles) * 8 interps for SPKSLOW
			// try changing this up/down
			
			int current_sample = (m_subcycle - m_subc_reload)+(m_PC*(3-m_subc_reload))+((m_subc_reload?25:38)*((m_IP-1)&7));

			zpar = OLD_FRAME_UNVOICED_FLAG;
			//fprintf(stderr, "CS: %03d", current_sample);
			// reset the current energy, pitch, etc to what it was at frame start
			m_current_energy = m_coeff->energytable[m_old_frame_energy_idx];
			m_current_pitch = m_coeff->pitchtable[m_old_frame_pitch_idx];
			for (i = 0; i < 4; i++)
				m_current_k[i] = m_coeff->ktable[i][m_old_frame_k_idx[i]];
			for (i = 4; i < m_coeff->num_k; i++)
				m_current_k[i] = (m_coeff->ktable[i][m_old_frame_k_idx[i]] * (1-zpar));
			// now adjust each value to be exactly correct for each of the samples per frame
			if (m_IP != 0) // if we're still interpolating...
			{
				m_current_energy += (((m_target_energy - m_current_energy)*(1-inhibit_state))*current_sample)/samples_per_frame;
				m_current_pitch += (((m_target_pitch - m_current_pitch)*(1-inhibit_state))*current_sample)/samples_per_frame;
				for (i = 0; i < m_coeff->num_k; i++)
					m_current_k[i] += (((m_target_k[i] - m_current_k[i])*(1-inhibit_state))*current_sample)/samples_per_frame;
			}
			else // we're done, play this frame for 1/8 frame.
			{
				m_current_energy = m_target_energy;
				m_current_pitch = m_target_pitch;
				for (i = 0; i < m_coeff->num_k; i++)
					m_current_k[i] = m_target_k[i];
			}
#else
			//Updates to parameters only happen on subcycle '2' (B cycle) of PCs.
			if (m_subcycle == 2)
			{
				switch(m_PC)
				{
					case 0: /* PC = 0, B cycle, write updated energy */
					m_current_energy += (((m_target_energy - m_current_energy)*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 1: /* PC = 1, B cycle, write updated pitch */
					m_current_pitch += (((m_target_pitch - m_current_pitch)*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 9: case 10: case 11:
					/* PC = 2 through 11, B cycle, write updated K1 through K10 */
					m_current_k[m_PC-2] += (((m_target_k[m_PC-2] - m_current_k[m_PC-2])*(1-inhibit_state)) INTERP_SHIFT);
					break;
					case 12: /* PC = 12, do nothing */
					break;
				}
			}
#endif
		}

		// calculate the output
			// generate voiced samples here
			/* US patent 4331836 Figure 14B shows, and logic would hold, that a pitch based chirp
			 * function has a chirp/peak and then a long chain of zeroes.
			 * The last entry of the chirp rom is at address 0b110011 (51d), the 52nd sample,
			 * and if the address reaches that point the ADDRESS incrementer is
			 * disabled, forcing all samples beyond 51d to be == 51d
			 */
			if (m_pitch_count >= 51)
				m_excitation_data = (INT8)m_coeff->chirptable[51];
			else /*m_pitch_count < 51*/
				m_excitation_data = (INT8)m_coeff->chirptable[m_pitch_count];
	/////////////////////////////

		this_sample = lattice_filter(); /* execute lattice filter */
		/* next, force result to 14 bits (since its possible that the addition at the final (k1) stage of the lattice overflowed) */
		while (this_sample > 16383) this_sample -= 32768;
		while (this_sample < -16384) this_sample += 32768;
		if (m_digital_select == 0) // analog SPK pin output is only 8 bits, with clipping
		  //			buffer[buf_count] = clip_analog(this_sample);
		  sample= clip_analog(this_sample);
		else // digital I/O pin output is 12 bits
		{
#ifdef ALLOW_4_LSB
			// input:  ssss ssss ssss ssss ssnn nnnn nnnn nnnn
			// N taps:                       ^                 = 0x2000;
			// output: ssss ssss ssss ssss snnn nnnn nnnn nnnN
		  //			buffer[buf_count] = (this_sample<<1)|((this_sample&0x2000)>>13);
		  sample=(this_sample<<1)|((this_sample&0x2000)>>13);
#else
			this_sample &= ~0xF;
			// input:  ssss ssss ssss ssss ssnn nnnn nnnn 0000
			// N taps:                       ^^ ^^^            = 0x3E00;
			// output: ssss ssss ssss ssss snnn nnnn nnnN NNNN
			//			buffer[buf_count] = (this_sample<<1)|((this_sample&0x3E00)>>9);
			sample=(this_sample<<1)|((this_sample&0x3E00)>>9);
#endif
		}
		// Update all counts

		m_subcycle++;
		if ((m_subcycle == 2) && (m_PC == 12))
		{
			/* Circuit 412 in the patent acts a reset, resetting the pitch counter to 0
			 * if INHIBIT was true during the most recent frame transition.
			 * The exact time this occurs is betwen IP=7, PC=12 sub=0, T=t12
			 * and m_IP = 0, PC=0 sub=0, T=t12, a period of exactly 20 cycles,
			 * which overlaps the time OLDE and OLDP are updated at IP=7 PC=12 T17
			 * (and hence INHIBIT itself 2 t-cycles later). We do it here because it is
			 * convenient and should make no difference in output.
			 */
			if ((m_IP == 7)&&(m_inhibit==1)) m_pitch_count = 0;
			m_subcycle = m_subc_reload;
			m_PC = 0;
			m_IP++;
			m_IP&=0x7;
		}
		else if (m_subcycle == 3)
		{
			m_subcycle = m_subc_reload;
			m_PC++;
		}
		m_pitch_count++;
		if (m_pitch_count >= m_current_pitch) m_pitch_count = 0; 

		m_pitch_count &= 0x1FF;
		//		buf_count++;
		//		size--;
		//	}
	return sample;
}


int16_t process5100raw_no_inter() // this is with no interpolation
{
  static u8 counter;
	int i, bitout, zpar;
	INT32 this_sample;
	int16_t sample, val;

			/* Parse a new frame into the new_target_energy, new_target_pitch and new_target_k[] */
			if (m_coeff->pitch_bits==5)			parse_frame_raw_5100();
			else parse_frame_raw_5200();

			m_current_energy = m_coeff->energytable[m_new_frame_energy_idx];
			m_current_pitch = m_coeff->pitchtable[m_new_frame_pitch_idx];
			zpar = NEW_FRAME_UNVOICED_FLAG; // find out if parameters k5-k10 should be zeroed
			for (i = 0; i < 4; i++)
				m_current_k[i] = m_coeff->ktable[i][m_new_frame_k_idx[i]];
			for (i = 4; i < m_coeff->num_k; i++)
				m_current_k[i] = (m_coeff->ktable[i][m_new_frame_k_idx[i]]);

			
			// generate voiced samples here
			/* US patent 4331836 Figure 14B shows, and logic would hold, that a pitch based chirp
			 * function has a chirp/peak and then a long chain of zeroes.
			 * The last entry of the chirp rom is at address 0b110011 (51d), the 52nd sample,
			 * and if the address reaches that point the ADDRESS incrementer is
			 * disabled, forcing all samples beyond 51d to be == 51d
			 */
			if (m_pitch_count >= 51)
				m_excitation_data = (INT8)m_coeff->chirptable[51];
			else /*m_pitch_count < 51*/
				m_excitation_data = (INT8)m_coeff->chirptable[m_pitch_count];


		this_sample = lattice_filter(); /* execute lattice filter */
		/* next, force result to 14 bits (since its possible that the addition at the final (k1) stage of the lattice overflowed) */
		while (this_sample > 16383) this_sample -= 32768;
		while (this_sample < -16384) this_sample += 32768;
		if (m_digital_select == 0) // analog SPK pin output is only 8 bits, with clipping
		  //			buffer[buf_count] = clip_analog(this_sample);
		  sample= clip_analog(this_sample);
		else // digital I/O pin output is 12 bits
		{
#ifdef ALLOW_4_LSB
			// input:  ssss ssss ssss ssss ssnn nnnn nnnn nnnn
			// N taps:                       ^                 = 0x2000;
			// output: ssss ssss ssss ssss snnn nnnn nnnn nnnN
		  //			buffer[buf_count] = (this_sample<<1)|((this_sample&0x2000)>>13);
		  sample=(this_sample<<1)|((this_sample&0x2000)>>13);
#else
			this_sample &= ~0xF;
			// input:  ssss ssss ssss ssss ssnn nnnn nnnn 0000
			// N taps:                       ^^ ^^^            = 0x3E00;
			// output: ssss ssss ssss ssss snnn nnnn nnnN NNNN
			//			buffer[buf_count] = (this_sample<<1)|((this_sample&0x3E00)>>9);
			sample=(this_sample<<1)|((this_sample&0x3E00)>>9);
#endif
		}
		// Update all counts

		m_subcycle++;
		m_pitch_count++;
		if (m_pitch_count >= m_current_pitch) m_pitch_count = 0; 

		m_pitch_count &= 0x1FF;
		//		buf_count++;
		//		size--;
		//	}
	return sample;
}



void tms_init()
{

  INT16 i,j;

  for(i=0;i<256;i++)
    {
      j = (i>>4) | (i<<4); // Swap in groups of 4
      j = ((j & 0xcc)>>2) | ((j & 0x33)<<2); // Swap in groups of 2
      byte_rev[i] = ((j & 0xaa)>>1) | ((j & 0x55)<<1); // Swap bit pairs
      }

  m_coeff=vocab_alphons.m_coeff; // init for first bank which is 5100
  
  // what needs to be set to start speaking?

  m_new_frame_energy_idx = 0;
  m_new_frame_pitch_idx = 0;
  for (i = 0; i < 4; i++)
    m_new_frame_k_idx[i] = 0;
  for (i = 4; i < 7; i++)
    m_new_frame_k_idx[i] = 0xF;
  for (i = 7; i < m_coeff->num_k; i++)    m_new_frame_k_idx[i] = 0x7;
  m_talk_status = 1;

  tms_newsay();
  tms_newsay_allphon();
  
};

//////////////////{{{{ NEWSAY

void tms_newsay(){
  m_digital_select=1;
  u8 whichbank;
  
  // selx is pitch bend, sely is bank, selz is phrase - trigger re_inits this
  whichbank=1;//_sely*TMS_VOCAB_F;
  //  MAXED(whichbank, TMS_VOCAB_TOP);
  //  whichbank=TMS_VOCAB_TOP-whichbank; // inversion
  //  if (whichbank<19 || whichbank==71) last5100=whichbank;
  //  else last5200=whichbank;
  //  whichbank=56;
  m_coeff=allTMSvocabs[whichbank]->m_coeff;
  
  m_new_frame_energy_idx = 0;
  m_new_frame_pitch_idx = 0;
  m_talk_status = 1;

  INT16 sel=_selz*allTMSvocabs[whichbank]->extentplus;
  MAXED(sel, allTMSvocabs[whichbank]->extent);
    sel=allTMSvocabs[whichbank]->extent-sel; // inversion
  ptrAddr=allTMSvocabs[whichbank]->wordlist[sel]; 
  ptrBit = 0;

  //  		if ((m_IP == 0) && (m_PC == 12) && (m_subcycle == 1))
  m_IP=0; m_PC=12;m_subcycle=1; // ???
  
};

void tms_newsay_lowbit(){
  m_digital_select=0;
  u8 whichbank;
  
  // selx is pitch bend, sely is bank, selz is phrase - trigger re_inits this
  whichbank=_sely*TMS_VOCAB_F; 
  MAXED(whichbank, TMS_VOCAB_TOP);
  whichbank=TMS_VOCAB_TOP-whichbank; // inversion
  if (whichbank<19 || whichbank==71) last5100=whichbank;
  else last5200=whichbank;

  m_coeff=allTMSvocabs[whichbank]->m_coeff;
  
  m_new_frame_energy_idx = 0;
  m_new_frame_pitch_idx = 0;
  m_talk_status = 1;

  INT16 sel=_selz*allTMSvocabs[whichbank]->extentplus; 
  MAXED(sel, allTMSvocabs[whichbank]->extent);
  sel=allTMSvocabs[whichbank]->extent-sel; // inversion
  ptrAddr=allTMSvocabs[whichbank]->wordlist[sel]; 
  ptrBit = 0;
    m_IP=0; m_PC=12;m_subcycle=1; // ???

};

void tms_newsay_allphon(){
  m_digital_select=1;
  u8 whichbank=ALLPHON_BANK; 
  last5200=whichbank;
  m_coeff=allTMSvocabs[whichbank]->m_coeff;
  m_new_frame_energy_idx = 0;
  m_new_frame_pitch_idx = 0;
  m_talk_status = 1;
  INT16 sel=(int)(_selz*allTMSvocabs[whichbank]->extentplus);
  MAXED(sel, allTMSvocabs[whichbank]->extent);
  //  sel=allTMSvocabs[whichbank]->extent-sel; // inversion
  ptrAddr=allTMSvocabs[whichbank]->wordlist[sel]; 
  ptrBit = 0;
  m_IP=0; m_PC=12;m_subcycle=1; // ???
};

void tms_newsay_TTS(){
  // this is our allphon code
  m_digital_select=1;
  m_coeff=allTMSvocabs[ALLPHON_BANK]->m_coeff;
  m_new_frame_energy_idx = 0;
  m_new_frame_pitch_idx = 0;
  m_talk_status = 1;

    ptrAddr=allTMSvocabs[ALLPHON_BANK]->wordlist[TTSoutarray[TTSindex]]; 
    ptrBit = 0;

    TTSindex++;
   if (TTSindex>=TTSlength) {
     TTSindex=0;
     TTSlength= text2speechforTMS(16,TTSinarray,TTSoutarray);
     }
   m_IP=0; m_PC=12;m_subcycle=1; // ???
};

void tms_retriggerTTS(){
  m_digital_select=1;
  m_coeff=allTMSvocabs[ALLPHON_BANK]->m_coeff;
  m_new_frame_energy_idx = 0;
  m_new_frame_pitch_idx = 0;
  m_talk_status = 1;

  TTSlength= text2speechforTMS(16,TTSinarray,TTSoutarray); 
  ptrAddr=allTMSvocabs[ALLPHON_BANK]->wordlist[TTSoutarray[0]]; 
  ptrBit = 0;
  TTSindex=1; // for next round
  m_IP=0; m_PC=12;m_subcycle=1; // ???
}

void tms_newsay_raw5100(){
  m_digital_select=1;
  u8 whichbank=0; // 5100 
  m_coeff=allTMSvocabs[whichbank]->m_coeff;
    m_new_frame_energy_idx = 0;
  m_new_frame_pitch_idx = 0;
  m_talk_status = 1;
  m_IP=0; m_PC=12;m_subcycle=1; // ??? but this always stays so...
};

void tms_newsay_raw5200(){
  m_digital_select=1;
  m_coeff=&T0285_2501E_coeff;
  m_new_frame_energy_idx = 0;
  m_new_frame_pitch_idx = 0;
  m_talk_status = 1;
  m_IP=0; m_PC=12;m_subcycle=1; // ???

};

void tms_newsay_raw5220(){
  m_digital_select=1;
  m_coeff=&tms5220_coeff;
  m_new_frame_energy_idx = 0;
  m_new_frame_pitch_idx = 0;
  m_talk_status = 1;
  //    m_IP=0; m_PC=12;m_subcycle=1; // ???

};

void tms_newsay_specific(u8 whichbank){
  m_digital_select=1;
  m_coeff=allTMSvocabs[whichbank]->m_coeff;
  
  m_new_frame_energy_idx = 0;
  m_new_frame_pitch_idx = 0;
  m_talk_status = 1;

  INT16 sel=_selz*allTMSvocabs[whichbank]->extentplus; 
  MAXED(sel, allTMSvocabs[whichbank]->extent);
  sel=allTMSvocabs[whichbank]->extent-sel; // inversion
  ptrAddr=allTMSvocabs[whichbank]->wordlist[sel]; 
  ptrBit = 0;
    m_IP=0; m_PC=12;m_subcycle=1; // ???

};

/// newsays for specific vocabs/chipsets

void tms_newsay_specifica(){
  tms_newsay_specific(ALLPHON_BANK);
}

void tms_newsay_specificx(){
  tms_newsay_specific(last5200);  //
}


void tms_newsay_specific5100(){ 
  tms_newsay_specific(last5100); // TODO: or these could randomly select from list or sel the last bank of 5100/5200
}

///// get_samples

int16_t tms_get_sample_allphon(){
  modus=10; // length and pitch
  int16_t sample; u8 ending=0;
  sample=  process(&ending);
  if (ending==1){
    tms_newsay_allphon();
  }
  return sample;
}

int16_t tms_get_sample_allphon_sing(){
  modus=10; // abs pitch and length
  int16_t sample; u8 ending=0;
  sample=  process(&ending);// 
  if (ending==1){
    tms_newsay_allphon();
  }
  return sample;
}


int16_t tms_get_sample(){
  int16_t sample; u8 ending=0;
  modus=1;
  sample=  process(&ending);
  if (ending==1){
    //    if (triggermode==0) tms_newsay(); // where we repeat --- if we last triggered then do not repeat? TESTY!
    tms_newsay();
  }
  return sample;
}

int16_t tms_get_sample_sing(){
  int16_t sample; u8 ending=0;
  modus=2;
  sample=  process(&ending);
  if (ending==1){
    tms_newsay();
  }
  return sample;
}

int16_t tms_get_sample_lowbit(){
  modus=1;
  int16_t sample; u8 ending=0;
  sample=  process(&ending);
  if (ending==1){
    tms_newsay_lowbit();
  }
  return sample;
}

int16_t tms_get_sample_TTS(){
  modus=1;
  int16_t sample; u8 ending=0;
  sample=  process(&ending);
  if (ending==1){
    tms_newsay_TTS();
  }
  return sample;
}

int16_t tms_get_sample_bendlength(){
  modus=4; // no pitch
  int16_t sample; u8 ending=0;
  sample=  process(&ending);
  if (ending==1){
    tms_newsay();
  }
  return sample;
}

int16_t tms_get_sample_raw5100(){
  modus=0;
  m_coeff=&T0280B_0281A_coeff;
  int16_t sample;
  sample=  process5100raw();
  return sample;
}

int16_t tms_get_sample_raw5200(){
  modus=0;
  m_coeff=&T0285_2501E_coeff;
  int16_t sample; 
  sample=  process5100raw(); // combined...
  return sample;
}

int16_t tms_get_sample_raw5220(){
  modus=0;
  m_coeff=&tms5220_coeff;
  int16_t sample; 
  sample=  process5100raw(); // combined
  return sample;
}

//// bends and vocab selection

int16_t tms_get_sample_bend5100(u8 fix){
  modus=0;
  m_coeff=&T0280B_0281A_coeff;
  int16_t sample; u8 ending=0;
  sample=processbend5100(&ending);
  if (ending==1){
    tms_newsay_specific(fix);
  }
  return sample;
}

int16_t tms_get_sample_bend5100w(){// for a 5100 vocab such as 0 to begin with
  tms_get_sample_bend5100(last5100);
}

int16_t tms_get_sample_bend5200(u8 fix){ 
  modus=0;
  m_coeff=&T0285_2501E_coeff;
  int16_t sample; u8 ending=0;
  sample=processbend5100(&ending); // combined
  if (ending==1){
    tms_newsay_specific(fix); 
  }
  return sample;
}

int16_t tms_get_sample_bend5200a(){ // for allphons 
  tms_get_sample_bend5200(ALLPHON_BANK);
}

int16_t tms_get_sample_bend5200x(){ 
  tms_get_sample_bend5200(last5200);
}


int16_t tms_get_sample_5100pitchtable(u8 fix){ 
  modus=0;
  m_coeff=&T0280B_0281A_coeff;
  int16_t sample; u8 ending=0;
  sample= process_pitch_tabled5100(&ending);
  if (ending==1){
    tms_newsay_specific(fix); 
  }
  return sample;
}

int16_t tms_get_sample_5100pitchtablew(){  
  tms_get_sample_5100pitchtable(last5100);
}

int16_t tms_get_sample_5200pitchtable(u8 fix){ 
  modus=0;
  m_coeff=&T0285_2501E_coeff;
  int16_t sample; u8 ending=0;
  sample= process_pitch_tabled5100(&ending);
  if (ending==1){
    tms_newsay_specific(fix);
  }
  return sample;
}

int16_t tms_get_sample_5200pitchtablea(){
  tms_get_sample_5200pitchtable(ALLPHON_BANK);
}

int16_t tms_get_sample_5200pitchtablex(){
  tms_get_sample_5200pitchtable(last5200);
}


int16_t tms_get_sample_5100ktable(u8 fix){
  modus=0;
  m_coeff=&T0280B_0281A_coeff;
  int16_t sample; u8 ending=0;
  sample= process_k_tabled5100(&ending);
  if (ending==1){
    tms_newsay_specific(fix); 
  }
  return sample;
}

int16_t tms_get_sample_5100ktablew(){  /// additional vocabs as extra functions - which ones?
  tms_get_sample_5100ktable(last5100);
}

int16_t tms_get_sample_5200ktable(u8 fix){
  modus=0;
  m_coeff=&T0285_2501E_coeff; //5200
  int16_t sample; u8 ending=0;
  sample= process_k_tabled5100(&ending);
  if (ending==1){
    tms_newsay_specific(fix);
  }
  return sample;
}

int16_t tms_get_sample_5200ktablea(){ 
tms_get_sample_5200ktable(ALLPHON_BANK);
}

int16_t tms_get_sample_5100kandpitchtable(u8 fix){
  modus=0;
  m_coeff=&T0280B_0281A_coeff; // 5100
  int16_t sample; u8 ending=0;
  sample= process_pitchk_tabled5100(&ending);
  if (ending==1){
    tms_newsay_specific(fix); 
  }
  return sample;
}

/*
int16_t tms_get_sample5100fulltablex(u8 fix){
  modus=0;
  m_coeff=&T0280B_0281A_coeff; // 5100
  int16_t sample; u8 ending=0;
  sample= process_full_5100(&ending);
  if (ending==1){
    tms_newsay_specific(fix); 
  }
  return sample;
}
*/

int16_t tms_get_sample_5100kandpitchtablew(){ 
  tms_get_sample_5100kandpitchtable(last5100);
}

int16_t tms_get_sample_5200kandpitchtable(u8 fix){
  modus=0;
  m_coeff=&T0285_2501E_coeff; // 5200
  int16_t sample; u8 ending=0;
  sample= process_pitchk_tabled5100(&ending); // now we recognise pitchbits 
  if (ending==1){
    tms_newsay_specific(fix); 
  }
  return sample;
}

int16_t tms_get_sample_5200kandpitchtablea(){ 
  tms_get_sample_5200kandpitchtable(ALLPHON_BANK);
}

/*
int16_t tms_get_sample5100fulltable(){
  tms_get_sample5100fulltablex(1);
}
*/

int16_t tms_get_sample_5200kandpitchtablex(){ 
  tms_get_sample_5200kandpitchtable(last5200);
}
