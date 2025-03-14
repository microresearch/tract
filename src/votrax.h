// license:BSD-3-Clause
// copyright-holders:Olivier Galibert
/***************************************************************************

    votrax.h

    Votrax SC01A simulation

***************************************************************************/

// PROBS:

// bitswap - is only one TODO

#define BIT(x,n) (((x)>>(n))&1)

/*

BIT:

template <typename T, typename U, typename... V> constexpr T bitswap(T val, U b, V... c)
{
	return (BIT(val, b) << sizeof...(c)) | bitswap(val, c...);
}

m_rom_f1  = bitswap(val,  0,  7, 14, 21);
			
BIT(val, 7) << 1 | BIT(val, 14) << 

#define BITSWAP8(val,B7,B6,B5,B4,B3,B2,B1,B0) \
	((BIT(val,B7) << 7) | (BIT(val,B6) << 6) | (BIT(val,B5) << 5) | (BIT(val,B4) << 4) | \
		(BIT(val,B3) << 3) | (BIT(val,B2) << 2) | (BIT(val,B1) << 1) | (BIT(val,B0) << 0))


eg.			m_rom_fc  = bitswap(val,  3, 10, 17, 24);

test out as:
	((BIT(val,3) << 3) | (BIT(val,10) << 2) | (BIT(val,17) << 1) | (BIT(val,24) << 0)

we have either 4 bits or 1 or 7 so...

*/

#define BITSWAP4(val,B3,B2,B1,B0) \
  ((BIT(val,B3) << 3) | (BIT(val,B2) << 2) | (BIT(val,B1) << 1) | (BIT(val,B0) << 0))

#define BITSWAP7(val,B6,B5,B4,B3,B2,B1,B0) \
	((BIT(val,B6) << 6) | (BIT(val,B5) << 5) | (BIT(val,B4) << 4) | \
		(BIT(val,B3) << 3) | (BIT(val,B2) << 2) | (BIT(val,B1) << 1) | (BIT(val,B0) << 0))

#define BITSWAP1(val,B0) \
  ((BIT(val,B0) << 0))


#ifndef VOTRAX_H
#define VOTRAX_H

#ifdef LAP
#define u8 unsigned char
#define u16 unsigned int
#define u32 unsigned int
#endif
#define u32 unsigned int
#define u64 unsigned long long 
#define bool unsigned char
#define true 1
#define false 0

//	sound_stream *m_stream;                         // Output stream
//	emu_timer *m_timer;                             // General timer
//	required_memory_region m_rom;                   // Internal ROM
	myfloat m_mainclock;                                // Current main clock
	myfloat m_sclock;                                // Stream sample clock (40KHz, main/18)
	myfloat m_cclock;                                // 20KHz capacitor switching clock (main/36)
	u32 m_sample_count;                             // Sample counter, to cadence chip updates

	// Inputs
	u8 m_inflection;                                // 2-bit inflection value
	u8 m_phone;                                     // 6-bit phone value

	// Outputs
//	devcb_write_line m_ar_cb;                       // Callback for ar
	bool m_ar_state;                                // Current ar state

	// "Unpacked" current rom values
u8 m_rom_duration, m_rom_duration_orig;                              // Duration in 5KHz units (main/144) of one tick, 16 ticks per phone, 7 bits
u8 m_rom_vd, m_rom_cld;                         // Duration in ticks of the "voice" and "closure" delays, 4 bits
u8 m_rom_fa, m_rom_fc, m_rom_va;                // Analog parameters, noise volume, noise freq cutoff and voice volume, 4 bits each
u8 m_rom_f1, m_rom_f2, m_rom_f2q, m_rom_f3;     // Analog parameters, formant frequencies and Q, 4 bits each
u8 m_rom_vd_orig, m_rom_cld_orig;                         // Duration in ticks of the "voice" and "closure" delays, 4 bits
u8 m_rom_fa_orig, m_rom_fc_orig, m_rom_va_orig;                // Analog parameters, noise volume, noise freq cutoff and voice volume, 4 bits each
u8 m_rom_f1_orig, m_rom_f2_orig, m_rom_f2q_orig, m_rom_f3_orig;     // Analog parameters, formant frequencies and Q, 4 bits each



bool m_rom_closure;                             // Closure bit, true = silence at cld
	bool m_rom_pause;                               // Pause bit

	// Current interpolated values (8 bits each)
u8 m_cur_fa, m_cur_fc, m_cur_va;
u8 m_cur_fa_orig, m_cur_fc_orig, m_cur_va_orig;
u8 m_cur_f1, m_cur_f2, m_cur_f2q, m_cur_f3;
u8 m_cur_f1_orig, m_cur_f2_orig, m_cur_f2q_orig, m_cur_f3_orig;
	// Current committed values
	u8 m_filt_fa, m_filt_fc, m_filt_va;             // Analog parameters, noise volume, noise freq cutoff and voice volume, 4 bits each
	u8 m_filt_f1, m_filt_f2, m_filt_f2q, m_filt_f3; // Analog parameters, formant frequencies/Q on 4 bits except f2 on 5 bits

	// Internal counters
	u16 m_phonetick;                                // 9-bits phone tick duration counter
	u8  m_ticks;                                    // 5-bits tick counter
	u8  m_pitch;                                    // 7-bits pitch counter
	u8  m_closure;                                  // 5-bits glottal closure counter
	u8  m_update_counter;                           // 6-bits counter for the 625Hz (main/1152) and 208Hz (main/3456) update timing generators

	// Internal state
	bool m_cur_closure;                             // Current internal closure state
	u16 m_noise;                                    // 15-bit noise shift register
	bool m_cur_noise;                               // Current noise output

	// Filter coefficients and level histories
	myfloat m_voice_1[4];
	myfloat m_voice_2[4];
	myfloat m_voice_3[4];

	myfloat m_noise_1[3];
	myfloat m_noise_2[3];
	myfloat m_noise_3[2];
	myfloat m_noise_4[2];

	myfloat m_vn_1[4];
	myfloat m_vn_2[4];
	myfloat m_vn_3[4];
	myfloat m_vn_4[4];
	myfloat m_vn_5[2];
	myfloat m_vn_6[2];

	myfloat m_f1_a[4],  m_f1_b[4];                   // F1 filtering
	myfloat m_f2v_a[4], m_f2v_b[4];                  // F2 voice filtering
	myfloat m_f2n_a[2], m_f2n_b[2];                  // F2 noise filtering
	myfloat m_f3_a[4],  m_f3_b[4];                   // F3 filtering
	myfloat m_f4_a[4],  m_f4_b[4];                   // F4 filtering
	myfloat m_fx_a[1],  m_fx_b[2];                   // Final filtering
	myfloat m_fn_a[3],  m_fn_b[3];                   // Noise shaping

	// Compute a total capacitor value based on which bits are currently active
static myfloat bits_to_caps(u32 value, u32 *caps_values, u8 howm) {
		myfloat total = 0.0f;
		myfloat d;
		u8 i;
		/*		for(myfloat d : caps_values) { // what is this doing?
			if(value & 1)
				total += d;
			value >>= 1;
			}*/
		for (i=0;i<howm;i++){
		  d=caps_values[i];
		  if(value & 1)
		    total += d;
		  value >>= 1;
			}
		return total;
	}

	// Shift a history of values by one and insert the new value at the front
static void shift_hist(myfloat val, myfloat *hist_array, u8 N) {
  for(u8 i=N-1; i>0; i--)
    hist_array[i] = hist_array[i-1];
  hist_array[0] = val;
	}

	// Shift a history of values by one and insert the new value at the front
/*	template<u32 N> static void shift_hist(double val, double (&hist_array)[N]) {
		for(u32 i=N-1; i>0; i--)
			hist_array[i] = hist_array[i-1];
		hist_array[0] = val;
		}*/



	// Apply a filter and compute the result. 'a' is applied to x (inputs) and 'b' to y (outputs)
static myfloat apply_filter(const myfloat *x, const myfloat *y, const myfloat *a, const myfloat *b, u8 Na, u8 Nb) {
		myfloat total = 0.0f;
		for(u8 i=0; i<Na; i++)
			total += x[i] * a[i];
		for(u8 i=1; i<Nb; i++)
			total -= y[i-1] * b[i];
		return total / b[0];
	}

/*older

	// Apply a filter and compute the result. 'a' is applied to x (inputs) and 'b' to y (outputs)
	template<u32 Nx, u32 Ny, u32 Na, u32 Nb> static double apply_filter(const double (&x)[Nx], const double (&y)[Ny], const double (&a)[Na], const double (&b)[Nb]) {
		double total = 0;
		for(u32 i=0; i<Na; i++)
			total += x[i] * a[i];
		for(u32 i=1; i<Nb; i++)
			total -= y[i-1] * b[i];
		return total / b[0];
	}


myfloat apply_filter(const myfloat *x, const myfloat *y, const myfloat *a, const myfloat *b)
{
  	return (x[0]*a[0] + x[1]*a[1] + x[2]*a[2] + x[3]*a[3] - y[0]*b[1] - y[1]*b[2] - y[2]*b[3]) / b[0];
}
*/

	void build_standard_filter(myfloat *a, myfloat *b,
							   myfloat c1t, // Unswitched cap, input, top
							   myfloat c1b, // Switched cap, input, bottom
							   myfloat c2t, // Unswitched cap, over first amp-op, top
							   myfloat c2b, // Switched cap, over first amp-op, bottom
							   myfloat c3,  // Cap between the two op-amps
							   myfloat c4); // Cap over second op-amp

	void build_noise_shaper_filter(myfloat *a, myfloat *b,
								   myfloat c1,  // Cap over first amp-op
								   myfloat c2t, // Unswitched cap between amp-ops, input, top
								   myfloat c2b, // Switched cap between amp-ops, input, bottom
								   myfloat c3,  // Cap over second amp-op
								   myfloat c4); // Switched cap after second amp-op

	void build_lowpass_filter(myfloat *a, myfloat *b,
							  myfloat c1t,  // Unswitched cap, over amp-op, top
							  myfloat c1b); // Switched cap, over amp-op, bottom

	void build_injection_filter(myfloat *a, myfloat *b,
								myfloat c1b, // Switched cap, input, bottom
								myfloat c2t, // Unswitched cap, over first amp-op, top
								myfloat c2b, // Switched cap, over first amp-op, bottom
								myfloat c3,  // Cap between the two op-amps
								myfloat c4); // Cap over second op-amp

	static u8 interpolate(u8 reg, u8 target);    // Do one interpolation step
	void chip_update();                             // Global update called at 20KHz (main/36)
	void filters_commit(bool force);                // Commit the currently computed interpolation values to the filters
	void phone_commit();                            // Commit the current phone id
	u32 analog_calc();                  // Compute one more sample


#endif /* VOTRAX_H */
