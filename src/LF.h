typedef struct LF_ {
  double t;
  double period; //= fxs->_period;
  int dataLength; //= fxs->_dataLength;
  double te; //= fxs->_te;
  double Eo; //= fxs->_Eo;
  double alpham; //= fxs->_alpha;
  double wg; //= fxs->_wg;
  double Ee; // = fxs->_Ee;
  double epsilon; //= fxs->_epsilon;
  double ta,tppp; //= fxs->_ta;
  double tc; //= fxs->_tc;
  int k; //= fxs->_k;
}LF;

void initLF(LF* lf);

int16_t rosenberg_get_sample();
int16_t LF_get_sample();

typedef struct PAR {
  float dur,		/* duration of the created file */
	jitter,
	cq,
	K,		/* speed of closure */
	Fg,		/* glottal formant; Fg = 1/(2.T2) */
	F0,		/* Fg > F0 */
	DC,		/* DC flow (% of max amplitude) */
	noise;		/* pow(10, arg.noise/10) */

  long fs;		/* sampling rate */
   int amp,             /* maximum amplitude */
       NoiseDistWidth;	/* noise = uniform(0,...,NoiseDistWidth) */
  float Kvar,Shimmer;           /* speed closure variation, Shimmer */
} PAR;


typedef struct SHIM_ {
int i, j, k,			/* general use */
      P,			/* nominal Pitch period */
      T,			/* real Pitch period (jittered) */
      T2,			/* = 1/(w.Fg) */
      T3,			/* point in the closing phase where flow = 0 */
      T4,			/* instant where the rising flow = DCflow */
  w[500],			/* vector for additive white noise */
  x[500];
float aux,	x_pow,	w_pow,  J,	DeltaPer[2];
float S, DeltaShimmer[2];
} SHIM;

void shimmer_init(PAR* param, SHIM* shim);
int16_t shimmer_get_sample();


typedef struct {
  int samplingRate;// # Hertz
  int freq,n;// # F0, given in Hertz (T0 = 1 / F0)
  float CQ;// # closed quotient
  float amp;//, # amplitude of voicing
  float a, b, period, T0, OQ, openPhaseCorrection, amplitudeCorrectionFactor,totalOffset;// what else?
} klglott;

void klglott88_init(klglott* self, int samplingRate, int freq, float CQ, float amp);
void setfreq(klglott* self, int freq);
int16_t klg_get_sample();

int16_t LF2_get_sample();
