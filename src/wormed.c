//xxxxx/wormed voice touch synth 2017 

// Modes: 0/1/TMS/speak and spell vocal, 2/square wave, 3->11
// wavetables, 12-pulse with length of pulse and frequency, 13-sp0256
// pulse, 14-votrax glottal impulse, 15-impulses from klatt

// license:GPL-2.0+
// copyright-holders: Martin Howse

// NOTES:  ADC is 0 with no finger, do we want noise volume on adcread10(1)

//                  ^TOP^
//  fingers: X- pitch   X-unused
//
//                 X-mode

#include <stdint.h>

#define FS 8000 // sample rate
#define CHIRP_SIZE 41
#define howmany 10

//sample, rate, counter;
uint16_t synthEnergy=14;
volatile uint16_t location=0;

static uint16_t pitch[howmany];
static float floatrate;
uint8_t pitchindex=0;
uint16_t total=0, average=0;
uint8_t pulsecounter, pulselength;

static uint8_t periodCounter;
static int16_t x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,nsq;
static uint16_t synthRand = 1;
static float floatlocation;
int16_t u0,u1,u2,u3,u4,u5,u6,u7,u8,u9,u10;

extern float _selz, _selx;

uint8_t chirp[CHIRP_SIZE] = {0x00,0x2a,0xd4,0x32,0xb2,0x12,0x25,0x14,0x02,0xe1,0xc5,0x02,0x5f,0x5a,0x05,0x0f,0x26,0xfc,0xa5,0xa5,0xd6,0xdd,0xdc,0xfc,0x25,0x2b,0x22,0x21,0x0f,0xff,0xf8,0xee,0xed,0xef,0xf7,0xf6,0xfa,0x00,0x03,0x02,0x01};// tms

uint8_t glottal_wave[9] = {128, 54, 0, 237, 219, 201, 182, 164, 146};// klatt
uint8_t impulsive_source[3] = {128,0,250};
uint8_t rates[16] = {6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
uint8_t inks[16] =  {1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5};

void wormed_newsay(void){
}

int16_t chirp_get_sample(){
  static uint8_t periodCounter;
  uint8_t synthPeriod=255-(_selz*255);
  int16_t u10;
    if (periodCounter < synthPeriod) {
      periodCounter++;
    } else {
      periodCounter = 0;
    }
    if (periodCounter < CHIRP_SIZE) {
      u10 = ((chirp[periodCounter]) * (uint32_t) synthEnergy) >> 4; // try to increase volume >>
    } else {
      u10 = 0;
    }
    if (u10 > 511) u10 = 511;
    if (u10 < -512) u10 = -512;
    uint8_t nextPwm = (u10>>2)+0x80;
    return (nextPwm<<8);
}

//////

int16_t votglot_get_sample(){
  static uint8_t periodCounter;
  uint8_t synthPeriod=255-(_selz*255);
  uint8_t nextPwm = periodCounter >= (9 << 2) ? 0 : glottal_wave[periodCounter >> 2];
  
  // how we change pitch:
  periodCounter = (periodCounter + 1) & 0x7f; // inc
  
  if(periodCounter == (0x7f ^ synthPeriod)) periodCounter = 0; // 7 bits 128
  return (nextPwm<<8);
}

int16_t impulsive_get_sample(){
  static uint8_t periodCounter;
  uint8_t nextPwm;
  uint8_t synthPeriod=255-(_selz*255);
	if (periodCounter < 3)
	{
		nextPwm = impulsive_source[periodCounter];
	}
	else
	{
		nextPwm = 128;
	}
	if (periodCounter++>synthPeriod) periodCounter=0;
  return (nextPwm<<8);
}


static const int16_t natural_samples[100]= 
  {
    -310,-400,530,356,224,89,23,-10,-58,-16,461,599,536,701,770,
    605,497,461,560,404,110,224,131,104,-97,155,278,-154,-1165,
    -598,737,125,-592,41,11,-247,-10,65,92,80,-304,71,167,-1,122,
    233,161,-43,278,479,485,407,266,650,134,80,236,68,260,269,179,
    53,140,275,293,296,104,257,152,311,182,263,245,125,314,140,44,
    203,230,-235,-286,23,107,92,-91,38,464,443,176,98,-784,-2449,
    -1891,-1045,-1600,-1462,-1384,-1261,-949,-730
};

int16_t sampled_source_getsample()
{
  static int nper;
  int itemp;
  float ftemp;
  float result;
  float diff_value;
  int current_value;
  int next_value;
  float temp_diff;
  float vol;

    ftemp = (float) nper;
    ftemp = ftemp / 4.0;
    ftemp = ftemp * 100.0 ;
    itemp = (int) ftemp;

    temp_diff = ftemp - (float) itemp;
  
    current_value = natural_samples[itemp];
    next_value = natural_samples[itemp+1];

    diff_value = (float) next_value - (float) current_value;
    diff_value = diff_value * temp_diff;

    result = natural_samples[itemp] + diff_value;
    result = result * 2.0f;
    nper+=1;//(1.0-_selz);
    if (nper>4.0) nper=0.0;
    vol=10240+(1024*_selz);
    return (int)(vol*result);
}

  static const short B0[224] = 
  {
    1200,1142,1088,1038, 991, 948, 907, 869, 833, 799, 768, 738, 710, 683, 658,
    634, 612, 590, 570, 551, 533, 515, 499, 483, 468, 454, 440, 427, 415, 403,
    391, 380, 370, 360, 350, 341, 332, 323, 315, 307, 300, 292, 285, 278, 272,
    265, 259, 253, 247, 242, 237, 231, 226, 221, 217, 212, 208, 204, 199, 195,
    192, 188, 184, 180, 177, 174, 170, 167, 164, 161, 158, 155, 153, 150, 147,
    145, 142, 140, 137, 135, 133, 131, 128, 126, 124, 122, 120, 119, 117, 115, 
    113,111, 110, 108, 106, 105, 103, 102, 100, 99, 97, 96, 95, 93, 92, 91, 90,
    88, 87, 86, 85, 84, 83, 82, 80, 79, 78, 77, 76, 75, 75, 74, 73, 72, 71,
    70, 69, 68, 68, 67, 66, 65, 64, 64, 63, 62, 61, 61, 60, 59, 59, 58, 57, 
    57, 56, 56, 55, 55, 54, 54, 53, 53, 52, 52, 51, 51, 50, 50, 49, 49, 48, 48,
    47, 47, 46, 46, 45, 45, 44, 44, 43, 43, 42, 42, 41, 41, 41, 41, 40, 40,
    39, 39, 38, 38, 38, 38, 37, 37, 36, 36, 36, 36, 35, 35, 35, 35, 34, 34,33,
    33, 33, 33, 32, 32, 32, 32, 31, 31, 31, 31, 30, 30, 30, 30, 29, 29, 29, 29,
    28, 28, 28, 28, 27, 27
  };

int16_t natural_source_getsample()
{
  static int16_t nper;
  float lgtemp;
  static float vwave;
  uint16_t nopen=40+(int)(220.0*_selz);
  float pulse_shape_a;  /* Makes waveshape of glottal pulse when open   */
  float pulse_shape_b;  /* Makes waveshape of glottal pulse when open   */

    pulse_shape_b = B0[nopen-40];
    pulse_shape_a = (pulse_shape_b * nopen) * 0.333f;

  pulse_shape_a = (pulse_shape_b * nopen) * 0.333f;


  nper++;
  if (nper>250) nper=0;
  if (nper < nopen) 
  {
    pulse_shape_a -= pulse_shape_b;
    vwave += pulse_shape_a;
    lgtemp=vwave * 0.028f;

    return(lgtemp*1024.0);
  }
  else 
  {
    vwave = 0.0f;
    return 0;
  }
}

typedef struct
 {
  char *name;
  float a;
  float b;
  float c;
  float p1;
  float p2;
 }
resonator_t, *resonator_ptr;

static resonator_t rgl =
{"crit-damped glot low-pass filter"};

#define PI               3.1415927f


 static void setabc(long f, long bw, resonator_ptr rp)
{
  float minus_pi_t = -PI / 32000.0;
  float two_pi_t = -2.0f * minus_pi_t;
  float arg = minus_pi_t * bw;
  float r = expf(arg);              /* Let r  =  exp(-pi bw t) */
  rp->c = -(r * r);                /* Let c  =  -r**2 */
  arg = two_pi_t * f;
  rp->b = r * cosf(arg) * 2.0f;      /* Let b = r * 2*cos(2 pi f t) */
  rp->a = 1.0f - rp->b - rp->c;     /* Let a = 1.0 - b - c */
}

/* Generic resonator function */
static float resonator(resonator_ptr r, float input)
{
	register float x = r->a * input + r->b * r->p1 + r->c * r->p2;
	r->p2 = r->p1;
	r->p1 = x;
	return x;
}


int16_t impulsivex_get_sample()
{
  long temp;
  static int16_t nper;
  static float doublet[] = {0.0,13000000.0f,-13000000.0f};
  static float vwave, val;
  uint16_t len;
  
  len=60-(60*_selz);
  temp = 32000/len;
  setabc(0L, temp, &rgl);
  nper++;
  if (nper>len) nper=0;
  
  if (nper < 3) 
  {
    vwave = doublet[nper];
  }
  else 
  {
    vwave = 0.0f;
  }
  
    val=resonator(&(rgl),vwave);
    return val;//*1024.0;
  //    return vwave;
}
