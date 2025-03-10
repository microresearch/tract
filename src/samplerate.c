// license:GPL-2.0+
// copyright-holders: Martin Howse
#include "stm32h7xx_hal.h"
#include "main.h"
#include <string.h>
#include <stdlib.h>
#include "wavetable.h"
#include "resources.h"
#include "tms5200x.h"
#include "process.h"
#include "audio.h"

#define LIMIN 32000 // or 32000?

#define WIDTH 16              
#define DELAY_SIZE 6 // was 40 --- 3*width=16 = 3*16=48-5=43 - use 7 for simplea

extern __IO uint16_t adc_buffer[10];
extern float _selx, _sely, _selz;
extern float smoothed_adc_value[5];
extern u8 _intmode;

typedef struct {
  float FIRData[49];
  float const *FIRCoef;
    int FIRPtr, numberTaps;
} TRMFIRFilter;


static TRMFIRFilter firfilt;

static const float filtertaps[49]={ 0.00000001f,  0.00000007f,  0.00000021f, -0.00000007f, -0.00000243f, -0.00000627f,  0.00000224f,  0.00004326f,  0.00006970f, -0.00008697f, -0.00043244f, -0.00028644f,  0.00113161f,  0.00232968f, -0.00061607f, -0.00709116f, -0.00592280f,  0.01120457f,  0.02472395f, -0.00136324f, -0.05604421f, -0.05124764f,  0.08785309f,  0.29650419f,  0.39847428f,  0.29650419f,  0.08785309f, -0.05124764f, -0.05604421f, -0.00136324f,  0.02472395f,  0.01120457f, -0.00592280f, -0.00709116f, -0.00061607f,  0.00232968f,  0.00113161f, -0.00028644f, -0.00043244f, -0.00008697f,  0.00006970f,  0.00004326f,  0.00000224f, -0.00000627f, -0.00000243f, -0.00000007f,  0.00000021f,  0.00000007f,  0.00000001f}; // filter coeffs for above settings

void samplerate_init(){
  firfilt.FIRCoef=filtertaps;
  firfilt.numberTaps=49; firfilt.FIRPtr=0;
}

inline int iincrement(int pointer, int modulus)
{
    if (++pointer >= modulus)
	return 0;

    return pointer;
}

inline int ddecrement(int pointer, int modulus)
{
    if (--pointer < 0)
return modulus - 1;

    return pointer;
}


static float doFIRFilter(TRMFIRFilter *filter, float input, u8 needOutput)
{
        if (needOutput) {
	float output = 0.0f;

	filter->FIRData[filter->FIRPtr] = input;

	for (u8 i = 0; i < filter->numberTaps; i++) {
	  output += filter->FIRData[filter->FIRPtr] * filter->FIRCoef[i];
	  //output=input;
	  filter->FIRPtr = iincrement(filter->FIRPtr, filter->numberTaps);
	}

	filter->FIRPtr = ddecrement(filter->FIRPtr, filter->numberTaps);
	return output;
    } else {
	  filter->FIRData[filter->FIRPtr] = input;
	  filter->FIRPtr = ddecrement(filter->FIRPtr, filter->numberTaps);
	return 0.0f;
	}
}


static int16_t delay_buffer[2] = { 0 }; // was 48 but it doesn't need to be so big

void new_data(int16_t data)
{
  //    for (u8 ii=0;ii<DELAY_SIZE-5;ii++)	
  delay_buffer[0] = delay_buffer[1];
    delay_buffer[1] = data;
}


void samplerate_simple(int16_t* out, float factor, u8 size, int16_t(*getsample)(void), void(*newsay)(void), float sampleratio, u8 triggered){
  float alpha;
  int16_t outt;
  static float time_now=0.0f;
  long last_time;
  static long int_time=0;
  //  static u8 triggered=0;
  factor*=sampleratio;
  //  factor=0.1;
    if (triggered) newsay();
  
    for (u8 ii=0;ii<size;ii++){
      
      if (time_now>32768){
	int_time=0; 
	time_now-=32768.0f;
      }

  //  out[ii]=getsample();
    alpha = time_now - (float)int_time;
    //    if (factor>1.0){
    //       out[ii]=InterpolateHermite4pt3oX(delay_buffer[DELAY_SIZE-8], delay_buffer[DELAY_SIZE-7], delay_buffer[DELAY_SIZE-6], delay_buffer[DELAY_SIZE-5], alpha);
      //      out[ii]=0;
    //    }

    //    factor=0.5f;
    /*    if (factor>1){
    // low pass for decimation - 32000 / factor = max 4 = 8khz = 16->4kHz
      // int to float
      float interpolatedValue=(float32_t)(((float)delay_buffer[DELAY_SIZE-5] * alpha) + ((float)delay_buffer[DELAY_SIZE-6] * (1.0f - alpha)))/32768.0f;
      float sample = doFIRFilter(&firfilt, interpolatedValue, 1); //sample is float
      // float to int
      int16_t tmp = sample * 32768.0f;
      tmp = (tmp <= -32768) ? -32768 : (tmp >= 32767) ? 32767 : tmp;

      out[ii]=sample;

      //      out[ii] = ((float)delay_buffer[DELAY_SIZE-5] * alpha) + ((float)delay_buffer[DELAY_SIZE-6] * (1.0f - alpha));
        }
	else */// no filter 
    outt = ((float)delay_buffer[DELAY_SIZE-5] * alpha) + ((float)delay_buffer[DELAY_SIZE-6] * (1.0f - alpha));
    if (outt>LIMIN) outt=LIMIN;
    else if (outt<-LIMIN) outt=-LIMIN;
    out[ii]=outt;
        
    //out[ii] = delay_buffer[DELAY_SIZE-5];


  time_now += factor;
  last_time = int_time;
  int_time = time_now;
  while(last_time<int_time)      {
    //    doadc();
    int16_t val=getsample();
    //    new_data(val);
    delay_buffer[0] = delay_buffer[1];
    delay_buffer[1] = val;//data;

    last_time += 1;
  }
 }
}
