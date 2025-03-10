#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <main.h>
#include "process.h"
#include "wavetable.h"
//#include "wvocoder.h"
//#include "svf.h"

extern Wavetable wavtable;


void processe(float *flinbuffer,float *frinbuffer,float *foutbuffer, int16_t sz)
{
  uint16_t x;
  float modulation[AUDIO_BLOCK_SIZE];
 
  // generate source
  //      dowavetable(foutbuffer, &wavtable, 200.0f, sz); 
  //  dowavetable(modulation, &wavtable, 200.0f, sz); 
  // process modulation source and/or incoming
  /*
void Vocoder_Process(
    const float* modulator,
    const float* carrier, // carrier is modulation source
    float* out,
    u8 size) {
  */
  //    Vocoder_Process(flinbuffer,modulation,foutbuffer,sz);
  //    runBANDStest_(flinbuffer, foutbuffer, sz); // workings
            for (x=0;x<sz;x++){
	      foutbuffer[x]=frinbuffer[x]; // copy works
	      //	      foutbuffer[x]=modulation[x]; // copy works
      }
}

