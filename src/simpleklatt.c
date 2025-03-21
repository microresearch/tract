/*
 * Copyright (C) 2017 Martin Howse
 * Copyright (C) 2013 Reece H. Dunn
 * (c) 1993,94 Jon Iles and Nick Ing-Simmons
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "audio.h"
#include "stdlib.h"
#include "stdint.h"
#include <stdio.h>
#include "parwave.h"
#include "resources.h"

/* for default sampled glottal excitation waveform */

#define NUMBER_OF_SAMPLES 100
//#define SAMPLE_FACTOR 0.00001 // this is too low - hardcoded in parwave

static   int16_t *iwave;
static   int16_t isam;
static   int16_t icount;
static   int16_t nmspf_def;
static   klatt_global_ptrr globals;
//  klatt_frame_ptrr frame;
static   unsigned char high_byte;
static   unsigned char low_byte;
static   flag raw_flag;
static   flag raw_type;

static int16_t frame[40];
extern float exy[64][240];
extern u8 _intmode;
extern float _selx, _sely, _selz;

static const int16_t natural_samples[NUMBER_OF_SAMPLES]= // unused source
  {
    -310,-400,530,356,224,89,23,-10,-58,-16,461,599,536,701,770,
    605,497,461,560,404,110,224,131,104,-97,155,278,-154,-1165,
    -598,737,125,-592,41,11,-247,-10,65,92,80,-304,71,167,-1,122,
    233,161,-43,278,479,485,407,266,650,134,80,236,68,260,269,179,
    53,140,275,293,296,104,257,152,311,182,263,245,125,314,140,44,
    203,230,-235,-286,23,107,92,-91,38,464,443,176,98,-784,-2449,
    -1891,-1045,-1600,-1462,-1384,-1261,-949,-730
  };

// check parameters: see klattparams - could be good to have contour morphing (snap on trigger)
// param0 is fund freq, param39 is volume

static const int16_t mins[40] __attribute__ ((section (".flash"))) = {200,  0, 200, 40, 550, 40, 1200, 40, 1200, 40, 1200, 40, 1200, 40, 248, 40, 248, 40, 0, 10, 0, 0, 0, 0, 0, 40, 0, 40, 0, 40, 0, 40, 0, 40, 0, 40, 0, 0, 0, 0};

//static const int16_t maxs[40] __attribute__ ((section (".flash"))) = {4000, 70, 1300, 1000, 3000, 1000, 4999, 1000, 4999, 1000, 4999, 1000, 4999, 2000, 528, 1000, 528, 1000, 70, 65, 80, 24, 80, 40, 80, 1000, 80, 1000, 80, 1000, 80, 1000, 80, 1000, 80, 2000, 80, 80, 70, 60};



static const int16_t range[40] __attribute__ ((section (".flash"))) ={3800, 70, 1100, 960, 2450, 960, 3799, 960, 3799, 960, 3799, 960, 3799, 1960, 280, 960, 280, 960, 40, 55, 40, 20, 40, 40, 80, 960, 80, 960, 80, 960, 80, 960, 80, 960, 80, 1960, 80, 40, 40, 40}; // changed some of range 

// after 280 is nasal bw, asp ampX-reduced, 

static klatt_global_tt globale;


void simpleklatt_init(void){

  //  straightwormy=addworm(10.0f,10.0f,100.0f, 100.0f, straightworm);

  globals=&globale;

  //  globals = (klatt_global_ptrr)malloc(sizeof(klatt_global_tt));
  //  frame = (klatt_frame_ptrr)malloc(sizeof(klatt_frame_tt));
  //  framer framezz[40];
  //  frame_init(globals,simpleklattset.val); 

  globals->synthesis_model = 0; // all_parallel=1 cascade=0
 globals->samrate = 16000;
 globals->glsource = 1; // 1=impulsive 2=glottal impulse 3=sampled as above =1
 globals->natural_samples = natural_samples;
 globals->num_samples = NUMBER_OF_SAMPLES;
//  globals->sample_factor = (float) SAMPLE_FACTOR;
 nmspf_def = 10;
  globals->nfcascade = 6;
  globals->f0_flutter = 0; // fix this!

unsigned char y;

// or this could be 32 for audio.c frame

 globals->nspfr = (globals->samrate * nmspf_def) / 100; // was / 1000// number of samples per frame = 320000 /1000 = 320
simple_parwave_init(globals);
}

void generate_exy_frame(int16_t* frame){
unsigned char y;

    for (y=0;y<39;y++){
      frame[y]=mins[y] + (range[y]*(1.0f-exy[0][y])); 
    }
    frame[39]=42; // FIXED volume - higher if we revert to all parallel
      }

static uint16_t samplenumber=0;


int16_t simpleklatt_get_sample(){
  u8 x=0;
  int16_t sampel;

  sampel=single_single_parwave(globals, frame);

  samplenumber++;
  if (samplenumber>globals->nspfr*_selz*4.0f) { // greater than what???? 320 samples - this can be our selz???? 
      // end of frame so...????
    samplenumber=0;
    simpleklatt_newsay();
    }
    return sampel;
}

void simpleklatt_newsay(){
  // generate the frame from our exy -> frame
  generate_exy_frame(frame);
  frame_init(globals,frame); 
    samplenumber=0;
  //    if (globals->f0_flutter != 0)
  //      flutter(globals,frame);  
    globals->ns=0;
}


