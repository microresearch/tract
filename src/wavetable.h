#ifndef __WAVE_H
#define __WAVE_H

#include "stdlib.h"
#include <stdint.h>


typedef struct _Wavetable {
    const float *wavetable;
    float basicIncrement;
    float currentPosition;
  int16_t length;
} Wavetable;

void wavetable_init(Wavetable* wavtable, const float *tableitself, int16_t length); // need to declare wavetable struct and and ourtable we use
void dowavetable(float* outgoing, Wavetable *wavetable, float frequency, unsigned char length);

//static const wormer waver={0, 1.0f, wave_get_sample, wave_newsay, 0, 0};
int16_t wave_get_sample1();
int16_t wave_get_sample2();
int16_t wave_get_sample3();
int16_t wave_get_sample4();
int16_t wave_get_sample5();
int16_t wave_get_sample6();
void wave_newsay();

#endif
