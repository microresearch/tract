/*
    Copyright (C) 2017 Martin Howse
    Copyright (c) 1994,2001-2003 Nick Ing-Simmons. All rights reserved.
 
    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public
    License along with this library; if not, write to the Free
    Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
    MA 02111-1307, USA

*/


//#include "stm32f4xx.h"

#include "stdint.h"
#include <stdio.h>

#include "config.h"
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include "useconfig.h"
#include <math.h>
#include "rsynth.h"
#include "resources.h"


typedef unsigned char UINT8;
typedef signed char INT8;
typedef unsigned char u8;
typedef uint16_t UINT16;
typedef uint16_t u16;
typedef int16_t INT16;
typedef uint32_t UINT32;
typedef int32_t INT32;


float ampl;
extern float _selx, _sely, _selz;
extern u8 test_elm_rsynthy[106]; // as is just phon code and length... with lat stop



#define MAXED(var, max) \
  if (var > (max)) {	\
    var = (max);	\
  }

typedef struct {
    float v;			/* boundary value */
    long t;			/* transition time */
} slope_t;

typedef struct {
    float a;
    float b;
    float v;
} filter_t, *filter_ptr;

static float
filter(filter_ptr p, float v)
{
    return p->v = (p->a * v + p->b * p->v);
}

/* 'a' is dominant element, 'b' is dominated
   ext is flag to say to use external times from 'a' rather
   than internal i.e. ext != 0 if 'a' is NOT current element.

 */

static void
set_trans(slope_t * t, int i, Elm_ptr a, Elm_ptr b, int ext, int e,
	  float speed)
{
    t[i].t = (long) (((ext) ? a->p[i].ed : a->p[i].id) * speed);
    if (t[i].t) {
	float afrac = a->p[i].prop * 0.01F;
	t[i].v = a->p[i].stdy * (1.0F - afrac) + (afrac * b->p[i].stdy);
    }
    else
	t[i].v = b->p[i].stdy;
}


/*
   ______________ b
   /
   /
   /
   a____________/
   0   d
   ---------------t---------------
 */

float
linear(float a, float b, long t, long d)
{
    if (t <= 0)
	return a;
    else if (t >= d)
	return b;
    else {
	float f = (float) t / (float) d;
	return a + (b - a) * f;
    }
}

float
interpolate(char *w, char *p, slope_t * s, slope_t * e, float mid, long t,
	    long d)
{
    float steady = d - (s->t + e->t);
    if (steady >= 0) {
	/* Value reaches stready state somewhere ... */
	if (t < s->t)
	    return linear(s->v, mid, t, s->t);	/* initial transition */
	else {
	    t -= s->t;
	    if (t <= steady)
		return mid;	/* steady state */
	    else
		return linear(mid, e->v, (int) (t - steady), e->t);
	    /* final transition */
	}
    }
    else {
	float f = (float) 1.0 - ((float) t / (float) d);
	float sp = linear(s->v, mid, t, s->t);
	float ep = linear(e->v, mid, d - t, e->t);
	return f * sp + ((float) 1.0 - f) * ep;
    }
}


// rsynth and rest as fixed:

//extern rsynth_t rsynthi;
extern rsynth_t rsynthi;
static rsynth_t * rsynth=&rsynthi; // where do we fix this
static float const *f0; // >>>???
static float *ff0;
//static unsigned char nf0; // >>>???
static filter_t flt[nEparm];
static float f0s, f0e;
static unsigned int tf0 = 0;
static unsigned int ntf0 = 0;
static unsigned char i = 0;
static unsigned char nextelement=1;

typedef struct{
float f0[3];
 unsigned char elements[52];
 unsigned char length;
} vocab_t_;

#include "rsynthy_vocab.h"

static unsigned char nelm=0;
static const unsigned char *elm;

// for single phonemes - 2 elements and always shift out???

static u8 element[4];
static float curpitch;


void rsynth_newsay_single(u8 once)
{
    ampl=0.7f;

  u8 selected=1;
  //  static u8 once=1;
  nelm=4; //0123

  u8 val=_selz*87.0f;
  MAXED(val,83);
  //  val=83-val; // there are 84 elements // DON't INVERT?

  element[0]=element[2];
  element[1]=element[3];
  element[2]=val;
  element[3]=((_sely)*32.0f)+1; // length
  
  //  unsigned i = 0;
  int j;

  if (once==1){ // question of this ONCE? - but we should reset or?
    //    once=0;
  Elm_ptr le = &Elementz[0]; 
  
  for (j = 0; j < nEparm; j++) {
	flt[j].v = le->p[j].stdy;
	flt[j].a = rsynth->smooth;
	flt[j].b = 1.0F - rsynth->smooth;
	}
    }
    i=0;
    tf0=0;
    ntf0=0;
    nextelement=1;
}

void rsynth_newsay_singlex(){
  rsynth_newsay_single(1);
}


void rsynth_newsay()
{
  //  u8 selected=1;
  ampl=1.0f;
  u8 selected=_selz*131.0f;
  MAXED(selected,127);
  //  selected=121;
  selected=127-selected;
  //  selected=126;
  nelm=rsynth_vocab[selected]->length;
  elm=rsynth_vocab[selected]->elements;
  f0=rsynth_vocab[selected]->f0;

  //  unsigned i = 0;
  f0s = rsynth->speaker->F0Hz;
  f0e = f0s;
  int j;
  const Elm_ptr le = &Elementz[0]; 
  
  f0e = f0s = *f0++;
    for (j = 0; j < nEparm; j++) {
	flt[j].v = le->p[j].stdy;
	flt[j].a = rsynth->smooth;
	flt[j].b = 1.0F - rsynth->smooth;
    }

    i=0;
    tf0=0;
    ntf0=0;
    nextelement=1;
}

float contour[3];//={146.300003, 108.000000, 133.339996};


void rsynth_newsay_elm() 
{

  nelm=106; // length
  ampl=1.0f;

  // do f0 contour?
  unsigned char ii;
  ff0 = contour;
  //  nf0 = 3;
  //  ff0[1]=0;
  /*  for (ii = 0; ii < nelm; ii += 2) {
    ff0[1] += test_elm_rsynthy[ii + 1]; // length in frames
    }*/
  //  ff0[1] = 0.6f * ff0[0];	/* bottom */

  //  ff0[0] = 1.1f * rsynth->speaker->F0Hz;	/* top */
  //    ff0[1] = 1.1f * rsynth->speaker->F0Hz;	/* top */
  //  ff0[2] = 0.6f * ff0[0];	/* bottom */

  f0s = rsynth->speaker->F0Hz;
  f0e = f0s;
  int j;
  const Elm_ptr le = &Elementz[0];
  
  //  f0e = f0s = *ff0++;
    for (j = 0; j < nEparm; j++) {
	flt[j].v = le->p[j].stdy;
	flt[j].a = rsynth->smooth;
	flt[j].b = 1.0f - rsynth->smooth;
    }

    i=0; // start of new elements
    //    tf0=0;
    //    ntf0=0;
    nextelement=1;
}


int16_t rsynth_get_sample(){

  static short samplenumber=0;
  static unsigned char newframe=0, dur;
  int16_t sample;
  static float ep[nEparm];
  static Elm_ptr le = &Elementz[0];
  float speed = rsynth->speed;
  static float F0Hz;
  static Elm_ptr ce;
  static slope_t start[nEparm];
  static slope_t end[nEparm];
  static unsigned t=0;
  
  if (i>nelm && nextelement==1){   // NEW utterance which means we hit nelm=0 in our cycling:
    rsynth_newsay();
  }


  if (nextelement==1){ // we have a new element
    samplenumber=0; 
    ce = &Elementz[elm[i++]];
    dur = elm[i++];
	/* Skip zero length elements which are only there to affect
	   boundary values of adjacent elements.
	   Note this mainly refers to "QQ" element in fricatives
	   as stops have a non-zero length central element.
	 */
	if (dur > 0) {
	    Elm_ptr ne = (i < nelm) ? &Elementz[elm[i]] : &Elementz[0];

	    int ii;
	    for (ii = 0; ii < nEparm; ii++) {

		if (ce->p[ii].rk > le->p[ii].rk) {
		    set_trans(start, ii, ce, le, 0, 's', speed);
		    /* we dominate last */
		}
		else {
 		    set_trans(start, ii, le, ce, 1, 's', speed);
		    /* last dominates us */
		}

		if (ne->p[ii].rk > ce->p[ii].rk) {
		    set_trans(end, ii, ne, ce, 1, 'e', speed);
		    /* next dominates us */
		}
		else {
		    set_trans(end, ii, ce, ne, 0, 'e', speed);
		    /* we dominate next */
		}
	    }
	    t=0;
	    newframe=1;
	}
	else {
	  // dur==0???
	  //	    newframe=1;
	  rsynth_get_sample();
	}
  } // end of next element stuff   
	//	    for (t = 0; t < dur; t++, tf0++) {

  if (newframe==1) { // this is a new frame - so we need new parameters
    newframe=0;
    // inc and are we at end of frames in which case we need next element?

    if (t<dur){ //
      int j;
      //      float peak = 0.25f;

		for (j = 0; j < nEparm; j++) {
		    ep[j] =
			filter(flt + j,
			       interpolate(ce->name, Ep_namez[j], &start[j],
					   &end[j], (float) ce->p[j].stdy,
					   t, dur));
		}

		while (tf0 == ntf0) {
		  tf0 = 0;
		  f0s = f0e;
		  ntf0 = (unsigned) *f0++;
		  f0e = *f0++;
		  }

		/* interpolate the f0 value */
		F0Hz = linear(f0s, f0e, tf0, ntf0);
		nextelement=0;
		t++; tf0++;
    } //dur

    else { // hit end of DUR number of frames...
      nextelement=1; 
      le = ce; // where we can put this?????? 
      rsynth_get_sample();
    }
  } // newframe ==1


  int val=_selx*1028.0f;
  MAXED(val,1023);
  val=1023-val;
  float newfreq=F0Hz* logspeed[val]; // was * 0.5f
    sample=rsynth_frame_single(rsynth, newfreq, ep);
  //  sample=rand()%32768;
      samplenumber++;

      val=_sely*131.0f;
      MAXED(val,127);
      //      val=127-val;
      
      if (samplenumber>=rsynth->samples_frame*logpitch[val]) { // how many in a frame??? 256 for 32000 samplerate
      // end of frame so...????
      newframe=1;
      samplenumber=0;
      }
//  }
      return (int)(sample);
}


int16_t rsynth_get_sample_sing(){

  static short samplenumber=0;
  static unsigned char newframe=0, dur;
  int16_t sample;
  static float ep[nEparm];
  Elm_ptr le = &Elementz[0];
  float speed = rsynth->speed;
  static float F0Hz;
  static Elm_ptr ce;
  static slope_t start[nEparm];
  static slope_t end[nEparm];
  static unsigned t=0;
  
  if (i>nelm && nextelement==1){   // NEW utterance which means we hit nelm=0 in our cycling:
    rsynth_newsay();
  }

  if (nextelement==1){ // we have a new element
    samplenumber=0; 

    ce = &Elementz[elm[i++]];
    dur = elm[i++];
	/* Skip zero length elements which are only there to affect
	   boundary values of adjacent elements.
	   Note this mainly refers to "QQ" element in fricatives
	   as stops have a non-zero length central element.
	 */
	if (dur > 0) {
	    Elm_ptr ne = (i < nelm) ? &Elementz[elm[i]] : &Elementz[0];

	    int ii;
	    for (ii = 0; ii < nEparm; ii++) {

		if (ce->p[ii].rk > le->p[ii].rk) {
		    set_trans(start, ii, ce, le, 0, 's', speed);
		    /* we dominate last */
		}
		else {
 		    set_trans(start, ii, le, ce, 1, 's', speed);
		    /* last dominates us */
		}

		if (ne->p[ii].rk > ce->p[ii].rk) {
		    set_trans(end, ii, ne, ce, 1, 'e', speed);
		    /* next dominates us */
		}
		else {
		    set_trans(end, ii, ce, ne, 0, 'e', speed);
		    /* we dominate next */
		}
	    }
	    t=0;
	    newframe=1;
	}
	else {
	  // dur==0???
	  //	    newframe=1;
	  rsynth_get_sample_sing();
	}
  } // end of next element stuff   
	//	    for (t = 0; t < dur; t++, tf0++) {

  if (newframe==1) { // this is a new frame - so we need new parameters
    newframe=0;
    // inc and are we at end of frames in which case we need next element?

    if (t<dur){ //
      int j;
      //      float peak = 0.25f;

		for (j = 0; j < nEparm; j++) {
		    ep[j] =
			filter(flt + j,
			       interpolate(ce->name, Ep_namez[j], &start[j],
					   &end[j], (float) ce->p[j].stdy,
					   t, dur));
		}

		/*		while (tf0 == ntf0) { // unused for singing
		  tf0 = 0;
		  f0s = f0e;
		  ntf0 = (unsigned) *f0++;
		  f0e = *f0++;
		  }*/

		// interpolate the f0 value 
		//		F0Hz = linear(f0s, f0e, tf0, ntf0);
		nextelement=0;
		t++; tf0++;
    } //dur

    else { // hit end of DUR number of frames...
      nextelement=1; 
      le = ce; // where we can put this??????
      rsynth_get_sample_sing();
    }
  } // newframe ==1


  int val=_selx*1028.0f;
  MAXED(val,1023);
  val=1023-val;
  //   float newfreq=F0Hz* logspeed[val] * 0.5f;
  float newfreq=rsynth->speaker->F0Hz * logspeed[val];// * 0.5f;
    sample=rsynth_frame_single(rsynth, newfreq, ep);
  //  sample=rand()%32768;
      samplenumber++;

      val=_sely*131.0f;
      MAXED(val,127);
      //      val=127-val;
      
      if (samplenumber>=rsynth->samples_frame*logpitch[val]) { // how many in a frame??? 256 for 32000 samplerate
      // end of frame so...????
      newframe=1;
      samplenumber=0;
      }
//  }
      return (int)(sample);
}

int16_t rsynth_get_sample_single(){

  static short samplenumber=0;
  static unsigned char newframe=0, dur;
  int16_t sample;
  static float ep[nEparm];
  Elm_ptr le = &Elementz[0];
  float speed = rsynth->speed;
  static float F0Hz;
  static Elm_ptr ce;
  static slope_t start[nEparm];
  static slope_t end[nEparm];
  static unsigned t=0;

  if (i>1 && nextelement==1){   // NEW utterance which means we hit nelm=0 in our cycling:
    rsynth_newsay_single(0);
  }

  if (nextelement==1){ // we have a new element
    samplenumber=0; 

    ce = &Elementz[element[i++]];
    dur = element[i++];
	/* Skip zero length elements which are only there to affect
	   boundary values of adjacent elements.
	   Note this mainly refers to "QQ" element in fricatives
	   as stops have a non-zero length central element.
	 */
	if (dur > 0) {
	    Elm_ptr ne = (i < nelm) ? &Elementz[element[i]] : &Elementz[0];

	    int ii;
	    for (ii = 0; ii < nEparm; ii++) {

		if (ce->p[ii].rk > le->p[ii].rk) {
		    set_trans(start, ii, ce, le, 0, 's', speed);
		    /* we dominate last */
		}
		else {
 		    set_trans(start, ii, le, ce, 1, 's', speed);
		    /* last dominates us */
		}

		if (ne->p[ii].rk > ce->p[ii].rk) {
		    set_trans(end, ii, ne, ce, 1, 'e', speed);
		    /* next dominates us */
		}
		else {
		    set_trans(end, ii, ce, ne, 0, 'e', speed);
		    /* we dominate next */
		}
	    }
	    t=0;
	    newframe=1;
	}
	else {
	  // dur==0???
	  //	    newframe=1;
	  rsynth_get_sample_single();
	}
  } // end of next element stuff   
	//	    for (t = 0; t < dur; t++, tf0++) {

  if (newframe==1) { // this is a new frame - so we need new parameters
    newframe=0;
    // inc and are we at end of frames in which case we need next element?

    if (t<dur){ //
      int j;
      //      float peak = 0.25f;

		for (j = 0; j < nEparm; j++) {
		    ep[j] =
			filter(flt + j,
			       interpolate(ce->name, Ep_namez[j], &start[j],
					   &end[j], (float) ce->p[j].stdy,
					   t, dur));
		}

		nextelement=0;
		t++;
    } //dur

    else { // hit end of DUR number of frames...
      nextelement=1; 
      le = ce; // where we can put this?????? 
      rsynth_get_sample_single();
    }
  } // newframe ==1


  int val=_selx*1027.0f;
  MAXED(val,1023);
  val=1023-val;
    float newfreq=rsynth->speaker->F0Hz* logspeed[val];//rsynth->speaker->F0Hz
  //float newfreq=F0Hz* logspeed[val] * 0.5f;

    sample=rsynth_frame_single(rsynth, newfreq, ep);
      samplenumber++;
      
      if (samplenumber>=rsynth->samples_frame) { // how many in a frame??? 256 for 32000 samplerate
      // end of frame so...????
      newframe=1;
      samplenumber=0;
      }
//  }
      return (int)(sample);
}

int16_t rsynth_get_sample_elm(){

  static short samplenumber=0;
  static unsigned char newframe=0, dur;
  int16_t sample;
  static float ep[nEparm];
  Elm_ptr le = &Elementz[0];
  float speed = rsynth->speed;
  static float F0Hz;
  static Elm_ptr ce;
  static slope_t start[nEparm];
  static slope_t end[nEparm];
  static unsigned t=0;
  //  extent nextent={51,53.0f};

  u8 xaxis=_selx*56.0f;
  MAXED(xaxis,52); 
  xaxis=52-xaxis;
  u8 val=_selz*87.0f;
  MAXED(val,83);
  //  val=83-val; // DON'T INVERT?
  test_elm_rsynthy[xaxis*2]=val+1; // xaxis must be 
  test_elm_rsynthy[(xaxis*2)+1]=(_sely*33.0f)+2; // length say max 32 - long to short left to right

  
  if (i>nelm && nextelement==1){   // NEW utterance which means we hit nelm=0 in our cycling:
    rsynth_newsay_elm();
  }


  if (nextelement==1){ // we have a new element
    samplenumber=0; 

    ce = &Elementz[test_elm_rsynthy[i++]];
    dur = test_elm_rsynthy[i++];
	/* Skip zero length elements which are only there to affect
	   boundary values of adjacent elements.
	   Note this mainly refers to "QQ" element in fricatives
	   as stops have a non-zero length central element.
	 */
	if (dur > 0) {
	  Elm_ptr ne = (i < nelm) ? &Elementz[test_elm_rsynthy[i]] : &Elementz[0];

	    int ii;
	    for (ii = 0; ii < nEparm; ii++) {

		if (ce->p[ii].rk > le->p[ii].rk) {
		    set_trans(start, ii, ce, le, 0, 's', speed);
		    /* we dominate last */
		}
		else {
 		    set_trans(start, ii, le, ce, 1, 's', speed);
		    /* last dominates us */
		}

		if (ne->p[ii].rk > ce->p[ii].rk) {
		    set_trans(end, ii, ne, ce, 1, 'e', speed);
		    /* next dominates us */
		}
		else {
		    set_trans(end, ii, ce, ne, 0, 'e', speed);
		    /* we dominate next */
		}
	    }
	    t=0;
	    newframe=1;
	}
	else {
	  // dur==0???
	  //	    newframe=1;
	  rsynth_get_sample_elm();
	}
  } // end of next element stuff   
	//	    for (t = 0; t < dur; t++, tf0++) {

  if (newframe==1) { // this is a new frame - so we need new parameters
    newframe=0;
    // inc and are we at end of frames in which case we need next element?

    if (t<dur){ //
      int j;
      //      float peak = 0.25f;

		for (j = 0; j < nEparm; j++) {
		    ep[j] =
			filter(flt + j,
			       interpolate(ce->name, Ep_namez[j], &start[j],
					   &end[j], (float) ce->p[j].stdy,
					   t, dur));
		}

		//		while (tf0 == ntf0) { // first time?
      // tf0 = 0;
		  //		  f0s = f0e;
		//	  ntf0 = (unsigned int) *ff0++;
		  //		  f0e = *ff0++;
		//  }

		/* interpolate the f0 value */
		//		F0Hz = linear(f0s, f0e, tf0, ntf0);
		F0Hz = rsynth->speaker->F0Hz;
		nextelement=0;
		t++;
		// tf0++;
    } //dur

    else { // hit end of DUR number of frames...
      nextelement=1; 
      le = ce; // where we can put this?????? 
      rsynth_get_sample_elm();
    }
  } // newframe ==1

  sample=rsynth_frame_single(rsynth, F0Hz, ep);
  //  sample=rand()%32768;
      samplenumber++;
      
      if (samplenumber>=rsynth->samples_frame) { // how many in a frame??? 256 for 32000 samplerate
      // end of frame so...????
      newframe=1;
      samplenumber=0;
      }
//  }
      return (int)(sample);
}
