#include <stdio.h>
#include "math.h"
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/times.h>
#include <sys/unistd.h>
#include "LF.h"

extern float _selx, _sely, _selz;

#define PI            3.1415927
#define MY_PI 3.14159265359
#define TRUE 1
#define FALSE 0

LF LFxx;

// basic voice
/*
double _tcVal = 1.0;
double _teVal = 0.780;
double _tpVal = 0.6;
double _taVal = 0.028;
double vocalTension = 0.0;
int noiseOn = FALSE;
*/

// FRY?

double _tcVal = 1.0;
double _teVal = 0.251;
double _tpVal = 0.19;
double _taVal = 0.008;
//double noiseOn = FALSE;
double vocalTension = 0.0;

/* others

FRY:
double _tcVal = 1.0;
double _teVal = 0.251;
double _tpVal = 0.19;
double _taVal = 0.008;
double _noiseOn = FALSE;
double _vocalTension = 0.0;

FALSETTO:
double _tcVal = 1.0;
double _teVal = 0.770;
double _tpVal = 0.570;
double _taVal = 0.133;
double _noiseOn = TRUE;
double _noiseAmount = 0.015;
double _noiseDuration = 0.5;
double _noiseStart = 0.75;
double _vocalTension = 0.0;

BREATHY:
double _noiseOn = TRUE;
double _noiseDuration = 0.5;
double _noiseStart = 0.75;
double _noiseAmount = 0.025;
double _teVal = 0.756;
double _tpVal = 0.529;
double _taVal = 0.082;

MODAL:
double _noiseOn = FALSE;
double _teVal = 0.575;
double _tpVal = 0.457;
double _taVal = 0.028;

*/

void initLF(LF* lf);

// THIS IS ROSENBERG:
int16_t rosenberg_get_sample(){
  double LFcurrentSample1;
  signed int s16;
  float sample;
  float T=1/(100.0+(1000.0*_selz));
  float fs=48000;
  float pulselength=T*fs;

 // %N2 is duty cycle of the pulse, from 0 to 1.
 //%N1 is the duration of the glottal opening as a fraction of the 
 //%total pulse, from 0 to 1.

 float N2=pulselength*0.9;
 float N1=0.9*N2; // was 0.9
 static int16_t kkk;
 kkk++;
 if (kkk<=N1){
   LFcurrentSample1=0.5*(1-cos(MY_PI*(kkk-1)/N1));
   s16=(signed int)(LFcurrentSample1*32000.0);
 }
 else if (kkk>N1 && kkk<=N2){
   LFcurrentSample1=cos(MY_PI*(kkk-N1)/(N2-N1)/2);
   s16=(signed int)(LFcurrentSample1*32000.0);
 }
 else if (kkk>N2){
   LFcurrentSample1=0.0;
   s16=(signed int)(LFcurrentSample1);
   if (kkk>=pulselength) kkk=0;
   }
return s16*2;
}


int16_t LF_calc(LF* LFx){
  double LFcurrentSample;

LFx->t = (double)LFx->k*LFx->period/(double)LFx->dataLength;

if (LFx->t<LFx->te) {
LFcurrentSample = LFx->Eo*(exp(LFx->alpham*LFx->t)) * sin(LFx->wg*LFx->t);
}

if (LFx->t>=LFx->te) {
  LFcurrentSample = -((LFx->Ee)/(LFx->epsilon*LFx->ta))*(exp(LFx->epsilon*(LFx->t-LFx->te)) - exp(-LFx->epsilon*(LFx->tc-LFx->te))); // crazy exponential
}

if (LFx->t>LFx->tc) {
LFcurrentSample = 0.0;
LFx->k=0;
}

if (LFcurrentSample>32768.0) LFcurrentSample=32768.0;
if (LFcurrentSample<-32768.0) LFcurrentSample=-32768.0;

//   printf("%d\n",s16);
 unsigned int s16=(unsigned int)(LFcurrentSample*1.0);
 LFx->k += 1;
 //   fwrite(&s16,2,1,fo);
 return s16;
}

// This function calculates LF equation coefficients based on currently selected voice type
void initLF(LF* LFx){
 double f0;
int Fs;
int overSample;
double areaSum, area1, area2;
double optimumArea;
double epsilonTemp, epsilonDiff, epsilonOptimumDiff;
double tn,tb;

Fs = 8000;
f0 = 100.0;
overSample = 1000;
LFx->period = 1/f0;
LFx->Ee = - 1.0;
LFx->tc = _tcVal;
LFx->te = _teVal;
LFx->tppp = _tpVal;
LFx->ta = _taVal;

LFx->te = LFx->te + LFx->te*vocalTension;
LFx->tppp = LFx->tppp + LFx->tppp*vocalTension;
LFx->ta = LFx->ta - LFx->ta*vocalTension;

LFx->tc=LFx->tc*LFx->period;
LFx->te=LFx->te*LFx->period;
LFx->tppp=LFx->tppp*LFx->period;
LFx->ta=LFx->ta*LFx->period;

//These if statements prevent timing parameter values from going outside of their range with respect to other values.
if (LFx->te <= LFx->tppp) {
LFx->te = LFx->tppp + LFx->tppp*0.01;
}
if (LFx->te >= (LFx->tc-LFx->ta)){
LFx->te = LFx->tc-LFx->ta - (LFx->tc-LFx->ta)*0.01;
 }

// over sample values (can omit this section if it impacts real time operation)
/*period = period/overSample;
tc = tc/overSample;
te = te/overSample;
tppp = tppp/overSample;
ta = ta/overSample;*/
// dependant timing parameters
tn = LFx->te - LFx->tppp;
tb = LFx->tc - LFx->te;
// angular frequency of sinusoid section
LFx->wg = MY_PI/LFx->tppp;
// maximum negative peak value
LFx->Eo = LFx->Ee;
// epsilon and alpha equation coefficients


areaSum = 1.0;
optimumArea = 1e-14;
epsilonDiff = 10000.0;
epsilonOptimumDiff = 0.1;
epsilonTemp = 1/LFx->ta;
// solve iteratively for epsilon
while (fabs(epsilonDiff)>epsilonOptimumDiff) {
LFx->epsilon = (1/LFx->ta)*(1-exp(-epsilonTemp*tb));
epsilonDiff = LFx->epsilon - epsilonTemp;
epsilonTemp = LFx->epsilon;
if (epsilonDiff<0) {
epsilonTemp = epsilonTemp + (fabs(epsilonDiff)/100);
}
if (epsilonDiff>0) {
epsilonTemp = epsilonTemp - (fabs(epsilonDiff)/100);
}
}

//printf("Epislon %f\n", epsilon);

// iterate through area balance to get Eo and alpha to give area1 + area2 = 0
//while ((areaSum< -optimumArea)||(areaSum>optimumArea)) {

while (areaSum > optimumArea){

  //printf("Ee %f Eo %f wg %f te %f Areasum %f\n",Ee,Eo,wg,te,area2);

LFx->alpham = ( log(-LFx->Ee/(LFx->Eo*sin(LFx->wg*LFx->te))))/LFx->te; // problem is here! NaN for log of negative -> clog
//alpha = ((-Ee/(Eo*sin(wg*te))))/te;

//printf("al %f xx %f xx %f\n",alpha, te,sin(wg*te) );

area1 = ( LFx->Eo*exp(LFx->alpham*LFx->te)/(sqrt(LFx->alpham*LFx->alpham+LFx->wg*LFx->wg))) * (sin(LFx->wg*LFx->te-atan(LFx->wg/LFx->alpham))) + (LFx->Eo*LFx->wg/(LFx->alpham*LFx->alpham+LFx->wg*LFx->wg));

area2 = ( -(LFx->Ee)/(LFx->epsilon*LFx->epsilon*(LFx->ta))) * (1 - exp(-LFx->epsilon*tb*(1+LFx->epsilon*tb)));
areaSum = area1 + area2;


if (areaSum>0.0) {
LFx->Eo = LFx->Eo - 1e5*areaSum;
}
if (areaSum<0.0) {
LFx->Eo = LFx->Eo + 1e5*areaSum;
}
}


// alpham=938077.0;
// epsilon=12500000.0;

// alpham=305670.0;
// epsilon=3355800.0;


//calculate length of waveform in samples
//if (_pitchSlide == FALSE) {
//dataLength = floor(overSample*Fs*period);
//}
//else{
//dataLength = _dataLength;
//}
//calculate length of waveform in samples without overSample
//
LFx->dataLength = floor(Fs*LFx->period);
// pass all variables needed for waveform calculation to effectState

//printf("period=%f, length= %d, alpham = %f epsilon= %f \n", period, dataLength, alpham,epsilon);

//printf("period=%f, length= %d, tc = %f, te = %f, tp = %f, ta = %f wg = %f \n", period, dataLength, tc, te, tp,ta, wg);
}

// This function calculates LF equation coefficients based on currently selected voice type
void initLFF(LF* LFx){
 double f0;
int Fs;
int overSample;
double areaSum, area1, area2;
double optimumArea;
double epsilonTemp, epsilonDiff, epsilonOptimumDiff;
double tn,tb;

Fs = 8000;
 f0 = 20.0 +(400.0*_selz);
overSample = 1000;
LFx->period = 1/f0;
LFx->Ee = - 1.0;
LFx->tc = _tcVal;
LFx->te = _teVal;
LFx->tppp = _tpVal;
LFx->ta = _taVal;

LFx->te = LFx->te + LFx->te*vocalTension;
LFx->tppp = LFx->tppp + LFx->tppp*vocalTension;
LFx->ta = LFx->ta - LFx->ta*vocalTension;

LFx->tc=LFx->tc*LFx->period;
LFx->te=LFx->te*LFx->period;
LFx->tppp=LFx->tppp*LFx->period;
LFx->ta=LFx->ta*LFx->period;

//These if statements prevent timing parameter values from going outside of their range with respect to other values.
if (LFx->te <= LFx->tppp) {
LFx->te = LFx->tppp + LFx->tppp*0.01;
}
if (LFx->te >= (LFx->tc-LFx->ta)){
LFx->te = LFx->tc-LFx->ta - (LFx->tc-LFx->ta)*0.01;
 }

// over sample values (can omit this section if it impacts real time operation)
/*period = period/overSample;
tc = tc/overSample;
te = te/overSample;
tppp = tppp/overSample;
ta = ta/overSample;*/
// dependant timing parameters
tn = LFx->te - LFx->tppp;
tb = LFx->tc - LFx->te;
// angular frequency of sinusoid section
LFx->wg = MY_PI/LFx->tppp;
// maximum negative peak value
LFx->Eo = LFx->Ee;
// epsilon and alpha equation coefficients


areaSum = 1.0;
optimumArea = 1e-14;
epsilonDiff = 10000.0;
epsilonOptimumDiff = 0.1;
epsilonTemp = 1/LFx->ta;
// solve iteratively for epsilon
while (fabs(epsilonDiff)>epsilonOptimumDiff) {
LFx->epsilon = (1/LFx->ta)*(1-exp(-epsilonTemp*tb));
epsilonDiff = LFx->epsilon - epsilonTemp;
epsilonTemp = LFx->epsilon;
if (epsilonDiff<0) {
epsilonTemp = epsilonTemp + (fabs(epsilonDiff)/100);
}
if (epsilonDiff>0) {
epsilonTemp = epsilonTemp - (fabs(epsilonDiff)/100);
}
}

//printf("Epislon %f\n", epsilon);

// iterate through area balance to get Eo and alpha to give area1 + area2 = 0
//while ((areaSum< -optimumArea)||(areaSum>optimumArea)) {

while (areaSum > optimumArea){

  //printf("Ee %f Eo %f wg %f te %f Areasum %f\n",Ee,Eo,wg,te,area2);

LFx->alpham = ( log(-LFx->Ee/(LFx->Eo*sin(LFx->wg*LFx->te))))/LFx->te; // problem is here! NaN for log of negative -> clog
//alpha = ((-Ee/(Eo*sin(wg*te))))/te;

//printf("al %f xx %f xx %f\n",alpha, te,sin(wg*te) );

area1 = ( LFx->Eo*exp(LFx->alpham*LFx->te)/(sqrt(LFx->alpham*LFx->alpham+LFx->wg*LFx->wg))) * (sin(LFx->wg*LFx->te-atan(LFx->wg/LFx->alpham))) + (LFx->Eo*LFx->wg/(LFx->alpham*LFx->alpham+LFx->wg*LFx->wg));

area2 = ( -(LFx->Ee)/(LFx->epsilon*LFx->epsilon*(LFx->ta))) * (1 - exp(-LFx->epsilon*tb*(1+LFx->epsilon*tb)));
areaSum = area1 + area2;


if (areaSum>0.0) {
LFx->Eo = LFx->Eo - 1e5*areaSum;
}
if (areaSum<0.0) {
LFx->Eo = LFx->Eo + 1e5*areaSum;
}
}


// alpham=938077.0;
// epsilon=12500000.0;

// alpham=305670.0;
// epsilon=3355800.0;


//calculate length of waveform in samples
//if (_pitchSlide == FALSE) {
//dataLength = floor(overSample*Fs*period);
//}
//else{
//dataLength = _dataLength;
//}
//calculate length of waveform in samples without overSample
//
LFx->dataLength = floor(Fs*LFx->period);
// pass all variables needed for waveform calculation to effectState

//printf("period=%f, length= %d, alpham = %f epsilon= %f \n", period, dataLength, alpham,epsilon);

//printf("period=%f, length= %d, tc = %f, te = %f, tp = %f, ta = %f wg = %f \n", period, dataLength, tc, te, tp,ta, wg);
}





int16_t LF_get_sample(){
  initLFF(&LFxx);
  int16_t tmp=LF_calc(&LFxx);
  return tmp;
}


/////
/*
%Rosenberg Pulse
%this function accepts fundamental frequency of the glottal signal and 
%the sampling frequency in hertz as input and returns one period of 
%the rosenberg pulse at the specified frequency.
%N2 is duty cycle of the pulse, from 0 to 1.
%N1 is the duration of the glottal opening as a fraction of the 
%total pulse, from 0 to 1.
function[gn]=rosenberg(N1,N2,f0,fs)
T=1/f0;     %period in seconds
pulselength=floor(T*fs);    %length of one period of pulse
%select N1 and N2 for duty cycle
N2=floor(pulselength*N2);
N1=floor(N1*N2);
gn=zeros(1,N2);
%calculate pulse samples
for n=1:N1-1
    gn(n)=0.5*(1-cos(pi*(n-1)/N1));
end
for n=N1:N2
    gn(n)=cos(pi*(n-N1)/(N2-N1)/2);
end
gn=[gn zeros(1,(pulselength-N2))];
*/
