/* 

TODO: finish shimmer, all needs cleaning so does arbitrary numbers of samples in each case

LF is from: gm.py

Rosenberg is from lfgen.c

KLGLOTT88 is from glottalair.py

shimmer is from err flowgen_shimmer.c

 */

#include <errno.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/times.h>
#include "LF.h"

typedef unsigned int u16;
extern float _selx, _sely, _selz;

#define PI 3.14159265359
#define TRUE 1
#define FALSE 0

//FILE* fo;// = fopen("testlfgen.pcm", "wb");

signed short truncate(float aux)
{
  signed short i;
  if(aux > 32767) i = 32767;
  else if(aux < -32767) i = -32767;
  else i = (signed short) ceil(aux);

  return i;

}

struct PAR parry=  {1.0, 10.5, 0.55, 0.65, 125, 100, 0.0, 0.0, 8000L, 12000, 0.0, 0.5, 10.5};
SHIM shimm;
klglott glotty;

void shimmer_init(PAR* param, SHIM* shim){
  //  srandom(time(NULL));
  shim->P = shim->T = (int) ( (float) param->fs/param->F0);
}

int do_shimmer(PAR* param, SHIM* shim){
  uint16_t i;
    if(param->jitter != -1) {

      shim->DeltaPer[1] = shim->DeltaPer[0];
      do {
	shim->J = (  rand() / (RAND_MAX * 10000.0) ) * 40000.0 * param->jitter -
	    2.0*param->jitter;

	shim->DeltaPer[0] = shim->DeltaPer[1]*(2.0 + shim->J)/(2.0 - shim->J) +
		      2.0*shim->P*shim->J/(2.0 - shim->J);

	shim->T = (signed short) ceil((float) shim->P + shim->DeltaPer[0]);
      } while ( (float) shim->T > (float) 1.2*shim->P || (float) shim->T < (float) 0.8*shim->P);
    }

    float Amplitude;

   if(param->Shimmer != -1 ) {
   shim->DeltaShimmer[1] = shim->DeltaShimmer[0];
      do {
	float epsilon =  ((float)rand()) / RAND_MAX;
	
	shim->S =  epsilon  * 4.0 * param->Shimmer - 2.0*param->Shimmer;

	shim->DeltaShimmer[0] = shim->DeltaShimmer[1]*(2.0 + shim->S)/(2.0 - shim->S) +
		      2.0*param->amp*shim->S/(2.0 - shim->S);

	Amplitude = ((float) param->amp + shim->DeltaShimmer[0]);
      } while (  (Amplitude > (float) 1.8*param->amp) || (Amplitude < (float) 0.2*param->amp));
      //	printf("%5.2f \n", S);
	
    }
    else
    {
	Amplitude = (float) param->amp;
    }

    /* generate glottal flow */
    shim->T2 = ceil(0.5*param->cq*shim->P);
    ////
    
    for(i=0; i<shim->T2; i++) {
      shim->x[i] = ceil( Amplitude * 0.5*(1.0 - cos(PI*i/shim->T2) ));
      if(shim->x[i] < param->DC) {
	shim->x[i] = param->DC;
	shim->T4 = i;
      }
    }
    float Knew = param->K * ( 1 + 2 * param->Kvar *  ( ( (1.0 * rand())/RAND_MAX )  - 0.5));
    for(i=shim->T2; i<2*shim->T2; i++) {
      shim->x[i] = ceil( (float) Amplitude*(Knew*cos(PI*(i-shim->T2)/shim->T2) - Knew + 1.0));
      if(shim->x[i] < param->DC) break;
    }

    shim->T3 = i;

    for(i=shim->T3; i<shim->T; i++) {
      shim->x[i] = param->DC;
    }

    if(param->noise != -1) {
      shim->aux = 0.0;
      for(i=shim->T4; i<shim->T3; i++) {
	shim->aux += (float) shim->x[i]*shim->x[i];
      }
      shim->x_pow = shim->aux/((float) shim->T3 - shim->T4);


      shim->aux = 1.0 + ((float) shim->T3 - shim->T4)/((float) shim->T);
      param->NoiseDistWidth = sqrt(12*shim->aux*(shim->x_pow)/param->noise);

      shim->aux = 0.0;
      for(i=0; i<shim->T4; i++) {

	shim->w[i] = (signed short) ceil( ((  1.0 * rand()) / RAND_MAX )* param->NoiseDistWidth - param->NoiseDistWidth/2.0);

	shim->aux += (float) shim->w[i]*shim->w[i];
	shim->x[i] = truncate( (float) shim->x[i] + shim->w[i]);
      }
      for(i=shim->T3; i<shim->T; i++) {

	shim->w[i] = (signed short) ceil( ( ( 1.0*rand() )/RAND_MAX )* param->NoiseDistWidth - param->NoiseDistWidth/2.0);
	shim->aux += (float) shim->w[i]*shim->w[i];
	shim->x[i] = truncate( (float) shim->x[i] + shim->w[i]);
      }
      shim->w_pow = shim->aux/( (float) shim->T);

      //      printf("SNRdb = %5.2f\n", 10.0*log10(x_pow/w_pow));

    }
    return i;
}

int16_t shimmer_get_sample(){
  static uint16_t xx;
  static uint16_t count=0;
  int16_t val;
  float freqq=100.0+(2000.0*_selz);
  shimm.P = shimm.T = (int) (8000.0/freqq);
  if (count==0){
    xx=do_shimmer(&parry, &shimm);
  }
  val=shimm.x[count]; 
  count++;
  if (count>xx) count=0;
  return val;
} 

/////////////////////////////////////

float doklglott(klglott* self){
  //      # puts out the 1st derivative of the glottal flow wave
  float valOut = 0.0f;

  if (self->n < (self->T0 * self->OQ * self->samplingRate)){
  //          # open phase
    float tmp = (float)self->n / (float)(self->samplingRate);
    valOut = 2.0 * self->a * (float)self->n / self->samplingRate - 3.0 * self->b * tmp * tmp;
    //          #valOut = self->a * self->n * self->n - self->b * self->n * self->n * self->n
    valOut -= self->openPhaseCorrection;
    valOut /= self->amplitudeCorrectionFactor;
  }
    self->n += 1;
    if (self->n >= (self->T0 * self->samplingRate))  self->n = 0;  
  return valOut;
}

void setfreq(klglott* self, int freq){
  self->freq = freq;
  self->period = 1.0 / self->freq;
  self->T0 = self->period;// * (float)(self->samplingRate);
    self->OQ = 1.0 - self->CQ;
    self->a = 0.0;
    self->b = 0.0;
	  //        self->calculateParams() ->>>>>>>>>>>>>>>
    self->a = 27.0 * self->amp / (4.0 * self->OQ * self->OQ * self->T0);
    self->b = 27.0 * self->amp / (4.0 * self->OQ * self->OQ * self->OQ * self->T0 * self->T0);

    self->n = 0.0;
      
    self->totalOffset = 0.0;
    self->openPhaseCorrection = 0.0;
    self->amplitudeCorrectionFactor = 1.0;
    int framesPerPeriod = (int)(self->period * (float)(self->samplingRate));
    //int i;
    //for (i=0;i<framesPerPeriod;i++) self->totalOffset += doklglott(self);
    self->openPhaseCorrection = self->totalOffset / (self->OQ * (float)(framesPerPeriod));
    self->n = 0.0;
    self->amplitudeCorrectionFactor = (float)(self->samplingRate) / (float)(freq);
}

void setfreqq(klglott* self, int freq){
  self->freq = freq;
  self->period = 1.0 / self->freq;
  self->T0 = self->period;// * (float)(self->samplingRate);
  self->OQ = 1.0 - self->CQ;
  self->a = 0.0;
  self->b = 0.0;
  //        self->calculateParams() ->>>>>>>>>>>>>>>
  self->a = 27.0 * self->amp / (4.0 * self->OQ * self->OQ * self->T0);
  self->b = 27.0 * self->amp / (4.0 * self->OQ * self->OQ * self->OQ * self->T0 * self->T0);

  //  self->n = 0.0;
      
  //  self->totalOffset = 0.0;
  //  self->openPhaseCorrection = 0.0;
  //  self->amplitudeCorrectionFactor = 1.0;
}


void klglott88_init(klglott* self, int samplingRate, int freq, float CQ, float amp){
  self->samplingRate=samplingRate;
  self->amp=amp;
  self->CQ=CQ;
  setfreq(self, freq);
} 

int16_t klg_get_sample(){
  int16_t s16, freq;
  freq=100+(_selz*1000);
  setfreqq(&glotty, freq);
  float y=doklglott(&glotty);
  s16=y*32768.0;
  return s16*8;
}


///lF to use or not

int16_t pulse_lf(double T0, double Te, double alpha, double omega, double eps){
  //    Assumes E_e = -1.
  static int i=0; int s16; double pulse;
  int n=32000*T0;
  if (i>=0 && i<(int)(Te/T0*n+0.5)) {
    //  for (i=0;i<(int)(Te/T0*n+0.5);i++){
    double t = (double)(i)/n*T0;
    pulse = -(exp(alpha*t)  * sin(omega*t)     /       exp(alpha*Te) / sin(omega*Te));
    s16=pulse*32768.0;
  }
  
  if (i>=(int)(Te/T0*n+0.5) && i<n){
//  for (i=(int)(Te/T0*howmany+0.5);i<n;i++){
    double t = (double)(i)/n*T0;
    pulse = -((exp(-(t-Te)*eps) - exp(-(T0-Te)*eps)) /            (1.0 - exp(-(T0-Te)*eps)));
    //    printf("%f\n", pulse[i]);
    s16=pulse*32768.0;
    //   fwrite(&s16,2,1,fo);
  }
  i++;
  if (i>n) i=0;
return s16;
}

double lf_alpha(double tp, double te, double epsilon, double T0){
/*    """
    Given the three timing parameters and a starting point, uses
    Newton-Raphson to find a value of alpha for the Liljencrants-Fant
    glottal pulse shape.
    """*/
  double tc = T0; double alpha=0.0f;
  double omega = PI / tp; int i;
  for (i=0;i<5;i++){
    double expy = exp(-epsilon*(tc-te));
    double esin = exp(alpha*te)*sin(omega*te);
        double f = (            alpha            - omega/tan(omega*te)            + omega/esin            - (alpha*alpha + omega*omega)*((tc-te)*expy/(1.0-expy) - 1.0/epsilon)	     );
        double fd = (            1.0            - omega*te/esin            - 2.0*alpha*((tc-te)*expy/(1.0-expy) - 1.0/epsilon)            );
        alpha -= f/fd;
  }
	return alpha;
}

double lf_epsilon(double te, double ta, double T0){
  /*    """
    Given Ta, uses Newton-Raphson to find epsilon in an LF model.
    """*/
  double tce = T0 - te; int i;
  double epsilon = 1.0/ta;
  for (i=0;i<5;i++){
    double f = 1.0 - exp(-epsilon * tce) - ta * epsilon;    
    double fd = tce * exp(-epsilon * tce) - ta;
    epsilon -= f/fd;
  }
  return epsilon;
      }

int16_t LF2_get_sample(){
  int16_t s16; 
  double Fa= 400.0*_selz;// ???    
  double Rg= 1.2;//
  double Rk= 0.9; //
  double T0=1.0/Fa; // fundamental period in seconds??? can't be less than 0.001 or??? - 500 Hz - is same as Fa or????

  double     Ta = 1.0/(2.0*PI*Fa);
  double             Tp = T0/(2.0*Rg);
  double       Te = Tp*(Rk+1.0);
  double eps = lf_epsilon(Te, Ta, T0);
  double             alpha = lf_alpha(Tp, Te, eps, T0);
  double             omega = PI / Tp;

  s16=pulse_lf(T0,Te, alpha, omega, eps);
  //  printf("alpha %f omega %f eps %f\n",alpha,omega, eps);
  //    for (z=0;z<1000;z++){
      //    pulse_lf(32000*T0, pulse, T0,Te, alpha, omega, eps);
  //    }
  return s16;
}

/*
void main(){
  int z;
  fo = fopen("testlfgen.pcm", "wb");
  double pulse[32000];
  int puls[32000];
  double Fa= 50;// ???    
  double Rg= 1.2;//
  double Rk= 0.9; //
  double T0=0.02; // fundamental period in seconds??? can't be less than 0.001 or??? - 500 Hz - is same as Fa or????

  double     Ta = 1.0/(2.0*PI*Fa);
  double             Tp = T0/(2.0*Rg);
  double       Te = Tp*(Rk+1.0);
  double eps = lf_epsilon(Te, Ta, T0);
  double             alpha = lf_alpha(Tp, Te, eps, T0);
  double             omega = PI / Tp;

  //  printf("alpha %f omega %f eps %f\n",alpha,omega, eps);
    for (z=0;z<1000;z++){
    pulse_lf(32000*T0, pulse, T0,Te, alpha, omega, eps);
    }
  
  
  float recentY;
  int s16;
  klglott glotty;
  klglott88_init(&glotty, 32000, 100, 0.15, 1.0);

    for (z=0;z<32000;z++){
  float val=doklglott(&glotty);
  float y= recentY + val;
  recentY = y;
  //  printf("%f\n", y);
  s16=y*32768.0;
  fwrite(&s16,2,1,fo);
  }
  shimmer_init(&parry);
  while(1){
  int xx=do_shimmer(puls,&parry,32000);
  for (z=0;z<xx;z++){
    printf("%c",puls[z]>>8);
  }
}  
  }

        elif ptype == 'lf':
            # These three are all in seconds
            Ta = 1.0/(2.0*np.pi*Fa)
            Tp = T0/(2.0*Rg)
            Te = Tp*(Rk+1.0)
            eps = lf_epsilon(Te, Ta, T0)
            alpha = lf_alpha(Tp, Te, eps, T0)
            omega = np.pi / Tp
            pulse_lf(pulse, T0, Te, alpha, omega, eps)
*/
