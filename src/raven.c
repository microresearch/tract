
// from raven_lap.c - compare to raven.c

/*
Mindlin
Gardner
*/

#include <stdio.h>
#include "math.h"
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/times.h>
#include <sys/unistd.h>

extern float _selx, _sely, _selz;

double integrate(double old_value, double new_value, double period);
double old_integrate(double old_value, double new_value, double period);

#define PI            3.1415927

//////////////////////// MINDLIN - from Laje, Gardner and Mindlin:

typedef struct mindlin {
double x, xprime,xdobleprime, k, b, c, f0, T, p0;
}Mindlin;

Mindlin syrinxM;

double integrate(double old_value, double new_value, double period);
//  init_mindlin(&syrinxM, 600, 40000, 1200); // guessing these... /raven 1000 12000 1000 0 > testplot // oscillates
void init_mindlin(Mindlin* mind, double b, double k, double c){
mind->b=b;
mind->k=k;
mind->c=c;
mind->x=0.01;
 mind->p0=100000.0f; // unused
mind->xdobleprime=0.0f;
mind->f0=2000.0f;
mind->T = 1.0/48000.0;

}

// paper/book p51
/*

xdouble= -kx + ((B-b)xprime) - (c*x^2*xprime) (no f0 there)

B is bounding motion = touching the walls

k-springiness
b-friction
c-non-linear dissipation
p0?
f0?
 */
double calc_xdobleprime_mindlin(Mindlin *mind){
// from pseudocode: result/0?/=xdobleprime+k*xprime + (c*x)^2 * xprime * xprime - b * xprime + f0;
// xdobleprime = - (k*x + c *x^2 * xprime * xprime - b * xprime + f0);

//return - (mind->k*mind->x + mind->c* pow(mind->x,2) * mind->xprime * mind->xprime - mind->b * mind->xprime + mind->f0);

// from paper: xdobleprime= -k*x-c*x^2*xprime+b*xprime-f0;

// other paper xxdobleprime= -k*x-(B1 + B2 * x^2 - p0/=air sac pressure/) * xprime - f0

return -mind->k*mind->x - mind->c*pow(mind->x,2) * mind->xprime+mind->b *mind->xprime - mind->f0; // exp overruns

//  return -( mind->k*mind->x + mind->c*pow(mind->x,2) * mind->xprime - mind->b *mind->xprime + mind->f0); // exp overruns - same as above which is same as equation
// from: http://www.scholarpedia.org/article/Models_of_birdsong_%28physics%29

//return -mind->k * mind->x - (mind->b + pow(mind->c + mind->x,2) - mind->p0)* mind->xprime - mind->f0;

}

void no_newsay(){
}

double integrate(double old_value, double new_value, double period){
  return (old_value + new_value)*(period/2.0);
}


int16_t mindlin_oscillate(Mindlin* mind){
int16_t iii;  
// mind->f0=_selz*20000.0f;
//   init_mindlin(&syrinxM, 600, 40000, 1200); // guessing these... /raven 1000 12000 1000 0 > testplot // oscillates b k c
 mind->c=200.0+(_selz*2000.0f);
double newxdobleprime = calc_xdobleprime_mindlin(mind);
double newxprime = mind->xprime + integrate(mind->xdobleprime,newxdobleprime,mind->T);
 mind->xdobleprime = newxdobleprime;
 mind->x += integrate(mind->xprime,newxprime,mind->T);
//printf("%f\n", mind->x);
 iii=32768.0f*mind->x;
//fwrite(&iii,2,1,fo);
 mind->xprime = newxprime;
 return iii;//mind->x;
}

int16_t mindlin_get_sample(){
  int16_t i;
  i=mindlin_oscillate(&syrinxM);
  return i;
}

////////////////////// GARDNER

typedef struct gardner {
  double x, xprime,oldxprime,xdobleprime, K, Pb, a0, b0, t, M, K_scale, K_scalex, Pb_scalex, D, D2, Pb_scale, T, freq, ofreq;
}Gardner;

Gardner syrinx;

void init_gardner(Gardner* gd, int16_t kk, int16_t pb){
gd->x=0;
gd->xprime=1.0;
gd->xdobleprime = 0.0f;
gd->K =  (double)kk; // between 11k and 12k something happens... - and now???
gd->Pb = (double)pb;

/*
	upper_labia = 0.02 #cm
	lower_labia = 0.04 #cm
	t_constant = .00015 #s # phenomenological constant Titze!
	mass = .005 #g/cm3
	K=restitution_constant = 200000 #g*cm/s2cm3 // this is K- IF= approx 200 kdyn/cm3 ??? // paper 0-8 N/cm3

8 Newtons is 800 kdyn = 800 000 dyne = 

	D_coefficient = 5 #dynes*s/cm3
	D2_coefficient = .01 #dyne*s/cm5 .001
	bronchial_pressure = 10000 #g/(cm*s2) # 0-3 kPa

1 pascal is 1 N/m2 = 1 kg /ms2 =  say 2 kPa = 2 N/m2 2000 kPa=20000 g/s
*/

//#constants
 gd->a0 = 0.02;
 gd->b0 = 0.08;
 gd->t = 0.00015;
 gd->M = 0.005;
 gd->K_scale = 1.0;
 gd->K_scalex = 1000.0;
 gd->D = 5.0;
 gd->D2 = 0.001;
 gd->Pb_scale = 1.0;
 gd->Pb_scalex = 1000.0;
 gd->T = 1/48000.0;
 gd->freq = 100.0;
 gd->ofreq= 200.0;
}

double calc_xdobleprime_gardner(Gardner *gd){
  double a = gd->a0 + gd->x + (gd->t*gd->xprime);
  double b = gd->b0 + gd->x - (gd->t*gd->xprime);
  double Pf = gd->Pb*(1 - (a/b));
  //  printf("%f\n",gd->x);
  //  int iii=(double)gd->x*327680.0f;
  //  fwrite(&iii,2,1,fo);
  //		return (Pf - (self.K*self.x) - (self.D2*math.pow(self.xprime,3)) - (self.D*self.xprime))/self.M
  //    return (Pf - (gd->K*gd->x) - (gd->D*gd->xprime))/gd->M;
  //  return (Pf - (gd->K*gd->x) - (gd->D2*pow(gd->oldxprime,3)) - (gd->D*gd->xprime))/gd->M;
return ((gd->Pb*((gd->a0 - gd->b0) +(2*gd->t*gd->xprime)/(gd->x+gd->b0+(gd->t*gd->xprime))))-((gd->K*gd->x) + (gd->D2*pow(gd->xprime,3)) + (gd->D*gd->xprime)))/gd->M;
}


int16_t gardner_oscillate(Gardner* gd){
double newxdobleprime = calc_xdobleprime_gardner(gd);
double newxprime = gd->xprime + integrate(gd->xdobleprime,newxdobleprime,gd->T);
 static int16_t i;
 i++;
 if (i>32000) i=0;
 gd->K=_selz*4000.0;
 // gd->Pb=_selz*2000.0;
gd->xdobleprime = newxdobleprime;
gd->x += integrate(gd->xprime,newxprime,gd->T);
//fwrite(&gd->x,2,1,fo);
gd->oldxprime=gd->xprime;
gd->xprime = newxprime;
//gd->K+= gd->K_scale;
// gd->Pb+= gd->Pb_scale;
// gd->K = gd->K_scale*sin((2*PI*gd->T*gd->freq*i)) + gd->K_scalex;
 //		gd->K = gd->K_scale*2.0
// gd->Pb = gd->Pb_scale*cos((2*PI*gd->T*gd->ofreq*i) + PI) + gd->Pb_scalex;
		  //               gd->Pb = gd->Pb_scale
  //   gd->freq += 0.01;
  //   gd->ofreq+=0.01;
 return (gd->x*3276800.0);
}

int16_t gardner_get_sample(){
  int16_t i;
  i=gardner_oscillate(&syrinx);
  return i;
}

 /*
  *  based on balloon1.cpp - IF model...

 See http://www.dei.unipd.it/~avanzini/downloads/paper/avanzini_eurosp01_revised.pdf - measurements

 *  Implements Ishizaka & Flanagan two mass model of the vocal folds
 *	Uses backward finite difference approximation as well as many others
 *	Output is the flow (u)


  */

 double computeSample(double pressure_in);
 void clearOld();

 double ps;		//subglottal pressure
 double r1;		//damping factor
 double r2;
 double m1;		//mass
 double m2;
 double k1;		//spring constant
 double k2;
 double k12;		//coupling spring constant
 double d1;		//glottal width
 double d2;
 double lg;		//glottal length
 double aida;		//nonlinearity coefficient
 double S;			//subglottal surface area
 double Ag01;		//nominal glottal area, with mass at rest position
 double Ag02;
 double pm1Prev;	//pressure at previous time step
 double pm2Prev;
 double x1Prev;	//displacement at previous time step
 double x1PrevPrev;//displacement at previous time step to the previous one
 double x2Prev;
 double x2PrevPrev;
 double gain;		//after-market gain
 double uPrev;		//previous flow value
 double Fs;		//calculation sampling rate, not actual audio output sample rate


void init_balloon(){

   // adapt these settings for potential raven voice

   // from MATLAB code: also there is pressure envelope there

   /*

 p=0;        %relative output pressure, 
 rho = 1.14; %kg/m^3 mass density
 v = 1.85e-5; %N*s/m^2 greek new: air shear viscosity
 lg = 1.63e-2; %m glottal length

 twod = 3e-5; %m, glottal width 2d1
 d1=twod/2; %1.5000e-005
 d2=d1;
 m = 4.4e-5/90; %kg, glottal mass // why /90 - otherwise accords for human glottis
 m1=m; %4.8889e-007
 m2=m;
 k12=0.04; %coupling spring constant
 k = 0.09; %N/m, spring constant
 k1=k;
 k2=k1;
 aida=1000000.01; %non-linearity factor
 r = 0.0001*sqrt(m*k); %damper constant, N*s/m
 r1=r*1;
 r2=r1; %2.0976e-008
 %Ag0 = 5e-6; %m^2 glottal area at rest = 5mm^2=5e-6
 Ag0 = 5e-9; %m^2 glottal area at rest
 S = 5e-5; %m^2 output area (vocal tract end)
 %S = 5e-4; %m^2 output area (vocal tract end)

   */

   /* raven details (see kahrs.pdf and zacarelli)

 kahrs: 

 from zacarelli we have:

 stiffness (g ms−2)	k1, k2	22.0×10−3
 damping constant (g ms−1)	r1, r2	1.2×10−3
 coupling constant (g ms−2)	kc	6.0×10−3

 but not sure how to convert between????

 m1/m2=glottal mass - 3.848451000647498e-6 - 0.00000384 /90

 k1/k2-spring constant N/m - 3.11 ???
 k12=coupling spring constant ???

 d1/d2=glottal width 2dl /2??? diameter is 7mm  say 2mm now or is this *thickness?* 1e-4 - from fletcher is 100 micrometer
 r1/r2 = 0.0001*sqrt(m*k); %damper constant, N*s/m - 1.386e-7 // but depends on K spring constant can vary 0.0001

 Ag0 = 5e-9; %m^2 glottal area at rest - 2mm say at rest= 3.14mm 3.14e-6
 S = 5e-5; %m^2 output area (vocal tract end) - 20mm diameter BEAK 314mm = 0.000314

 lg= 1.63e-2; %m glottal length - say 7mm=7e-3

    */
d1 =0.0008;
d2 =0.0008;

r1=0.00000980665; // damper constant depends on k and mass r = 0.0001*sqrt(m*k); %damper constant, N*s/m Avanzini has 0.1 * sqrt(m*k)  m0.000044 * 20 = 0.002966
// above 0.001 is no sound
	r2 =r1;
m1 =1e-6; // -7 or -5 for -5 we would have for r1 and as k1=0.09 = 
	m2 =m1;
	k1 =0.014;
	k2 =k1;
	k12=0.00005;
	aida =10000000.0;
//	aida =10.0;
	d1 =1.5e-2;
	d2 =d1;

	lg =0.007;
	gain=400.0;
	S=0.0005;
	Ag01=3e-7;
	Ag02=3e-7; 

	x1Prev=0.0;
	x1PrevPrev=0.0;
	x2Prev=0.0;
	x2PrevPrev=0.0;
	pm1Prev=0.0;
	pm2Prev=0.0;
	uPrev=0.0;
	ps=0.0;
	Fs=32000.0;
  
}


int16_t balloon_get_sample() {

/*	bsynth->setM1(5e-8*(double)[(id)dataPointer m1In]);
	bsynth->setM2(5e-8*(double)[(id)dataPointer m2In]);
	bsynth->setR1(5e-9*(double)[(id)dataPointer r1In]);
	bsynth->setR2(5e-9*(double)[(id)dataPointer r2In]);
	bsynth->setK1(1e-3*(double)[(id)dataPointer k1In]);
	bsynth->setK2(1e-3*(double)[(id)dataPointer k2In]);
	bsynth->setD1(1.5e-7*(double)[(id)dataPointer d1In]);
	bsynth->setD2(1.5e-7*(double)[(id)dataPointer d2In]);
	bsynth->setK12(1e-4*(double)[(id)dataPointer k12In]);
	bsynth->setLg(1.3e-4*(double)[(id)dataPointer lgIn]);
	bsynth->setAida(1.1e-4*(double)[(id)dataPointer aidaIn]);
	bsynth->setS(5.5e-7*(double)[(id)dataPointer SIn]);
	bsynth->setAg01(5.1e-10*(double)[(id)dataPointer Ag01In]);
	bsynth->setAg02(5.1e-10*(double)[(id)dataPointer Ag02In]);
	bsynth->setGain((double)[(id)dataPointer gainIn]);
	bsynth->setFs((double)[(id)dataPointer FsIn]);
*/

/* pressure waveform:

%Set pressure waveform envelope
for n=1:N
if n<5
%ps(n)=MAXps;
ps(n)=0;
elseif n<7
ps(n)=0;
%ps(n)=MAXps;
elseif n<T1 // T1 is N/90 - N is number of iterations/samples
ps(n)=ps(n-1)+MAXps/T1; // MAXps is 400
elseif n<=T2 // T2 is N/80 
    ps(n)=ps(n-1);
else
    ps(n)=ps(n-1)+(MINps-MAXps)/(N-T2); // MINps is 30
*/

  static int16_t i; static double lastp;
  double pressureIn=_selz*1000.0;
  //	for (i=0; i<bufferSize; i++ ) {
	  //	  *samples++ = computeSample(pressureIn);
	  //	  printf("%f  ",computeSample(pressureIn));
  i++;
  //  if (i<5) pressureIn=0;
  //  else if (i< (32000/90)) pressureIn=lastp+ 400/(32000/90);
  //  else if (i<= (32000/80)) pressureIn=lastp;
  //  else pressureIn=lastp-360.0f/(32000.0f-(32000.0/80.0));
  if (i>32000) i=0;
  
  // constant pressure
//	  pressureIn=300; // 0.3 kPa after Fletcher

  signed int s16=(signed int)(computeSample(pressureIn)*32768.0);
  //	  *samples++=(double)(computeSample(pressureIn));
  //      printf("%d\n",s16);
  //      fwrite(&s16,2,1,fo);
  lastp=pressureIn;
  return s16;	
};

double computeSample(double pressure_in){
  double T=1/Fs;
  double rho = 1.14; 
  double rhosn = rho*0.69;
  double hfrho=rho/2;
  double v = 1.85e-5;
  double twvd1lg=12*v*d1*lg*lg;
  double twvd2lg=12*v*d2*lg*lg;
  double Ag012lg=Ag01/2/lg;
  double Ag022lg=Ag02/2/lg;
  double lgd1=lg*d1;
  double lgd2=lg*d2;
  double m1T=m1/T/T;
  double m2T=m2/T/T;
  double r1T=r1/T;
  double r2T=r2/T;
  double C11=k1*(1+aida*x1Prev*x1Prev);
  double C12=k2*(1+aida*x2Prev*x2Prev);
  double C21=k1*(1+aida*(x1Prev+Ag012lg)*(x1Prev+Ag012lg));
  double C22=k2*(1+aida*(x2Prev+Ag022lg)*(x2Prev+Ag022lg));
  double alpha1=lgd1*pm1Prev;
  double alpha2=lgd2*pm2Prev;
  double beta1=m1T*(x1PrevPrev-2*x1Prev);
  double beta2=m2T*(x2PrevPrev-2*x2Prev);
  double gamma1=-r1T*x1Prev;
  double gamma2=-r2T*x2Prev;
  double delta1=Ag012lg*C21;
  double delta2=Ag022lg*C22;
  double lambda1=-k12*x2Prev;
  double lambda2=-k12*x1Prev;
  double x1=0.0;
  double x2=0.0;
  double pm1=0.0;
  double pm2=0.0;
  double A1=0.0;
  double A2=0.0;
  double A1n2=0.0;
  double A1n3=0.0;
  double A2n2=0.0;
  double A2n3=0.0;
  double a=0.0;
  double b=0.0;
  double c=0.0;
  double det=0.0;
  double flow1=0.0;
  double flow2=0.0;
  double udif1=0.0;
  double udif2=0.0;
  double u=0.0;
  double g1=0.0;
  double g2=0.0;
  double g4=0.0;
  double g5=0.0;
  double pm1b=0.0;
  double pm2b=0.0;


  if  (x1Prev>=-Ag012lg){
    x1=(alpha1-beta1-gamma1-lambda1)/(m1T+r1T+C11+k12);
  }
  else {
    x1=(alpha1-beta1-gamma1-lambda1-delta1)/(m1T+r1T+C21+k12);
  }

  if  (x2Prev>=-Ag022lg){
    x2=(alpha2-beta2-gamma2-lambda2)/(m2T+r2T+C12+k12);
  }
  else{
    x2=(alpha2-beta2-gamma2-lambda2-delta2)/(m2T+r2T+C22+k12);
  }

  A1=Ag01+lg*x1;
  A2=Ag02+lg*x2;
    
  if (A1<=0)
    {A1=0.1e-25;}
    

  if (A2<=0)
    {A2=0.1e-25;}
   
	
	
  A1n2=A1*A1; 
  A1n3=A1n2*A1;

  A2n2=A2*A2; 
  A2n3=A2n2*A2;
	
  a= (rhosn/A1n2)+hfrho*(1/A2n2-1/A1n2)+hfrho/A2n2*(2*A2/S*(1-A2/S));
  b= twvd1lg/A1n3+twvd2lg/A2n3;
  c= -pressure_in;

  det=b*b-4.0*a*c;

  if (det>=0){
    flow1=(-b+sqrt(det))/(2.0*a);
    flow2=(-b-sqrt(det))/(2.0*a);
  }
  else{
    flow1=(-b)/(2.0*a);
    flow2=(-b)/(2.0*a);
  }

  udif1=fabsf(flow1-uPrev);
  udif2=fabsf(flow2-uPrev);

  if (udif1<udif2){
    u=flow1;
  }
  else{
    u=flow2;
  }

  //u=max(flow1,flow2);



  g1=rhosn*u*u/A1n2;
  g2=twvd1lg*u/A1n3;
  g4=twvd2lg*u/A2n3;
  g5=hfrho*u*u/A2n2*(2*A2/S*(1-A2/S));

  pm1=pressure_in-g1-g2/2;
  pm2=g5+g4/2;

  if (x1>=-Ag012lg){
    
    pm1b=(m1T*(x1-2*x1Prev+x1PrevPrev)+r1T*(x1-x1Prev)+ k1*x1*(1+aida*x1*x1) +k12*(x1-x2))/(lgd1);
    pm1=pm1/2+pm1b/2;
        
	
    if (x2>=-Ag022lg){

      pm2b=(m2T*(x2-2*x2Prev+x2PrevPrev)+r2T*(x2-x2Prev)+ k2*x2*(1+aida*x2*x2) -k12*(x1-x2))/(lgd2);
      pm2=pm2/2+pm2b/2;
    }
    else{
                 pm2=pressure_in;
               pm1=pressure_in;
    }
  }     
  else{
       pm1=pressure_in;
       pm2=0.0;
  }

    if (pm1>pressure_in)
    	{pm1=pressure_in;}
    
    if (pm2>pressure_in)
           {pm2=pressure_in;}
    
      if (pm1<0.0)
	//	   pm1(n)=abs(pm1(n));
	{pm1=fabsf(pm1);}
 

  if (pm2<0.0)
    {pm2=0.0;}
    
  //  if (u<0.0)
  //    {u=0.0;}
		
    pm1Prev=pm1;
  pm2Prev=pm2;
  x1PrevPrev=x1Prev;
  x1Prev=x1;
  x2PrevPrev=x2Prev;
  x2Prev=x2;
  uPrev=u;
  ps=pressure_in;
  //printf("UUUU %f pressure-in %f\n", u, pressure_in);
  return gain*u;
}

void clearOld()
{
x1Prev=0;
x1PrevPrev=0;
x2Prev=0;
x2PrevPrev=0;
pm1Prev=0;
pm2Prev=0;
uPrev=0;
ps=0;
}
