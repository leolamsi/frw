#include "MathAux.h"

#ifdef _WIN32
  #define MKL 0
#endif 

#if MKL
#include <mkl_vsl.h>
//#include <omp.h>
const int NRndStreams=16;
VSLStreamStatePtr RndStream[NRndStreams];

#endif



//#include <timeval.h>
#include <time.h>

//#include <sys/time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//#include "win32.h"

/* Random Numbers ***************************************************/
int ran2_seed=0; /* seed for ran2() */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-20
#define RNMX (1.0-EPS)

double Ran2()
{
  static long ran2_idum2=123456789;
  static long ran2_iy=0;
  static long ran2_iv[NTAB];
  static long ran2_idum=-12345;

	int j;
	long k;
	double temp;

	if (ran2_idum <= 0) {
	        if (ran2_seed!=0) ran2_idum=-ran2_seed;
		if (-(ran2_idum) < 1) ran2_idum=1;
		else ran2_idum = -(ran2_idum);
		ran2_idum2=(ran2_idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(ran2_idum)/IQ1;
			ran2_idum=IA1*(ran2_idum-k*IQ1)-k*IR1;
			if (ran2_idum < 0) ran2_idum += IM1;
			if (j < NTAB) ran2_iv[j] = ran2_idum;
		}
		ran2_iy=ran2_iv[0];
	}
	k=(ran2_idum)/IQ1;
	ran2_idum=IA1*(ran2_idum-k*IQ1)-k*IR1;
	if (ran2_idum < 0) ran2_idum += IM1;
	k=ran2_idum2/IQ2;
	ran2_idum2=IA2*(ran2_idum2-k*IQ2)-k*IR2;
	if (ran2_idum2 < 0) ran2_idum2 += IM2;
	j=ran2_iy/NDIV;
	ran2_iy=ran2_iv[j]-ran2_idum2;
	ran2_iv[j] = ran2_idum;
	if (ran2_iy < 1) ran2_iy += IMM1;
	if ((temp=AM*ran2_iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/****************************************************/

double Rand1(){
#if _WIN32
  return Ran2();
#else
  return drand48();
#endif
}

double Rand62() 
/* in [0,1), 62 bits from 4 values from ran2 */
{
  long high,low;
  high=(((unsigned int)(Rand1()*32768))<<16) |(((unsigned int)(Rand1()*65536)));
  low =(((unsigned int)(Rand1()*32768))<<16) |(((unsigned int)(Rand1()*65536)));
  return(  (((double)high)+((double)low)/2147483648.) / 2147483648. );
}

int Randbit()
{
  if (Rand1()<0.5) 
    return(0);
  else
    return(1);
}

double Gaudev()
{
	static int iset=0;
	static double gset;
	double fac,r,v1,v2;
	if  (iset==0) {
		do {
		        v1=Rand1()*2.-1;
		        v2=Rand1()*2.-1;
			r=v1*v1+v2*v2;
		} while (r >= 1.0 || r == 0.0);
		fac=sqrt(-2.0*log(r)/r);
		gset=v1*fac;
		iset=1; 
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

void RandArr(double *arrptr, int size, int stream)
// default to stream 0
{
  #if MKL
    //vsRngUniform(VSL_METHOD_SUNIFORM_STD, RndStream, size, arrptr, 0., 1.); //single precision
    vdRngUniform(VSL_METHOD_SUNIFORM_STD, RndStream[stream], size, arrptr, 0., 1.);
    //status = vsrnguniform( method, stream, n, r, a, b )

  #else
    for (int i=0; i<size; ++i) arrptr[i]=Rand1();
  #endif
}


void Randomize(int seed)
/* give seed to random from time of day if seed=0*/
{
   if (seed==0){
	   /*
     struct timeval  t;
     struct timezone tz;
     gettimeofday(&t, &tz);
     seed = t.tv_sec + t.tv_usec;
	 */
	   seed = (unsigned)time(NULL);
   }
   srand((unsigned int) seed);
   ran2_seed=seed;
#if (!_WIN32)
   srand48(seed);
#endif
#if MKL
   for (int s=0; s<NRndStreams; ++s)
     vslNewStream( &RndStream[s], VSL_BRNG_MT2203, (int)(1e6*Rand1()) ); 
   // status = vslNewStream( &stream, brng, seed );
   // A set of 1024 Mersenne Twister pseudorandom number generators.
#endif
}

//uint64_t CPUTimer::tmin=GetCPUTick();
//uint64_t CPUTimer::tmin=0;
