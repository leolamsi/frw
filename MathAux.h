#include <assert.h>
#include <iostream>
#include <fstream>
using namespace std;

void Randomize(int seed=0);
double Rand1();
double Ran2();
double Rand62();
int Randbit();
double Gaudev();
void RandArr(double *arrptr, int size, int stream=0);
//unsigned long long GetCPUTick();

/* Timer ***************************************************/

#if defined(_MSC_VER) // Or any intel style ASM code
	typedef unsigned __int32 uint32_t;
	typedef unsigned __int64 uint64_t;
#	define RDTSC(low, high)	\
	__asm rdtsc				\
	__asm mov low, eax		\
	__asm mov high, edx
#elif defined(__GNUC__)
#	include <inttypes.h>
#	define RDTSC(low, high)	\
	__asm__ __volatile__("rdtsc" : "=a" (low), "=d" (high))
#else
#	error "Not supported platform"
#endif

// Returns the current CPU clock tick. You have to divide the clock cycle by the CPU frequency in order to
// convert the value into seconds.
// But it's already enough to measure the percentage time spend of a function over the whole program
static uint64_t GetCPUTick() {
	union Counter64 {
		uint32_t n32[2];
		uint64_t n64;
	};

	uint32_t l, h;
	RDTSC(l, h);
	Counter64 tmp;
	tmp.n32[0] = l;
	tmp.n32[1] = h;
	return tmp.n64;
}

#if 0
unsigned long long GetCPUTick() {
	union Counter64 {
		uint32_t n32[2];
		uint64_t n64;
	};
	static Counter64 cnt;
	uint64_t cnt_old = cnt.n64;
	
	uint32_t l, h;
	RDTSC(l, h);
	cnt.n32[0] = l;
	cnt.n32[1] = h;
	//return cnt.n64;

	return (unsigned long long) (cnt.n64-cnt_old);
}
#endif

class CPUTimer
{
  uint64_t tsum, dt, t_start, ncall;  
  bool timing;
  static uint64_t tmin;
public:

  CPUTimer(){
    tsum=dt=t_start=ncall=0;
    timing=false;
  };

  void start(int line=0){
    if (timing!=false)
      cout << "WARNING: timing started twice (line: " << line << " )\n"; 
    t_start = GetCPUTick();
    timing=true;
    if (tmin==0) tmin = t_start; // initialize tmin as late as now to skip target program initializations
  }

  void stop(int line=0){
    if (timing!=true) 
      cout << "WARNING: timing stoped before started (line: " << line << " )\n"; 
    dt = GetCPUTick() - t_start;
    tsum += dt;
    ++ncall;
    timing=false;
  }

  double pctg(){
    uint64_t t=GetCPUTick();
    return 100. * tsum / ( t-tmin);
  }

  friend ostream& operator<< (ostream&os, CPUTimer T){
    double dtmean;
    if (T.ncall==0) dtmean=0; else dtmean=(double)T.tsum/T.ncall;
    return os << T.pctg() << "% " << dtmean/1000. << " kcycs "
	      << T.dt << " cycs " << T.ncall << " calls ";
  }
};

#undef RDTSC
