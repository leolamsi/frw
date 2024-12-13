/*
syntax:  [exe] [L] [N] [name of input file] [name of output file]
*/

#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <iomanip>
#include <sstream>
using namespace std;

int L;
char name[100];
char newname[100];
ofstream outfile;
int N, dt;

int main(int argc, char** argv){
  
  if (argc >= 2) {
    istringstream a(argv[1]);
    a >> L;
    istringstream aa(argv[2]);
    aa >> N;
    istringstream aaa(argv[3]);
    aaa >> dt;
    string str_name = argv[4];
    for (int i = 0; i < str_name.length(); ++i)
      name[i] = str_name[i];
    string str_new_name = argv[5];
    for (int i = 0; i < str_new_name.length(); ++i)
      newname[i] = str_new_name[i];

  }
  ifstream infile(name);
  if (!infile) {
    cout << "Failed to open file" << endl;
    exit(0);
  }
  outfile.open(newname);
  
  int d = 0, L_count = -L/2;
  long long int sum = 0, sqsum = 0;
  double n_total = N*10000;
  int ndt = 0;
  while (infile >> d){
    //    cout << L_count << ":" << d << " ";
    if (d != 0){
      sum += (L_count * d);
      sqsum += (L_count * L_count * d);
    }
    ++L_count;
    if (L_count == L/2){
      double mean = (double)sum/n_total;
      double variance = sqsum/n_total - mean*mean;
      //cout << sum << ", " << sqsum << endl;
      //cout << mean << ", " << variance << endl;
      double msd = sqrt(variance);
      outfile << ndt*dt << " " <<  msd << endl;
      ++ndt;
      sum = 0; sqsum = 0; L_count = -L/2;
      //cout << endl;
    }
  }
  infile.close();
  outfile.close();
}
