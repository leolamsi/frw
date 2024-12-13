/*
syntax: frw.exe [q, N, L, t_max, nfr, dir name, file name]
*/
#include <iostream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include "FileIO.h"
#include "MathAux.h"
#define GRAPHICS 0
#define pair_init_condition 0
#define pair_deviation 1
#define single_init_condition 0
#if __CYGWIN__
#include <GL/glut.h>
#undef GRAPHICS
#define GRAPHICS 0
#endif
using namespace std;
using namespace FileIO;

#if !pair_deviation
const int max_tree_lv = 20; //maximum lattice size = 2^20 = 1048576
#else
const int max_tree_lv = 14;
#endif
const int tree_widest = pow(2, max_tree_lv);
int **bin_tree;
int tree_h;

const int L_MAX = tree_widest, N_MAX = tree_widest;

const double pi = 3.14159265;
const double M = pow(2., 32);

int t_max = 10000; //default value
int outn = 20; //default value
double dt;
double t = 0;
int t_cnt = 0;

int N_hop_ok = 0;

char file_name[100]; // output file name

//double q = 1. / sqrt(3.);
double q = 0.6;

#if pair_deviation
int pair_trials = 100;       // number of pairs simulated
int N;
const int t_out_MAX = 1000;
ofstream deviation;
#else
int N = 10; // default number of particles
#endif

#if pair_init_condition || single_init_condition
int separation = 10;
#endif

#if pair_init_condition
int L = (N / 2) * 2 * separation;
#elif pair_deviation
int L = 100;
#elif single_init_condition
int L = N * (separation + 1);
#else
int L = 1000; 
#endif

#if pair_deviation
int histdata[L_MAX][t_out_MAX];
#endif

#if GRAPHICS
#include <glut.h>
GLfloat color[N_MAX][3];

void timer(int e) {
    glutPostRedisplay();
    glutTimerFunc(10, timer, 0);
}

void initGL() {
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glutInitDisplayMode(2.);
}
const int y_height = 8;
int pctl_R[L_MAX][y_height];
const int window_width = 1500, window_height = 500;
#endif


class vec {
public:
    unsigned int x, y;
    vec(unsigned int xi, unsigned int yi) {
	x = xi;
	y = yi;
    }
    vec() {}
    vec operator+(vec b) {
	return vec(x + b.x, y + b.y);
    }
    vec operator-(vec b) {
	return vec(x - b.x, y - b.y);
    }
};

ostream& operator<< (ostream& os, vec A) {
    return os << A.x << " " << A.y << " ";
};

vec Phi[L_MAX];

unsigned int dot(vec a1, vec b1) {
    return a1.x * b1.x + a1.y * b1.y;
}

class matrix22 {
public:
    unsigned int p[2][2];
    matrix22(unsigned int a, unsigned int b, unsigned int c, unsigned int d) {
	p[0][0] = a;
	p[0][1] = b;
	p[1][0] = c;
	p[1][1] = d;
    }
    matrix22() {}

    vec operator*(vec v) {
	return vec(dot(vec(p[0][0], p[0][1]), v), dot(vec(p[1][0], p[1][1]), v));
    }

    matrix22 operator*(matrix22 m) {
	matrix22 g;
	g.p[0][0] = dot(vec(p[0][0], p[0][1]), vec(m.p[0][0], m.p[1][0]));
	g.p[0][1] = dot(vec(p[0][0], p[0][1]), vec(m.p[0][1], m.p[1][1]));
	g.p[1][0] = dot(vec(p[1][0], p[1][1]), vec(m.p[0][0], m.p[1][0]));
	g.p[1][1] = dot(vec(p[1][0], p[1][1]), vec(m.p[0][1], m.p[1][1]));
	return g;
    }

    matrix22 transpose() {
	return matrix22(p[0][0], p[1][0], p[0][1], p[1][1]);
    }
};

ostream& operator<< (ostream& os, matrix22 A) {
    return os << A.p[0][0] << " " << A.p[0][1] << endl << A.p[1][0] << " " << A.p[1][1] << " ";
};

class Ptcl {
public:
    int x;
    int number;
    int y_cood;
    int init_pos;
    vec Psi;

    int theta(double b) {
	return b >= 0;
    }

    unsigned int S(unsigned int n) {
	return (n >> 16) | (n << 16);
	//return n;
    }

    int f(int xi) { // 0 is left, 1 is right
	unsigned int n = dot(Psi, Phi[(x + xi) % L]);
	return theta(q - S(n) / M);
    }

    void hop(int xi) {
	unsigned int n = dot(Psi, Phi[(x + xi) % L]);
	matrix22 U[2]{ matrix22(n + 1, n, 1, 1), matrix22(1, -(int)n, -1, n + 1) };
		
	// change particle state
	Psi = U[1 - xi].transpose() * Psi;

	// change bond state
	Phi[(x + xi) % L] = U[xi] * Phi[(x + xi) % L];
    }
};
Ptcl ptcl[N_MAX];

#if GRAPHICS
void drop(int k) {
    for (int i = ptcl[k].y_cood; i < y_height - 1; ++i) {
	pctl_R[ptcl[k].x][i] = pctl_R[ptcl[k].x][i + 1];
	if (pctl_R[ptcl[k].x][i] != -1)
	    ptcl[pctl_R[ptcl[k].x][i]].y_cood -= 1;
	else
	    break;
    }
    pctl_R[ptcl[k].x][y_height - 1] = -1;
}

void stack(int k) {
    for (int i = 0; i < y_height; ++i) {
	if (pctl_R[ptcl[k].x][i] == -1) {
	    pctl_R[ptcl[k].x][i] = ptcl[k].number;
	    ptcl[k].y_cood = i;
	    break;
	}
    }
}
#endif

void update_position(int xi, int k) {
#if GRAPHICS
    drop(k);	// gravity for top particles
#endif
    ptcl[k].x = (ptcl[k].x - 1 + xi * 2 + L) % L; //change position
#if GRAPHICS
    stack(k);
#endif
}

void update_N_hop_ok() {
    N_hop_ok = bin_tree[max_tree_lv - 1][0] + bin_tree[max_tree_lv - 1][1];
    //if (N_hop_ok == 0)
    //cout << "no hops initially available" << endl;
    //cout << "nhop " << N_hop_ok << endl;
}

void refresh_tree(int x, int level) {
    if (level != max_tree_lv-1) {
	x -= x & 1;
	int next_x = x / 2;
	bin_tree[level + 1][next_x] = bin_tree[level][x] + bin_tree[level][x + 1];
	refresh_tree(next_x, level + 1);
    }
}

int t0_N_hop_ok = 0;
void next_hop() {	// choose particle k to hop
    if (t == 0)  // save initial available hops
	t0_N_hop_ok = N_hop_ok;
    if (t != 0 && t0_N_hop_ok != 0 && N_hop_ok == 0) {  // error case: N_hop_ok suddenly turned to 0
	cout << "Error at line " << __LINE__ << " -- particles are suddenly stuck" << endl;
	exit(0); // break program
    }
    if (N_hop_ok == 0) {  // skip all calculations if no hops available 
	t += dt;
	//cout << "dt: " << dt << "  t: " << t << endl;
    }
    else {
	//cout << "rand ptcl" << endl;
	int ran_ptcl = (int)(Ran2() * N_hop_ok); // randomize ptcl
	t += 1. / N_hop_ok;    // add time
    
	//search tree
	int index = 0;
	for (int i = max_tree_lv - 1; i >= 0; --i) {
	    index *= 2;
	    if (ran_ptcl >= bin_tree[i][index]) {
		ran_ptcl -= bin_tree[i][index];
		index += 1;
	    }
	}
	int xi = index & 1;
	int k = (index - xi) / 2;

	int oldpos = ptcl[k].x;
    
	// perform hop
	ptcl[k].hop(xi);
	update_position(xi, k);
    	
	for (int i = 0; i < N; ++i) {
	    if (ptcl[i].x == ptcl[k].x) {
		bin_tree[0][i * 2 + 1 - xi] = ptcl[i].f(1-xi);
		refresh_tree(i * 2, 0);
	    }
	    if (ptcl[i].x == oldpos) {
		bin_tree[0][i * 2 + xi] = ptcl[i].f(xi);
		refresh_tree(i * 2, 0);
	    }
	}
	bin_tree[0][index] = ptcl[k].f(xi);
	refresh_tree(index, 0);
	//cout << "tree refreshed" << endl;
	update_N_hop_ok();
	//cout << "hopped" << endl;
    }
}

int out_cnt = 0;
int t_inc = 0;  // should be 0 during runs

void out_to_file() {
    //cout << "in out to file: " << t << endl;
    if ((t-t_inc >= (double)dt * t_cnt && t > t_inc) || t == 0) { //output to file 
	t_cnt += 1;
	int p_xpos[N_MAX];
	for (int i = 0; i < N; ++i) {
#if single_init_condition 
	    p_xpos[i] = ptcl[i].x - ptcl[i].init_pos;
#else
	    p_xpos[i] = ptcl[i].x;
#endif
	}
#if !pair_deviation
	AppendArrs(file_name, N, 1, p_xpos);

	/*unsigned int Phix[L], Phiy[L];
	  for (int i = 0; i < L; ++i) {
	  Phix[i] = Phi[i].x;
	  Phiy[i] = Phi[i].y;
	  }
	  AppendArrs("Phi", L, 2, Phix, Phiy);*/ 
	if ((double)t / t_max > out_cnt/100.) {
	    cout << out_cnt << "%" << endl;
	    out_cnt += 10;
	}
	//cout << t << endl;
#else
//	cout << t << "     " << p_xpos[0] << "    " << p_xpos[1] << endl;
	for (int i = 0; i < N; ++i) {
	    histdata[ptcl[i].x][t_cnt - 1] += 1;
	    //cout << "print" << endl;
	}
	
	char fname[100];
	int tt = sprintf(fname, "%s_traj", file_name);
	AppendArrs(fname, N, 1, p_xpos);
//	cout << "t " << t << endl << endl;	
/*
	char fname2[100];
	tt = sprintf(fname2, "%s_psi0", file_name);
	unsigned int psix[1], psiy[1], t_temp[1];
	psix[0] = ptcl[0].Psi.x;
	psiy[0] = ptcl[0].Psi.y;
	t_temp[0] = (int)t;
	AppendArrs(fname2, 1, 3, t_temp, psix, psiy);

	char fname3[100];
	tt = sprintf(fname3, "%s_psi1", file_name);
	psix[0] = ptcl[1].Psi.x;
	psiy[0] = ptcl[1].Psi.y;
	t_temp[0] = (int)t;
	AppendArrs(fname3, 1, 3, t_temp, psix, psiy);
	*/	
#endif
    }
}

void display() {
#if GRAPHICS 
    glClear(GL_COLOR_BUFFER_BIT);
    int z = max_tree_lv;  // no. of circle vertex
    double r = 300 / L;
    for (int i = 0; i < L; ++i) {
	for (int j = 0; j < y_height; ++j) {
	    if (pctl_R[i][j] != -1) {
		glColor3f(color[pctl_R[i][j]][0], color[pctl_R[i][j]][1], color[pctl_R[i][j]][2]);
		glBegin(GL_POLYGON);
		for (int vert = 0; vert < z; ++vert) {
		    glVertex3f(r * 0.5 + r * i + 0.5 * r * cos(vert * 2. * pi / z), r * 0.5 + r * j + 0.5 * r * sin(vert * 2. * pi / z), 0);
		}
		glEnd();
	    }
	}
    }
    glutSwapBuffers();
#endif
    out_to_file();
    next_hop();
#if GRAPHICS
    if (t >= t_max) { //end gl main loop
	exit(0);
    }
#endif
}

void bin_tree_init() {
    for (int i = 0; i < max_tree_lv-1; ++i) {
	for (int j = 0; j < pow(2, max_tree_lv - i - 1); ++j) {
	    bin_tree[i + 1][j] = bin_tree[i][2 * j] + bin_tree[i][2 * j + 1];
	}
    }
    update_N_hop_ok();
}

void main_init() {
#if GRAPHICS
    for (int i = 0; i < L; ++i) {
	for (int j = 0; j < y_height; ++j)
	    pctl_R[i][j] = -1;
    }
#endif
#if pair_init_condition && !pair_deviation
    for (int i = 0; i < N; i+=2) {
	ptcl[i].x = separation * (i + 1);
	ptcl[i + 1].x = ptcl[i].x + 1;
    }
#endif
#if single_init_condition
    for (int i = 0; i < N; ++i)
	ptcl[i].x = separation * (i + 1);
#endif
    for (int i = 0; i < N; ++i) {
#if !pair_init_condition && !single_init_condition
	ptcl[i].x = (int)(L * Ran2()); // randomize initial position
#endif
#if pair_deviation
	ptcl[i].x = (int)(L / 2)-1;
#endif
	ptcl[i].Psi = vec((unsigned int)(M * Ran2()), (unsigned int)(M * Ran2())); // randomize initial particle states
	ptcl[i].number = i; // particle numbering
#if GRAPHICS
	for (int j = 0; j < y_height; ++j) { // input position to array
	    if (pctl_R[ptcl[i].x][j] == -1) {
		pctl_R[ptcl[i].x][j] = i;
		ptcl[i].y_cood = j;
		break;
	    }
	}
	int k = 0;
	while (k == 0) {
	    double lightcolor = 0;
	    for (int j = 0; j < 3; ++j) {
		color[i][j] = Ran2(); // randomize particle color
		lightcolor += color[i][j];
	    }
	    if (lightcolor > 1)
		k = 1;
	}
#endif
    }
    for (int i = 0; i < L; ++i)
	Phi[i] = vec((unsigned int)(M * Ran2()), (unsigned int)(M * Ran2())); // randomize initial bond states

    //binary tree memory allocation
    bin_tree = new int* [max_tree_lv];
    for (int i = 0; i < max_tree_lv; ++i) {
	int pow_temp = pow(2, max_tree_lv - i);
	bin_tree[i] = new int[pow_temp];
    }
    for (int i = 0; i < pow(2, max_tree_lv); ++i) {
	bin_tree[0][i] = (int)0;
    }

    for (int i = 0; i < N; ++i) {       
	bin_tree[0][2*i] = ptcl[i].f(0);    // refresh availability
	bin_tree[0][2*i+1] = ptcl[i].f(1);    // refresh availability
	ptcl[i].init_pos = ptcl[i].x; // get init position
    }
    bin_tree_init();
}

int main(int argc, char** argv) {
  
    int rand_seed = (unsigned)time(NULL);
    Randomize(0);   // Randomize(rand_seed) should be best
    //remove("X");
    if (argc >= 2) {
	istringstream a(argv[1]);
	a >> q;
	istringstream aa(argv[2]);
	aa >> N;
	istringstream aaa(argv[3]);
	aaa >> L;
	istringstream aaaa(argv[4]);
	aaaa >> t_max;
	istringstream aaaaa(argv[5]);
	aaaaa >> outn;
	string str_dir_name = argv[6];
	for (int i = 0; i < str_dir_name.length(); ++i)
	    file_name[i] = str_dir_name[i];
	string str_file_name = argv[7];
	for (int i = 0; i < str_file_name.length(); ++i)
	    file_name[str_dir_name.length()+i] = str_file_name[i];

	ofstream seedfile;
	string seed_filename = str_dir_name + (string)"seed";
	seedfile.open(seed_filename);
	seedfile << rand_seed;
	seedfile.close();
    }


    dt = (double)t_max / outn;
    //cout << "dt " << dt << endl;

	

    t_max += t_inc; // for testing

	
#if pair_deviation
    for (int i = 0; i < L; ++i) {
	for (int j = 0; j < outn; ++j)
	    histdata[i][j] = 0;
    }
#endif

#if pair_deviation
    int count = 0;
    for (int i = 0; i < pair_trials; ++i) {
	if (i >= count * 0.1 * pair_trials) { //print progress
	    cout << i*100/pair_trials << "%" << endl;
	    ++count;
	}
	t = 0;
	//cout << "new run" << endl;
	t_cnt = 0;
#endif
	main_init(); // initialize everything

#if GRAPHICS
	glutInitWindowSize(window_width, window_height);
	glutInit(&argc, argv);
	glutCreateWindow("frw animation");
	glOrtho(0, 300, 0, 100, 0, 1);  // window cover region: ( x1, to x2, y1, to y2, z1, to z2)
	glutDisplayFunc(display);
	initGL();
	timer(0);
	glutMainLoop();
#else
	while (t < t_max) {
  //          cout << "start t " << t << endl;
	    display();
	    //cout << "display t: " << t << endl;
	    //out_to_file();
	}

#if single_init_condition
	double count = 0;
	for (int i = 0; i < N; ++i) {
	    count += fabs(ptcl[i].x - ptcl[i].init_pos) / N;
	}
	//double average = (double)count / N;
	cout << count << endl;
#endif
#endif

#if pair_deviation
    }
    deviation.open(file_name);
    for (int j = 0; j < outn; ++j) {
	for (int jj = 0; jj < L; ++jj) {
	    deviation << histdata[jj][j] << " ";
	}
	deviation << endl;
    }
#endif

    //delete dynamically sized array
    for (int i = 0; i < max_tree_lv; ++i) {
	delete[] bin_tree[i];
    }
    delete[] bin_tree;

}
