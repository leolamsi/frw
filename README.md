# Emergent Facilitation by Random Constraints in a Facilitated Random Walk Model of Glass  
[PDF of paper](https://arxiv.org/pdf/2412.08986)

## Authors  
Leo S.I. Lam, Hai-Yao Deng, Wei-Bing Zhang, Udoka Nwankwo, Chu Xiao, Cho-Tung Yip, Chun-Shing Lee, Haihui Ruan and Chi-Hang Lam  

## TLDR
- This paper introduces a new glass lattice model called the **Facilitated Random Walk (FRW)**, which incorporates random local kinetic constraints.  
- FRW is a coarse-grained version of the distinguishable particle lattice model (DPLM). FRW exhibits essential glassy properties without explicit energetic interactions.
- The code in this repository is used to run the main simulation (frw.cpp, with .py or .sh drivers) and to run the analysis (.py or matlab or .sh) 

## Abstract
The physics of glass has been a significant topic of interest for decades. Dynamical facilitation is widely believed to be an important characteristic of glassy dynamics, but the precise mechanism is still under debate. We propose a lattice model of glass called the facilitated random walk (FRW). Each particle performs continuous time random walk in the presence of its own random local kinetic constraints. The particles do not interact energetically. Instead, they interact kinetically with a hopping rate resampling rule under which motions of a particle can randomly perturb the local kinetic constraints of other particles. This dynamic interaction is reversible, following a rate restoration rule. A step-by-step reversal of the particle motions exactly restore the previous constraints, modeling randomness quenched in the configuration space of glass. The model exhibits  stretched exponential relaxation and dynamical heterogeneity typical of glasses. Despite the lack of explicit facilitation rule, the FRW shows facilitation behaviors closely analogous to those of the kinetically constrained  models (KCM).  The FRW is  a coarse-grained version of the distinguishable particle lattice model (DPLM) and this  exemplifies that compatible defect and atomistic models can complement each other on the study of glass.

## Project overview
### Key Concepts of Model  
- **Kinetically Constrained Models (KCM)**: Models that incorporate local constraints to simulate glassy behavior.  
- **Continuous-Time Random Walk (CTRW)**: A framework for modeling particle movements in a lattice with random constraints.  
- **Random Constraints**:  Each particle performs a random walk influenced by local kinetic constraints.
- **Energetically trivial**: Availability of particle hops are determined by a probability q, neither energy nor temperature is considered in this model.
- **Particle Facilitation**: Hopping rates can be resampled based on the movements of other particles, introducing a form of dynamical interaction.  
- **Glassy Material**: The model exhibits behaviors typical of glasses, such as stretched exponential relaxation and dynamical heterogeneity.  

### Simulation and Results  
- The FRW model was simulated in one dimension under various conditions (density rho and unblocking probability q).  
- Observations included:  
  - **Dynamical Facilitation**: Emergence of mobile groups of particles that can traverse the lattice, enhancing mobility.  
  - **Mean Square Displacement (MSD)**: Characteristic plateaus indicating temporary trapping of particles, typical in glassy systems.  
  - **Self-Intermediate Scattering Function (SISF)**: Fitted to the Kohlrausch-Williams-Watts (KWW) stretched exponential function, confirming glassy dynamics.
 
## Code
### Compiling main program frw.cpp
Make sure the following files are within your compilation directory: frw.cpp, MathAux.cpp, MathAux.h, FileIO.h. Use the command "g++ -O3 frw.cpp MathAux.cpp -o <executable name>". We usually use frw.exe for our executable name.
### Running frw.exe (or your executable file name)
We use the python script msd_run.py to run our scripts. (Note: Many file names start with "msd" . This is only because the msd graph was the first set of results investigated. These programs work for all applications. Please rename to your own liking. Also redesign your own data saving architecture to your liking.)
### perco_frw.cpp
Only a few modifications made from frw.cpp to focus on working on the percolation threshold. The main difference is that one of the #define variables within the frw.cpp are changed (pair_init_condition). It is mainly used for debugging and more straightforward compiling.
### #define variables in frw.cpp
There is a display mode where the frw is visually represented (OpenGL is needed). There is also a single particle condition where one particle is initiated in the middle of the lattice. There is also a pair condition (explained above in perco_frw.cpp section). The pair deviation mode outputs a histogram of particle positions over time for outputting two-particle standard deviation over time.

## Acknowledgments  
- This work was supported by General Research Fund ofHong Kong (Grants 15210622 and 15303220). 
