#include <iostream>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include "headers/math33.H"
#include "headers/single.H"
#include "headers/cup.H"


/*
g++ -Wall -I/usr/local/include -L/usr/local/lib mie.cxx -o mie -lgsl -lgslcblas -lm -larmadillo
*/

using namespace std;

int main(int argc, char** argv){
    double   omeeV, omemi, omema, eps_b, E0, T, tpump, pG;
    complex<double> eps1, eps2, alph, alph_num, alph_anl;
    double tau2, gamd, ome0, ome, *p0;
    complex<double> GG, OmeH, OmeP, GamP, **A;
    
    char mtl[16], mdl[16], sol[16], active[16];
    if (argv[1]==0){
        cout<<endl<<"  Usage: "<<argv[0]<<" <omega in eV>"<<endl<<endl;
        exit(0);
        }
    omeeV=atof(argv[1]);
    
    nanosphere ns;    
    fstream nano, time;

    nano.open("in/nanosphere_eV.dat", ios::in);
    time.open("in/time.dat", ios::in);


    nano>>ns.r1>>ns.Dome>>ns.ome_0>>ns.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0;
    time>>T>>tpump;    
    
    pG = fabs(ns.G);
    
    ns.init();
    ns.set_metal(mtl,mdl,1);
    eps_b=ns.set_host(sol);
    ns.set_active(active);
    
    tau2 = 2.*ns.Ome_p/ns.Dome;
    GG   = -ns.img*pG/tau2;
    gamd = .5*ns.Gam_d/ns.Ome_p;
    ome0 = ns.ome_0/ns.Ome_p;
    
    ome = omeeV/ns.Ome_p;
    
    OmeH = ns.img*(ome-ome0)-1./tau2;
    GamP = 1./(2.*(gamd-ns.img*ome));
    OmeP = (ome*(ome+2.*ns.img*gamd))*GamP;
                
    p0 = pcfc0(ns.eps_inf, eps_b);
    
    
    A = new std::complex<double>*[3];
    for(int i = 0; i <= 2; i++)
        A[i] = new std::complex<double>[3];
  
    A=coefficients(OmeH, OmeP,  GG, GamP, p0);
    
    cout<<endl;
    for(int i = 0; i <= 2; i++){
          for(int j = 0; j <= 2; j++) cout<<std::setw(26)<<std::setiosflags (std::ios::left)<<A[i][j]
                                          <<std::setw(8)<<std::setiosflags (std::ios::left)<<" ";
          cout<<endl;
          }
    cout<<endl;
  return 0;
  }
    

