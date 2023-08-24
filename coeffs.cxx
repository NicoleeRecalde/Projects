#include <iostream>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <complex_bessel.h>
#include "headers/myBessel.H"
#include "cup_eV.H"
#include "headers/mathNN.H"
#include "headers/gimme_p1.H"

/** Compila con:
g++ coeffs.cxx      -o cffs    -lgsl -lgslcblas -lm -lcomplex_bessel -larmadillo
**/

using namespace std;

std::complex<double>** coeffsMatriz(int order, double ome, std::complex <double> k, std::complex <double>  x, std::complex <double> m, std::complex <double> GG, std::complex <double>     E0, std::complex<double> OmeH, std::complex<double> OmeP, std::complex<double> GamP, double cnorm){

      std::complex<double>** M = 0;
      M = new std::complex<double>*[6];
      for (int j = 0; j < 6; j++)  M[j] = new std::complex<double>[6];

      std::complex<double>** p;
      p = new std::complex<double>*[4];
      for (int j = 0; j < 4; j++)  p[j] = new std::complex<double>[7];   
      
      p=gimme_p(order, ome, k,  x, m, GG, E0, OmeH, OmeP, GamP, cnorm);
      
      M[0][0] =   GG*p[0][0] + OmeH;
      M[0][1] =   GG*p[0][1];
      M[0][2] =   GG*p[0][2];
      M[0][3] =   GG*p[0][3];
      M[0][4] =   GG*p[0][4];
      M[0][5] =   GG*p[0][5];

      M[1][0] =   GG*p[1][0];
      M[1][1] =   GG*p[1][1] + OmeH;
      M[1][2] =   GG*p[1][2];
      M[1][3] =   GG*p[1][3];
      M[1][4] =   GG*p[1][4];
      M[1][5] =   GG*p[1][5];

      M[2][0] =   GamP*p[2][0];
      M[2][1] =   GamP*p[2][1];
      M[2][2] =   GamP*p[2][2] + OmeP;
      M[2][3] =   GamP*p[2][3];
      M[2][4] =   GamP*p[2][4];
      M[2][5] =   GamP*p[2][5];

      M[3][0] =   GamP*p[3][0];
      M[3][1] =   GamP*p[3][1];
      M[3][2] =   GamP*p[3][2];
      M[3][3] =   GamP*p[3][3] + OmeP;
      M[3][4] =   GamP*p[3][4];
      M[3][5] =   GamP*p[3][5];

      M[4][0] =   0;
      M[4][1] =   0;
      M[4][2] =   0;
      M[4][3] =   0;
      M[4][4] =   OmeH;
      M[4][5] =   0;

      M[5][0] =   0;
      M[5][1] =   0;
      M[5][2] =   0;
      M[5][3] =   0;
      M[5][4] =   0;
      M[5][5] =   OmeH;

      return M;
    }

std::complex<double> gimme_coef(int coef, int order, double ome, std::complex <double> k, std::complex<double> x, std::complex<double> m, std::complex<double> *q,std::complex <double> GG, std::complex <double> E0, std::complex<double> OmeH, std::complex<double> OmeP, std::complex<double> GamP, double cnorm){
    std::complex<double> a;

    std::complex<double>** p;
    p = new std::complex<double>*[7];
    for (int j = 0; j < 7; j++)  p[j] = new std::complex<double>[4];
    p=gimme_p(order,ome,k,  x, m, GG, E0, OmeH, OmeP, GamP, cnorm);

    a=p[coef][0]*q[0]+p[coef][1]*q[1]+p[coef][2]*q[2]+p[coef][3]*q[3]+p[coef][4]*q[4]+p[coef][5]*q[5]+p[coef][6]*E0;
    return a;
    }


std::complex<double>* inhomogeneousB( int order, double ome, std::complex <double> k, std::complex <double>  x, std::complex <double> m,  std::complex <double> GG, std::complex <double> E0, std::complex<double> OmeH, std::complex<double> OmeP, std::complex<double> GamP, double cnorm){
      std::complex<double>* B = 0;
      B = new std::complex<double>[6];

      std::complex<double>** p;
      p = new std::complex<double>*[7];
      for (int j = 0; j < 7; j++)  p[j] = new std::complex<double>[4];
      
      p=gimme_p(order,ome,k,  x, m, GG, E0, OmeH, OmeP, GamP, cnorm);

      B[0] = p[0][6]*GG*E0;
      B[1] = p[1][6]*GG*E0;
      B[2] = p[2][6]*GamP*E0;
      B[3] = p[3][6]*GamP*E0;
      B[4] = GG*E0;
      B[5] = GG*E0;

      return B;
    }


int main (int argc, char** argv){
    if (argv[1]==0){
        cout<<endl<<"  Usage: "<<argv[0]<<" <omega in eV>"<<endl<<endl;
        exit(0);
        }
    //constants
    double   ome, ome_21, omemi, omema, T2, gamd, lam, lam_nm, eps_b, eps_inf, G=0, ome_eV=3.2, cnorm;
    complex<double> eps1, eps2, m, x, n1, n2, k,k1, Nf=1, N=-1, E0=1;
    complex <double> OmeH, OmeP, GamP, GG, alph_mie;
    int order=1;
    char mtl[16], mdl[16], sol[16], active[16];

    ome_eV=atof(argv[1]);

    nanosphere ns;
    ns.init();

    ifstream nano("in/nanosphere_eV.dat");
    ofstream egva("out/eigenvalues.dat");
    ofstream fnct("out/anlfunc.dat");
    ofstream miec("out/anlmiec.dat");
    ofstream miea("out/anlmiea.dat");
    
    nano>>ns.r1>>ns.Dome>>ns.ome_0>>G>>omemi>>omema>>mtl>>mdl>>active>>sol;

    ns.r1=ns.r1*1.e-9;
    ns.set_metal(mtl,mdl,1);
    ns.set_active(active);

    eps_b=ns.set_host(sol);
    eps_inf=ns.eps_inf;
    eps1 = ns.metal(ome_eV);
    ns.G = G;
    eps2 = ns.active(ome_eV, eps_b);
    ome=ome_eV/ns.Ome_p;
    ome_21=ns.ome_0/ns.Ome_p;
    
    gamd=.5*ns.Gam_d/ns.Ome_p;
    T2=2.*ns.Ome_p/ns.Dome;

    n1=sqrt(eps1);
    n2=sqrt(eps2);
    m=n1/n2;

    lam = h*cc/(ome_eV*eV2j);
    lam_nm=lam*1.e9;
    cout<<"lam\t=\t"<<lam_nm<<" nm"<<endl;
    cout<<"r1\t=\t"<<ns.r1*1.e9<<" nm"<<endl;
//     cnorm = sqrt(lam)/(ns.r1*1.e9);
    lam =lam/ns.r1;
    cout<<"lam\t=\t"<<lam<<"*r1"<<endl;
    k = 2.*M_PI*n2/lam;
    cout<<"k\t=\t"<<k<<endl<<endl;    
    k1=m*k;
    x=k;

    OmeH= img*(ome-ome_21)-1./T2;
    GamP= 1./(2.*(gamd-img*ome));
    OmeP= ome*(ome+2.*img*gamd)*GamP;

    GG=Gwiggly(order, ns.G, E0, T2);
    
    cnorm = 1.; //cc*ns.Ome_p*eV2j/(h*ns.r1);
    
    complex <double> **coefis;
    coefis = new std::complex<double>*[6];    
    for (int j = 0; j < 6; j++)  coefis[j] = new std::complex<double>[6];    

    coefis=coeffsMatriz(order, ome, k, x, m, GG, E0, OmeH, OmeP, GamP, cnorm);

    complex <double> *inhomog;
    inhomog = new std::complex<double>[6];
  
    inhomog=inhomogeneousB(order, ome, k, x, m, GG, E0, OmeH, OmeP, GamP, cnorm);

    complex<double> *kap;
    kap = new std::complex<double>[6];

    
    kap = eigenvalues(coefis,6);

    egva<<"  "<<setw(8)<<setiosflags (ios::left)<<ome_eV<<                 // 1 ome
        "\t"<<setw(11)<<setiosflags (ios::left)<<real(kap[0])<< //ns.Ome_p*real(kap[0])<<  // 2 Re(kap1)
        "\t"<<setw(11)<<setiosflags (ios::left)<<imag(kap[0])<< //ns.Ome_p*imag(kap[0])<<  // 3 Im(kap1)
        "\t"<<setw(11)<<setiosflags (ios::left)<<real(kap[1])<< //ns.Ome_p*real(kap[1])<<  // 4 Re(kap2)
        "\t"<<setw(11)<<setiosflags (ios::left)<<imag(kap[1])<< //ns.Ome_p*imag(kap[1])<<  // 5 Im(kap2)
        "\t"<<setw(11)<<setiosflags (ios::left)<<real(kap[2])<< //ns.Ome_p*real(kap[2])<<  // 6 Re(kap3)
        "\t"<<setw(11)<<setiosflags (ios::left)<<imag(kap[2])<< //ns.Ome_p*imag(kap[2])<<  // 7 Im(kap3)
        "\t"<<setw(11)<<setiosflags (ios::left)<<real(kap[3])<< //ns.Ome_p*real(kap[3])<<  // 6 Re(kap4)
        "\t"<<setw(11)<<setiosflags (ios::left)<<imag(kap[3])<< //ns.Ome_p*imag(kap[3])<<  // 7 Im(kap4)
        "\t"<<setw(11)<<setiosflags (ios::left)<<real(kap[4])<< //ns.Ome_p*real(kap[4])<<  // 6 Re(kap5)
        "\t"<<setw(11)<<setiosflags (ios::left)<<imag(kap[4])<< //ns.Ome_p*imag(kap[4])<<  // 7 Im(kap5)
        "\t"<<setw(11)<<setiosflags (ios::left)<<real(kap[5])<< //ns.Ome_p*real(kap[5])<<  // 6 Re(kap6)
        "\t"<<setw(11)<<setiosflags (ios::left)<<imag(kap[5])<< //ns.Ome_p*imag(kap[5])<<  // 7 Im(kap6)
    endl;

    
    double omep = eV2j*ns.Ome_p/h;  //converto in Hz
    double t, T = 10., dt=.01; // tiempo total en picosegundos
    complex<double> *qss, *q, **EVE, *C, a1, d1;
    T=T*omep*1.e-12; // in ome_p
    dt=dt*omep*1.e-12;
    int i, Nt=T/dt;
    qss = new std::complex<double>[6];
    EVE = new std::complex<double>*[6];
    for(int i = 0; i < 6; i++)
        EVE[i] = new std::complex<double>[6];
    C   = new std::complex<double>[6];

    q   = new std::complex<double>[6];
    for(int i = 0; i < 6; i++) q[i] = std::complex<double> (0., 0.);

    qss = steady_state_solution(coefis, inhomog, 6);
    kap = eigenvalues(coefis,6);
    EVE = eigenvectors(coefis, 6);
    C   = constantes(coefis, inhomog, q, 6);
    i=0;
    while (i<=Nt){
        t=i*dt;
        i++;
        for(int ii = 0; ii < 6; ii++)
                    q[ii]   = qss[ii] + C[0]*EVE[ii][0]*exp(kap[0]*t)    //*(2.*M_PI*n2/(x*ns.r1*1.e9)) //*lam/pow(ns.r1*1.e9, 2)
                                      + C[1]*EVE[ii][1]*exp(kap[1]*t)    //*(2.*M_PI*n2/(x*ns.r1*1.e9))
                                      + C[2]*EVE[ii][2]*exp(kap[2]*t)    //*(2.*M_PI*n2/(x*ns.r1*1.e9))
                                      + C[3]*EVE[ii][3]*exp(kap[3]*t)    //*(2.*M_PI*n2/(x*ns.r1*1.e9))
                                      + C[4]*EVE[ii][4]*exp(kap[4]*t)    //*(2.*M_PI*n2/(x*ns.r1*1.e9))
                                      + C[5]*EVE[ii][5]*exp(kap[5]*t);   //*(2.*M_PI*n2/(x*ns.r1*1.e9))

            fnct<<"  "<<setw(8)<<setiosflags (ios::left)<<t*1.e+12/omep<< //  1 time (ps)
                "\t"<<setw(13)<<setiosflags (ios::left)<<real(q[0])<<         //  2 Re(alph)
                "\t"<<setw(13)<<setiosflags (ios::left)<<imag(q[0])<<         //  3 Im(alph)
                "\t"<<setw(13)<<setiosflags (ios::left)<<real(q[1])<<         //  4 Re(beta)
                "\t"<<setw(13)<<setiosflags (ios::left)<<imag(q[1])<<         //  5 Im(beta)
                "\t"<<setw(13)<<setiosflags (ios::left)<<real(q[2])<<         //  6 Re(kappa)
                "\t"<<setw(13)<<setiosflags (ios::left)<<imag(q[2])<<         //  7 Im(kappa)
                "\t"<<setw(13)<<setiosflags (ios::left)<<real(q[3])<<         //  8 Re(delta)
                "\t"<<setw(13)<<setiosflags (ios::left)<<imag(q[3])<<         //  9 Im(delta)
                "\t"<<setw(13)<<setiosflags (ios::left)<<real(q[4])<<         // 10 Re(eta)
                "\t"<<setw(13)<<setiosflags (ios::left)<<imag(q[4])<<         // 11 Im(eta)
                "\t"<<setw(13)<<setiosflags (ios::left)<<real(q[5])<<         // 12 Re(zeta)
                "\t"<<setw(13)<<setiosflags (ios::left)<<imag(q[5])<<         // 13 Im(zeta)
                endl;
            
            a1 = gimme_coef(0,order,ome, k, x, m, q, GG, E0, OmeH, OmeP, GamP, cnorm)/E0;
            miec<<"  "<<setw(8)<<setiosflags (ios::left)<<t/omep*1.e+12<< //  1 time (ps)
                "\t"<<setw(13)<<setiosflags (ios::left)<<real(a1)<<        //  2 Re(a1)
                "\t"<<setw(13)<<setiosflags (ios::left)<<imag(a1)<<        //  3 Im(a1)
                endl;
            
            alph_mie = 6.*M_PI*img*a1/pow(k,3);
            
            miea<<"  "<<setw(8)<<setiosflags (ios::left)<<t/omep*1.e+12<< //  1 time (ps)
                "\t"<<setw(13)<<setiosflags (ios::left)<<real(alph_mie)/(4*M_PI)<<        //  2 Re(alpha)
                "\t"<<setw(13)<<setiosflags (ios::left)<<imag(alph_mie)/(4*M_PI)<<        //  3 Im(alpha)
                endl;
            
        }
//     cout<<"final anltcl a["<<order<<"] = "<<a1<<endl<<endl;
    cout<<ome_eV<<" "<<alph_mie.real()/(4*M_PI)<<" "<<alph_mie.imag()/(4*M_PI)<<endl;
    a1 = gimme_coef(0,order,ome, k, x, m, qss, GG, E0, OmeH, OmeP, GamP, cnorm)/E0;

//     cout<<"steady state a["<<order<<"] = "<<a1<<endl;
// 

//     
//     cout<<"OmeH = "<<OmeH<<endl;
//     cout<<"OmeP = "<<OmeP<<endl<<endl;
    
//     for(int ii = 0; ii < 6; ii++) cout<<"qss["<<ii<<"] = "
//                                       <<std::setw(28)<<std::setiosflags (std::ios::left)<<qss[ii]
//                                       <<" ---   q["<<ii<<"] = "
//                                       <<std::setw(28)<<std::setiosflags (std::ios::left)<<q[ii]<<endl;
    return 0;
    }
