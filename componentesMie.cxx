#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <sys/types.h>
#include <algorithm>
#include <complex_bessel.h>
#include <ctime>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include "cup_eV.H"
#include <pwd.h>
#define h  6.626068e-34
#define j2eV 6.24150636309e18
#define eV2j 1.60217733000103e-19

/*Compila con:
g++ -Wall -I/usr/local/include -L/usr/local/lib componentesMie.cxx -o Mfie -lgsl -lgslcblas -lm -lcomplex_bessel
*/
using namespace std;
using namespace sp_bessel;
double euler=2.718281828459045235360;
double eps_0=8.8541878176e-12;
nanosphere ns;

complex<double> x, m;

// Spherical unit vectors
double * hat_ere(double the, double phi){
    double* here;
    double her[3]={sin(the)*cos(phi), sin(the)*sin(phi), cos(the)};
    here=her;
    return here;
    }

double * hat_the(double the, double phi){
    double* hthe;
    double hth[3]={cos(the)*cos(phi), cos(the)*sin(phi), -sin(the)};
    hthe=hth;
    return hthe;
    }

double * hat_phi(double phi){
    double* hphi;
    double hph[3]={-sin(phi), cos(phi), 0};
    hphi=hph;
    return hphi;
    }

double dot(double a[3], double b[3]){
    double d=0.;
    for (int i=0; i<3; i++)
        d=d+a[i]*b[i];
    return d;
    }

// Bessel Functions.
complex<double> j (double order, complex<double> x){
    return sph_besselJ(order,x);
    }
complex<double> h1 (double order, complex<double> x){
    return sph_hankelH1(order, x);
    }

// Mie Angular Functions
double Pi (int n, double the){
    return boost::math::legendre_p(n, 1, cos(the))/sin(the);
    }
double Tu (int n, double the){
    return n*cos(the)*Pi(n,the)-(n+1)*Pi(n-1,the);
    }

// Spherical Harmonics.
complex<double> Y (int J, int M, double theta, double phi){
    return boost::math::spherical_harmonic(J,M,theta,phi);
    }

// Components of the vector spherical harmonics.
complex<double> Yt1(int J, int M, double theta, double phi){
    return 0.5*sqrt((J-M)*(J+M+1)/(J*(J+1)))*pow(euler,-ns.img*phi)*Y(J,M+1,theta,phi)-0.5*sqrt((J+M)*(J-M+1)/(J*(J+1)))*pow(euler,ns.img*phi)*Y(J,M-1,theta,phi);
    }

complex<double> Yp1(int J, int M, double theta, double phi){
    return ns.img*(double)M*Y(J,M,theta,phi)/(sqrt(J*(J+1))*sin(theta));
    }

complex<double> Yt0(int J, int M, double theta, double phi){
    return -(double)M*Y(J,M,theta,phi)/(sqrt(J*(J+1))*sin(theta));
    }

complex<double> Yp0(int J, int M, double theta, double phi){
    return  -ns.img*0.5*sqrt((J-M)*(J+M+1)/(J*(J+1)))*pow(euler, -ns.img*phi)*Y(J,M+1,theta,phi)+ns.img*0.5*sqrt((J+M)*(J-M+1)/(J*(J+1)))*pow(euler,ns.img*phi)*Y(J,M-1,theta,phi);
    }

// Dot product.
complex<double> Y1_Y1(int J1, int M1, int J2, int M2, double theta, double phi){
    return Yt1(J1,M1,theta,phi)*conj(Yt1(J2,M2,theta,phi)) + Yp1(J1,M1,theta,phi)*conj(Yp1(J2,M2,theta,phi));
    }
complex<double> Y1_Y0(int J1, int M1, int J2, int M2, double theta, double phi){
    return Yt1(J1,M1,theta,phi)*conj(Yt0(J2,M2,theta,phi)) + Yp1(J1,M1,theta,phi)*conj(Yp0(J2,M2,theta,phi));
    }
complex<double> Y0_Y1(int J1, int M1, int J2, int M2, double theta, double phi){
    return Yt0(J1,M1,theta,phi)*conj(Yt1(J2,M2,theta,phi)) + Yp0(J1,M1,theta,phi)*conj(Yp1(J2,M2,theta,phi));
    }
complex<double> Y0_Y0(int J1, int M1, int J2, int M2, double theta, double phi){
    return Yt0(J1,M1,theta,phi)*conj(Yt0(J2,M2,theta,phi)) + Yp0(J1,M1,theta,phi)*conj(Yp0(J2,M2,theta,phi));
    }
complex<double> Yn_Yn(int J1, int M1, int J2, int M2, double theta, double phi){
    return Y(J1,M1,theta,phi)*conj(Y(J2,M2,theta,phi));
    }

// Riccati-Bessel Functions.
complex<double> RBj (double order, complex<double> x){
    return x*sph_besselJ(order,x);
    }
complex<double> RBj_prime (double order, complex<double> x){
    return (x*sph_besselJ(order-1,x)-order*sph_besselJ(order,x));
    }

complex<double> RBh (double order, complex<double> x){
    return x*sph_hankelH1(order, x);
    }
complex<double> RBh_prime (double order, complex<double> x){
    return (x*sph_hankelH1(order-1,x)-order*sph_hankelH1(order,x));
    }



int main (int argc, char* argv[]){
    if (argc<2){
        cout<<"Usage "<<argv[0]<<"<omega en eV>"<<endl;
        exit(-1);
        }
    double omeeV=strtod(argv[1], NULL);
    double Esca[3], Einc[3];
    cout<<"ome = "<<omeeV<<" eV"<<endl;
    double lam, sp_r, eps2_0, G, omemi, omema;
    complex<double> k0, k1, k2, n1, n2, alpha1, alpha2;
    complex<double> k, eps1, eps2, alph, alph_mie,a[5], b[5], E0[5];
    char mtl[16], mdl[16], sol[16], active[16];


    fstream nano, comp, resu, cffc, flds, fdcs;

    nano.open("in/nanosphere_eV.dat",ios::in);
    flds.open("out/fields.dat", ios::out);
    fdcs.open("out/fldsmap.csv", ios::out);
    comp.open("out/compounds.dat", ios::out);
    //nano>>sp_r>>nada>>nada>>nada>>nada>>nada>>mtl>>mdl>>active>>sol;
    nano>>ns.r1>>ns.Dome>>ns.ome_0>>G>>omemi>>omema>>mtl>>mdl>>active>>sol;



    ns.r1=ns.r1*1.e-9;
    ns.init();
    ns.set_metal(mtl,mdl,1);
    eps2_0=ns.set_host(sol);

    ns.set_active(active);
    ns.ieps2_min=G;

 

    sp_r   = ns.r1;
    lam    = h*cc/(omeeV*eV2j);
    lam    = lam/sp_r;          // Normalized wavelength.
    cout<<sp_r<<endl;
    
    eps1=ns.metal(omeeV);
    eps2=ns.active(omeeV,eps2_0);
    

    n1=sqrt(eps1);            // Internal refractive index.
    n2=sqrt(eps2);            // External refractive index.
    m=n1/n2;                  // Relative refractive index.


    k0 = 2.*ns.pi/lam;
    k1 = 2.*ns.pi*n1/lam;
    k2 = 2.*ns.pi*n2/lam;


    x=k2; // k2*a con a=1


  /** Mie **/
      n1=sqrt(eps1);
      n2=sqrt(eps2);
      m=n1/n2;


      for (int n=1; n<=4;n++){
          a[n]=(m*RBj(n,m*x)*RBj_prime(n,x)-RBj(n,x)*RBj_prime(n,m*x))/(m*RBj(n,m*x)*RBh_prime(n,x)-RBh(n,x)*RBj_prime(n,m*x)); //Coefficient from Mie scattering theory.
          b[n]=(RBj(n,m*x)*RBj_prime(n,x)-m*RBj(n,x)*RBj_prime(n,m*x))/(RBj(n,m*x)*RBh_prime(n,x)-m*RBh(n,x)*RBj_prime(n,m*x)); //Coefficient from Mie scattering theory.
          
          cout<<n<<" "<<a[n]<<" "<<b[n]<<endl;
          }

    double xmax=5., ymax=5., zmax=5., x1, y1, z1, dx, dy, dz, erre, theta, phi;
    complex<double> CC, A_in, B_sc, A_sc;
    int Ncart=50;
//     int Ncart=10;


    double Es[3], Ei[3]; // Cartesian fields

    dx=2.*xmax/Ncart;
    dy=2.*ymax/Ncart;
    dz=2.*zmax/Ncart;
    
    fdcs<<"x,y,z,Ex,Ey,Ez"<<endl;

    for(int ii=0; ii<=Ncart; ii++){
        x1=-xmax+ii*dx;
        for(int ij=0; ij<=Ncart; ij++){
            y1=-ymax+ij*dy;
             for(int ik=0; ik<=Ncart; ik++){
                z1=-zmax+ik*dz;
                erre=sqrt(pow(x1,2)+pow(y1,2)+pow(z1,2));
                x=k2*erre;
                if (erre!=0.) {
                    theta=acos(z1/erre);
                    phi=atan2(y1,x1);
                    }
                    else {
                        theta=0;
                        phi=0.;
                        }

                complex<double> PRODUCTEr_sca = 0, PRODUCTEt_sca = 0, PRODUCTEp_sca = 0,
                                PRODUCTEr_inc = 0, PRODUCTEt_inc = 0, PRODUCTEp_inc = 0;

                for(int jj1 = 1; jj1<=4; jj1++){
                    CC=0.5*pow(ns.img,jj1)*sqrt(4.*ns.pi*(2*jj1+1));


 /*                   int mm1 = 1;

                    A_in = a_in(x, jj1, mm1);

                    B_sc=-a[jj1]*b_inc(jj1,1);
                    A_sc=-b[jj1]*a_inc(jj1,1);*/

                    E0[jj1]=pow(ns.img,jj1)*(2.*jj1+1)/(jj1*(jj1+1.));
                    // E_sca
//                     PRODUCTEr_sca+=CC*(alpha1*(ns.img/sp_r)*sqrt(jj1*(jj1+1))*B_sc*j(jj1,m*x))*Y(jj1,mm1,theta,phi);
//                     PRODUCTEt_sca+=CC*(alpha1*(ns.img/sp_r)*B_sc*RBj_prime(jj1,m*x)*Yt1(jj1,mm1,theta,phi)+a_in(x,jj1,mm1)*j(jj1,m*x)*Yt0(jj1,mm1,theta,phi));
//                     PRODUCTEp_sca+=CC*(alpha1*(ns.img/sp_r)*B_sc*RBj_prime(jj1,m*x)*Yp1(jj1,mm1,theta,phi)+a_in(x,jj1,mm1)*j(jj1,m*x)*Yp0(jj1,mm1,theta,phi));

                    PRODUCTEr_sca+=cos(phi)*sin(theta)*jj1*(jj1+1)*ns.img*E0[jj1]*a[jj1]*(h1(jj1,x)/x)*Pi(jj1,theta);
                    PRODUCTEt_sca+=cos(phi)*E0[jj1]*(ns.img*a[jj1]*(RBh_prime(jj1,x)/x)*Tu(jj1,theta)-b[jj1]*h1(jj1,x)*Pi(jj1,theta));
                    PRODUCTEp_sca+=sin(phi)*E0[jj1]*(b[jj1]*h1(jj1,x)*Tu(jj1,theta)-ns.img*a[jj1]*(RBh_prime(jj1,x)/x)*Pi(jj1,theta));

                    // E_inc
//                     PRODUCTEr_inc+=CC*alpha2*(ns.img/sp_r)*sqrt(jj1*(jj1+1))*b_inc(jj1,mm1)*j(jj1,x)*Y(jj1,mm1,theta,phi);
//                     PRODUCTEt_inc+=CC*(alpha2*(ns.img/sp_r)*b_inc(jj1,mm1)*RBj_prime(jj1,x)*Yt1(jj1,mm1,theta,phi)+a_inc(jj1,mm1)*j(jj1,x)*Yt0(jj1,mm1,theta,phi));
//                     PRODUCTEp_inc+=CC*(alpha2*(ns.img/sp_r)*b_inc(jj1,mm1)*RBj_prime(jj1,x)*Yp1(jj1,mm1,theta,phi)+a_inc(jj1,mm1)*j(jj1,x)*Yp0(jj1,mm1,theta,phi));



                    }
                if (erre>=1){
                    Einc[0]=PRODUCTEr_inc.real(); // radial component
                    Einc[1]=PRODUCTEt_inc.real(); // polar component
                    Einc[2]=PRODUCTEp_inc.real(); // azimuthal component

                    Esca[0]=PRODUCTEr_sca.real(); // radial component
                    Esca[1]=PRODUCTEt_sca.real(); // polar component
                    Esca[2]=PRODUCTEp_sca.real(); // azimuthal component
                    
                    for(int is = 0; is < 3; is++)
                        if (isnan(Esca[is])) Esca[is]=0.0;
                    
                    
                    } else {
                        Einc[0]=0; // radial component
                        Einc[1]=0; // polar component
                        Einc[2]=0; // azimuthal component

                        Esca[0]=0; // radial component
                        Esca[1]=0; // polar component
                        Esca[2]=0; // azimuthal component
                        }

                for (int ic = 0; ic < 3; ic++){
                    Es[ic] =  Esca[0]*hat_ere(theta,phi)[ic]+Esca[1]*hat_the(theta,phi)[ic]+Esca[2]*hat_phi(phi)[ic];
                    Ei[ic] =  Einc[0]*hat_ere(theta,phi)[ic]+Einc[1]*hat_the(theta,phi)[ic]+Einc[2]*hat_phi(phi)[ic];
                    }
                flds<<x1<<" "<<y1<<" "<<z1<<" "<<Es[0]<<" "<<Es[1]<<" "<<Es[2]<<" "<<Ei[0]<<" "<<Ei[1]<<" "<<Ei[2]<<endl;
                
                
                fdcs<<x1<<","<<y1<<","<<z1<<","<<Es[0]<<","<<Es[1]<<","<<Es[2]<<","<<erre<<endl;
                }
            flds<<endl;
            }
        flds<<endl;
        }
    return 0;
    }
