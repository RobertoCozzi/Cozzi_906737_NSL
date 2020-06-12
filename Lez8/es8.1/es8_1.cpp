#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "random.h"
using namespace std;

//random variables
int seed[4];
Random rnd;

//observables
const int m_props=100;
int n_props, iv, ik, ie;
double walker[m_props];

//averages
double blk_average[m_props], glob_av[m_props], glob_av2[m_props], blk_norm;
int accept, accept_opt;
double stima_epot, stima_ekin, stima_etot;
double err_epot, err_ekin, err_etot;
double H, H_old;

//position
double x, x_new, delta_x;
ofstream prob;

//parameters
double mu, sigma, mu_old, sigma_old, delta_m, delta_s;

//simulation
int nstep, nblock, nopt, print;

//functions
void Initialize(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void Measure(void);
double Error(double, double, int);
double psi(double);
double K_x(double);
double V_x(double);


int main(){
	Initialize();
	print=0;
	for(int iopt=0;iopt<nopt;iopt++){
		if(iopt%50==0) cout<<"optimization step: "<<iopt<<endl;
		if(iopt==nopt-1) print=1;
		if(iopt!=0){
			mu_old=mu;
			sigma_old=sigma;
			mu=rnd.Rannyu(mu-delta_m,mu+delta_m);
			sigma=rnd.Rannyu(sigma-delta_s,sigma+delta_s);
		}
		for(int iblk=1;iblk<=nblock;iblk++){
			Reset(iblk);
			for(int istep=1;istep<=nstep;istep++){
				Move();
				Measure();
				Accumulate();
			}
			Averages(iblk);
		}
		H=glob_av[ie]/nblock;
		if(iopt==0) H_old=H;
		if(iopt!=0){
			if(H<H_old){
				accept_opt++;
				H_old=H;
				delta_m/=2.;
				delta_s/=2.;
			}
			else{
				mu=mu_old;
				sigma=sigma_old;
			}
		}
	}
	cout<<"Optimized parameters"<<endl;
	cout<<"   mu   ="<<mu<<endl;
	cout<<"   sigma="<<sigma<<endl;
	cout<<"accepted optimization steps "<<accept_opt<<" of "<<nopt<<endl;
	
	return 0;
}

void Initialize(){
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
      Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	
	cout<<"Variational Monte Carlo simulation for the Ground State 1D particle with external potential:"<<endl;
	cout<<"     V(x)=x^4-2.5x^2"<<endl;
	cout<<"using the trial function:"<<endl;
	cout<<"     Psi_T(x)=exp(-(x-mu)^2/2sigma^2)+exp(-(x+mu)^2/2sigma^2)"<<endl<<endl;
	
	nblock=200;
	nstep=1000;
	nopt=500;
	x=1;
	delta_x=2;
	mu=0.8;
	delta_m=1.;
	sigma=1;
	delta_s=1.;
	accept=0;
	accept_opt=0;
	
	cout<<"Number of blocks "<<nblock<<endl;
	cout<<"Number of steps per block "<<nstep<<endl;
	cout<<"Initial parameters: mu="<<mu<<" sigma="<<sigma<<endl;
	cout<<"Optimization steps "<<nopt<<endl;
	cout<<"of lenght: dmu="<<delta_m<<" dsigma="<<delta_s<<endl;
	cout<<"Starting position x="<<x<<" with step lenght dx="<<delta_x<<endl<<endl;
	
	iv=0;
	ik=1;
	ie=2;
	n_props=3;
	
	return;
}

void Move(){
	x_new=rnd.Rannyu(x-delta_x,x+delta_x);
	double p_old=pow(psi(x),2);
	double p_new=pow(psi(x_new),2);
	double p=p_new/p_old;
	double b=rnd.Rannyu();
	if(p>b){
		x=x_new;
		accept++;
	}
	if(print==1){
		prob.open("pos.0",ios::app);
		prob<<x<<endl;
		prob.close();
	}
}

void Reset(int iblk){ //Reset block averages
   if(iblk == 1){
		for(int i=0; i<n_props; ++i){
			glob_av[i] = 0;
			glob_av2[i] = 0;
		}
	}

	for(int i=0; i<n_props; ++i){
		blk_average[i] = 0;
	}
	blk_norm = 0;
	accept = 0;
}

void Accumulate(void){ //Update block averages
	for(int i=0; i<n_props; ++i){
		blk_average[i] = blk_average[i] + walker[i];
	}
	blk_norm = blk_norm + 1.0;
}

void Measure(){
	double v=V_x(x);
	double k=K_x(x);
	
	walker[iv]=v;
	walker[ik]=k;
	walker[ie]=v+k;
}

void Averages(int iblk){
	ofstream Epot, Ekin, Etot;
	
	stima_epot=blk_average[iv]/blk_norm;
	glob_av[iv]+=stima_epot;
	glob_av2[iv]+=stima_epot*stima_epot;
	err_epot=Error(glob_av[iv],glob_av2[iv],iblk);
	
	stima_ekin=blk_average[ik]/blk_norm;
	glob_av[ik]+=stima_ekin;
	glob_av2[ik]+=stima_ekin*stima_ekin;
	err_ekin=Error(glob_av[ik],glob_av2[ik],iblk);
	//for(int i=0; i<100; i++) cout<<iblk<<" "<<stima_ekin<<" "<<mu<<" "<<sigma<<endl;
	
	stima_etot=blk_average[ie]/blk_norm;
	glob_av[ie]+=stima_etot;
	glob_av2[ie]+=stima_etot*stima_etot;
	err_etot=Error(glob_av[ie],glob_av2[ie],iblk);
	
	if(print==1){
		Epot.open("output.epot.0",ios::app);
		Ekin.open("output.ekin.0",ios::app);
		Etot.open("output.etot.0",ios::app);
		
		Epot<<iblk<<" "<<stima_epot<<" "<<glob_av[iv]/(double)iblk<<" "<<err_epot<<endl;
		Ekin<<iblk<<" "<<stima_ekin<<" "<<glob_av[ik]/(double)iblk<<" "<<err_ekin<<endl;
		Etot<<iblk<<" "<<stima_etot<<" "<<glob_av[ie]/(double)iblk<<" "<<err_etot<<endl;
		
		Epot.close();
		Ekin.close();
		Etot.close();
	}
}

double Error(double sum, double sum2, int iblk){
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

double psi(double y){
	double t1=pow((y-mu)/sigma,2);
	double t2=pow((y+mu)/sigma,2);
	return exp(-t1/2)+exp(-t2/2);
}

double K_x(double y){
	double t1=pow((y-mu)/sigma,2);
	double t2=pow((y+mu)/sigma,2);
	return -0.5*((t1-1)*exp(-t1/2)+(t2-1)*exp(-t2/2))/(sigma*sigma);
}

double V_x(double y){
	return pow(y,4)-2.5*pow(y,2);
}