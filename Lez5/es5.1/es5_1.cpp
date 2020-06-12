#include "functions.h"
#include "random.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <tuple>
using namespace std;

double psi100(double x, double y, double z){
	double rad=sqrt(x*x+y*y+z*z);
	double psi=exp(-rad)/sqrt(M_PI);
	return pow(psi,2);
}

double psi210(double x, double y, double z){
	double rad=sqrt(x*x+y*y+z*z);
	double psi=exp(-rad/2.)*sqrt(2./M_PI)*z/8;
	return pow(psi,2);
}

double A_P(double p, double p_prev){
	return min(1.,p/p_prev);
}

int main(){
	Random rnd;
	int seed[4];
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
	
	int M=1000000, blocks=500;
	int L=M/blocks;
	
	//uniform
	//psi_1,0,0
	vector<double> ru100, ru100_2, ru100_ave, ru100_err;
	double x_start=20, y_start=20, z_start=20;
	double x_prev=x_start, y_prev=y_start, z_prev=z_start;
	double x,y,z;
	int ac_u100=0;
	double step_u100=1;
	for(int i=0;i<10000;i++){												//equilibration
		double dx, dy, dz;
		ofstream uni100;
		do{
			x=rnd.Rannyu(x_prev-step_u100,x_prev+step_u100);
			y=rnd.Rannyu(y_prev-step_u100,y_prev+step_u100);
			z=rnd.Rannyu(z_prev-step_u100,z_prev+step_u100);
			dx=x-x_prev;
			dy=y-y_prev;
			dz=z-z_prev;
		}while((dx*dx+dy*dy+dz*dz)<step_u100*step_u100);
		double a=A_P(psi100(x,y,z),psi100(x_prev,y_prev,z_prev));
		double b=rnd.Rannyu();
		if(b<a){
			x_prev=x;
			y_prev=y;
			z_prev=z;
		}
		uni100.open("uni100.txt",ios::app);
		uni100<<x_prev<<" "<<y_prev<<" "<<z_prev<<endl;
		uni100.close();
	}
	for(int i=0;i<blocks;i++){												//actual simulation
		double sum=0;
		for(int j=0;j<L;j++){
			double dx, dy, dz;
			do{
				x=rnd.Rannyu(x_prev-step_u100,x_prev+step_u100);
				y=rnd.Rannyu(y_prev-step_u100,y_prev+step_u100);
				z=rnd.Rannyu(z_prev-step_u100,z_prev+step_u100);
				dx=x-x_prev;
				dy=y-y_prev;
				dz=z-z_prev;
			}while((dx*dx+dy*dy+dz*dz)<step_u100*step_u100);
			double a=A_P(psi100(x,y,z),psi100(x_prev,y_prev,z_prev));
			double b=rnd.Rannyu();
			if(b<a){
				x_prev=x;
				y_prev=y;
				z_prev=z;
				ac_u100++;
			}
			sum+=sqrt(x_prev*x_prev+y_prev*y_prev+z_prev*z_prev);
		}
		sum/=L;
		ru100.push_back(sum);
		ru100_2.push_back(sum*sum);
	}
	ru100_ave=get<0>(BlockStat(ru100,ru100_2,blocks));
	ru100_err=get<1>(BlockStat(ru100,ru100_2,blocks));
	PrintError("PSI100_uniform.txt",ru100_ave,ru100_err);
	cout<<"Accepted moves "<<ac_u100<<" of 10e6"<<endl;
	
	//psi_2,1,0
	vector<double> ru210, ru210_2, ru210_ave, ru210_err;
	x_start=20, y_start=20, z_start=20;
	x_prev=x_start, y_prev=y_start, z_prev=z_start;
	int ac_u210=0;
	double step_u210=2;
	for(int i=0;i<10000;i++){												//equilibration
		double dx, dy, dz;
		ofstream uni210;
		do{
			x=rnd.Rannyu(x_prev-step_u210,x_prev+step_u210);
			y=rnd.Rannyu(y_prev-step_u210,y_prev+step_u210);
			z=rnd.Rannyu(z_prev-step_u210,z_prev+step_u210);
			dx=x-x_prev;
			dy=y-y_prev;
			dz=z-z_prev;
		}while((dx*dx+dy*dy+dz*dz)<step_u210*step_u210);
		double a=A_P(psi210(x,y,z),psi210(x_prev,y_prev,z_prev));
		double b=rnd.Rannyu();
		if(b<a){
			x_prev=x;
			y_prev=y;
			z_prev=z;
		}
		uni210.open("uni210.txt",ios::app);
		uni210<<x_prev<<" "<<y_prev<<" "<<z_prev<<endl;
		uni210.close();
	}
	for(int i=0;i<blocks;i++){												//actual simulation
		double sum=0;
		for(int j=0;j<L;j++){
			double dx, dy, dz;
			do{
				x=rnd.Rannyu(x_prev-step_u210,x_prev+step_u210);
				y=rnd.Rannyu(y_prev-step_u210,y_prev+step_u210);
				z=rnd.Rannyu(z_prev-step_u210,z_prev+step_u210);
				dx=x-x_prev;
				dy=y-y_prev;
				dz=z-z_prev;
			}while((dx*dx+dy*dy+dz*dz)<step_u210*step_u210);
			double a=A_P(psi210(x,y,z),psi210(x_prev,y_prev,z_prev));
			double b=rnd.Rannyu();
			if(b<a){
				x_prev=x;
				y_prev=y;
				z_prev=z;
				ac_u210++;
			}
			sum+=sqrt(x_prev*x_prev+y_prev*y_prev+z_prev*z_prev);
		}
		sum/=L;
		ru210.push_back(sum);
		ru210_2.push_back(sum*sum);
	}
	ru210_ave=get<0>(BlockStat(ru210,ru210_2,blocks));
	ru210_err=get<1>(BlockStat(ru210,ru210_2,blocks));
	PrintError("PSI210_uniform.txt",ru210_ave,ru210_err);
	cout<<"Accepted moves "<<ac_u210<<" of 10e6"<<endl;
	
	//Gaussian
	//psi_1,0,0
	vector<double> rg100, rg100_2, rg100_ave, rg100_err;
	x_start=20, y_start=20, z_start=20;
	x_prev=x_start, y_prev=y_start, z_prev=z_start;
	int ac_g100=0;
	double step_g100=1;
	for(int i=0;i<10000;i++){												//equilibration
		double dx, dy, dz;
		ofstream gauss100;
		do{
			x=rnd.Gauss(x_prev,0.5*step_g100);
			y=rnd.Gauss(y_prev,0.5*step_g100);
			z=rnd.Gauss(z_prev,0.5*step_g100);
			dx=x-x_prev;
			dy=y-y_prev;
			dz=z-z_prev;
		}while((dx*dx+dy*dy+dz*dz)<step_g100*step_g100);
		double a=A_P(psi100(x,y,z),psi100(x_prev,y_prev,z_prev));
		double b=rnd.Rannyu();
		if(b<a){
			x_prev=x;
			y_prev=y;
			z_prev=z;
		}
		gauss100.open("gauss100.txt",ios::app);
		gauss100<<x_prev<<" "<<y_prev<<" "<<z_prev<<endl;
		gauss100.close();
	}
	for(int i=0;i<blocks;i++){												//actual simulation
		double sum=0;
		for(int j=0;j<L;j++){
			double dx, dy, dz;
			do{
				x=rnd.Gauss(x_prev,0.5*step_g100);
				y=rnd.Gauss(y_prev,0.5*step_g100);
				z=rnd.Gauss(z_prev,0.5*step_g100);
				dx=x-x_prev;
				dy=y-y_prev;
				dz=z-z_prev;
			}while((dx*dx+dy*dy+dz*dz)<step_g100*step_g100);
			double a=A_P(psi100(x,y,z),psi100(x_prev,y_prev,z_prev));;
			double b=rnd.Rannyu();
			if(b<a){
				x_prev=x;
				y_prev=y;
				z_prev=z;
				ac_g100++;
			}
			sum+=sqrt(x_prev*x_prev+y_prev*y_prev+z_prev*z_prev);
		}
		sum/=L;
		rg100.push_back(sum);
		rg100_2.push_back(sum*sum);
	}
	rg100_ave=get<0>(BlockStat(rg100,rg100_2,blocks));
	rg100_err=get<1>(BlockStat(rg100,rg100_2,blocks));
	PrintError("PSI100_gauss.txt",rg100_ave,rg100_err);
	cout<<"Accepted moves "<<ac_g100<<" of 10e6"<<endl;
	
	//psi_2,1,0
	vector<double> rg210, rg210_2, rg210_ave, rg210_err;
	x_start=20, y_start=20, z_start=20;
	x_prev=x_start, y_prev=y_start, z_prev=z_start;
	int ac_g210=0;
	double step_g210=2;
	for(int i=0;i<10000;i++){												//equilibration
		double dx, dy, dz;
		ofstream gauss210;
		do{
			x=rnd.Gauss(x_prev,0.5*step_g210);
			y=rnd.Gauss(y_prev,0.5*step_g210);
			z=rnd.Gauss(z_prev,0.5*step_g210);
			dx=x-x_prev;
			dy=y-y_prev;
			dz=z-z_prev;
		}while((dx*dx+dy*dy+dz*dz)<step_g210*step_g210);
		double a=A_P(psi210(x,y,z),psi210(x_prev,y_prev,z_prev));
		double b=rnd.Rannyu();
		if(b<a){
			x_prev=x;
			y_prev=y;
			z_prev=z;
		}
		gauss210.open("gauss210.txt",ios::app);
		gauss210<<x_prev<<" "<<y_prev<<" "<<z_prev<<endl;
		gauss210.close();
	}
	for(int i=0;i<blocks;i++){
		double sum=0;
		for(int j=0;j<L;j++){
			double dx, dy, dz;
			do{
				x=rnd.Gauss(x_prev,0.5*step_g210);
				y=rnd.Gauss(y_prev,0.5*step_g210);
				z=rnd.Gauss(z_prev,0.5*step_g210);
				dx=x-x_prev;
				dy=y-y_prev;
				dz=z-z_prev;
			}while((dx*dx+dy*dy+dz*dz)<step_g210*step_g210);
			double a=A_P(psi210(x,y,z),psi210(x_prev,y_prev,z_prev));
			double b=rnd.Rannyu();
			if(b<a){
				x_prev=x;
				y_prev=y;
				z_prev=z;
				ac_g210++;
			}
			sum+=sqrt(x_prev*x_prev+y_prev*y_prev+z_prev*z_prev);
		}
		sum/=L;
		rg210.push_back(sum);
		rg210_2.push_back(sum*sum);
	}
	rg210_ave=get<0>(BlockStat(rg210,rg210_2,blocks));
	rg210_err=get<1>(BlockStat(rg210,rg210_2,blocks));
	PrintError("PSI210_gauss.txt",rg210_ave,rg210_err);
	cout<<"Accepted moves "<<ac_g210<<" of 10e6"<<endl;
	
	return 0;
}