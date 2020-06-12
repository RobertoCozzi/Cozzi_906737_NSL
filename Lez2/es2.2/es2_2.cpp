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

double error(vector<double> & v, vector<double> & v2, int k){
	if(k==0) return 0;
	else return sqrt((v2[k]-v[k]*v[k])/double(k));
}

int main(int argc, char** argv){
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
	
	int steps=100, M=10000;
	double lat_pass=0;
	
	//discrete
	ofstream outdis("discrete.txt");
	double xd[M][3];
	vector<double> rd, rd2, rd_err;
	for(int i=0;i<M;i++){
		for(int j=0;j<3;j++){
			xd[i][j]=0;
		}
	}
	outdis<<0<<" "<<0<<" "<<0<<endl;
	rd.push_back(0);
	rd2.push_back(0);
	rd_err.push_back(0);
	for(int i=1;i<steps;i++){
		double r=0;
		double r2=0;
		for(int j=0;j<M;j++){
			double direction=rnd.Rannyu(-1,1);
			int component=rnd.Rannyu(1,4);
			if(direction<0) lat_pass=-1;
			if(direction>=0) lat_pass=1;
			xd[j][component]+=lat_pass;
			r+=xd[j][0]*xd[j][0]+xd[j][1]*xd[j][1]+xd[j][2]*xd[j][2];
			r2+=pow(xd[j][0]*xd[j][0]+xd[j][1]*xd[j][1]+xd[j][2]*xd[j][2],2);
		}
		r/=M;
		r2/=M;
		rd.push_back(sqrt(r));
		rd2.push_back(sqrt(r2));
		rd_err.push_back(error(rd,rd2,i));
		outdis<<i<<" "<<rd[i]<<" "<<rd_err[i]<<endl;
	}
	outdis.close();
	
	//continuum
	ofstream outcon("continuum.txt");
	double xc[M][3];
	vector<double> rc, rc2, rc_err;
	for(int i=0;i<M;i++){
		for(int j=0;j<3;j++){
			xc[i][j]=0;
		}
	}
	outcon<<0<<" "<<0<<" "<<0<<endl;
	rc.push_back(0);
	rc2.push_back(0);
	rc_err.push_back(0);
	for(int i=1;i<steps;i++){
		double r=0;
		double r2=0;
		for(int j=0;j<M;j++){
			double x,y,z;
			do{
				x=rnd.Rannyu(-1,1);
				y=rnd.Rannyu(-1,1);
				z=rnd.Rannyu(-1,1);	
			}while(x*x+y*y+z*z>1);
			double phi=atan2(y,x);
			double theta=acos(z);
			xc[j][0]+=sin(theta)*cos(phi);
			xc[j][1]+=sin(theta)*sin(phi);
			xc[j][2]+=cos(theta);
			r+=xc[j][0]*xc[j][0]+xc[j][1]*xc[j][1]+xc[j][2]*xc[j][2];
			r2+=pow(xc[j][0]*xc[j][0]+xc[j][1]*xc[j][1]+xc[j][2]*xc[j][2],2);
		}
		r/=M;
		r2/=M;
		rc.push_back(sqrt(r));
		rc2.push_back(sqrt(r2));
		rc_err.push_back(error(rc,rc2,i));
		outcon<<i<<" "<<rc[i]<<" "<<rc_err[i]<<endl;
	}
	outcon.close();
	
	rnd.SaveSeed();
	return 0;
}