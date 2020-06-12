#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <tuple>
#include "random.h"
#include "functions.h"
using namespace std;

int main (int argc, char *argv[]){
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
	
	double L=0.8;
	double d=1;
	int N=100;
	int M=10000;
	vector<double> pi, pi2, pi_sum, pi2_sum, pi_err;
	for(int i=0;i<N;i++){
		int hit=0;
		int tot=0;
		for(int j=0;j<M;j++){
			double x=rnd.Rannyu(0,d);
			double xend, yend;
			do{
				xend=rnd.Rannyu();
				yend=rnd.Rannyu();
			}while((xend*xend+yend*yend)>1);
			double angle=atan2(yend,xend);
			if((x+L*sin(angle)/2)>d or (x-L*sin(angle)/2)<0) hit++;
			tot++;
		}
		double PI=2*L*tot/(d*hit);
		pi.push_back(PI);
		pi2.push_back(PI*PI);
	}
	pi_sum=get<0>(BlockStat(pi,pi2,N));
	pi_err=get<1>(BlockStat(pi,pi2,N));
	PrintError("pi.txt",pi_sum,pi_err);
	
	return 0;
}