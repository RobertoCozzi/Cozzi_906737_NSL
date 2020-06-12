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
	
	int M=10000, N=100;
	double L=M/N;
	//uniform
	vector<double> r;
	for(int i=0;i<M;i++){
		double x=rnd.Rannyu();
		r.push_back(x);
	}
	//average
	vector<double> ave, ave2, ave_sum, ave_err;
	for(int i=0;i<N;i++){
		double sum=0;
		for(int j=0;j<L;j++){
			int k=j+i*L;
			sum+=(M_PI/2)*cos(M_PI*r[k]/2);
		}
		sum/=L;
		ave.push_back(sum);								//block average
		ave2.push_back(sum*sum);
	}
	ave_sum=get<0>(BlockStat(ave,ave2,N));
	ave_err=get<1>(BlockStat(ave,ave2,N));
	//variance
	vector<double> var, var2, var_sum, var_err;
	for(int i=0;i<N;i++){
		double sum=0;
		for(int j=0;j<L;j++){
			int k=j+i*L;
			sum+=pow((M_PI/2)*cos(M_PI*r[k]/2)-1,2);
		}
		sum/=L;
		var.push_back(sum);								//block variance
		var2.push_back(sum*sum);
	}
	var_sum=get<0>(BlockStat(var,var2,N));
	var_err=get<1>(BlockStat(var,var2,N));
	
	PrintError("Uniform.txt",ave_sum,ave_err);
	
	//importance sampling  p(x)=2*(1-x)
	vector<double> ri;
	for(int i=0;i<M;i++){
		double x=1+sqrt(1-rnd.Rannyu());
		ri.push_back(x);
	}
	vector<double> avei, avei2, avei_sum, avei_err;
	for(int i=0;i<N;i++){
		double sum=0;
		for(int j=0;j<L;j++){
			int k=j+i*L;
			sum+=(M_PI/2)*cos((M_PI*ri[k]/2))/(2*(1-ri[k]));
		}
		sum/=L;
		avei.push_back(sum);								//block average
		avei2.push_back(sum*sum);
	}
	avei_sum=get<0>(BlockStat(avei,avei2,N));
	avei_err=get<1>(BlockStat(avei,avei2,N));
	//variance
	vector<double> vari, vari2, vari_sum, vari_err;
	for(int i=0;i<N;i++){
		double sum=0;
		for(int j=0;j<L;j++){
			int k=j+i*L;
			sum+=pow((M_PI/2)*cos(M_PI*ri[k]/2)/(2*(1-ri[k]))-1,2);
		}
		sum/=L;
		vari.push_back(sum);								//block variance
		vari2.push_back(sum*sum);
	}
	vari_sum=get<0>(BlockStat(vari,vari2,N));
	vari_err=get<1>(BlockStat(vari,vari2,N));
	
	PrintError("ImportanceSampling.txt",avei_sum,avei_err);
	
	rnd.SaveSeed();
	return 0;
}