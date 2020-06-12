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
	
	double S0=100., t=0., T=1., K=100., r=0.1, vol=0.25;
	int M=10000, N=100;
	double L=M/N;
	
	//final price
	vector<double> Cf, Pf, Call_final, Call_final2, Put_final, Put_final2, Call_final_sum, Put_final_sum, Call_final_err, Put_final_err;
	for(int i=0;i<M;i++){
		double W=rnd.Gauss(0,T);
		double vol2=pow(vol,2.)/2;
		double s_f=S0*exp((r-vol2)*T+vol*W);								//final price
		Cf.push_back(max(0.,s_f-K)*exp(-r*T));
		Pf.push_back(max(0.,K-s_f)*exp(-r*T));
	}
	for(int i=0;i<N;i++){
		double sumC=0;
		double sumP=0;
		for(int j=0;j<L;j++){
			int k=j+i*L;
			sumC+=Cf[k];
			sumP+=Pf[k];
		}
		sumC/=L;
		sumP/=L;
		Call_final.push_back(sumC);
		Call_final2.push_back(sumC*sumC);
		Put_final.push_back(sumP);
		Put_final2.push_back(sumP*sumP);
	}
	Call_final_sum=get<0>(BlockStat(Call_final,Call_final2,N));
	Call_final_err=get<1>(BlockStat(Call_final,Call_final2,N));
	Put_final_sum=get<0>(BlockStat(Put_final,Put_final2,N));
	Put_final_err=get<1>(BlockStat(Put_final,Put_final2,N));
	
	PrintError("Callfinal.txt",Call_final_sum,Call_final_err);
	PrintError("Putfinal.txt",Put_final_sum,Put_final_err);
	
	//discrete sampling
	vector<double> Cd, Pd, Call_discrete, Call_discrete2, Put_discrete, Put_discrete2, Call_discrete_sum, Put_discrete_sum, Call_discrete_err, Put_discrete_err;
	for(int i=0;i<M;i++){
		double dt=(T-t)/100;
		double Si=S0, Si_1=0;
		double vol2=pow(vol,2)/2;
		for(int j=1;j<100;j++){
			double Z=rnd.Gauss(0,1);
			double Si_1=Si*exp((r-vol2)*dt+vol*Z*sqrt(dt));
			Si=Si_1;
		}
		Cd.push_back(max(0.,Si-K)*exp(-r*T));
		Pd.push_back(max(0.,K-Si)*exp(-r*T));
	}
	for(int i=0;i<N;i++){
		double sumC=0;
		double sumP=0;
		for(int j=0;j<L;j++){
			int k=j+i*L;
			sumC+=Cd[k];
			sumP+=Pd[k];
		}
		sumC/=L;
		sumP/=L;
		Call_discrete.push_back(sumC);
		Call_discrete2.push_back(sumC*sumC);
		Put_discrete.push_back(sumP);
		Put_discrete2.push_back(sumP*sumP);
	}
	Call_discrete_sum=get<0>(BlockStat(Call_discrete,Call_discrete2,N));
	Call_discrete_err=get<1>(BlockStat(Call_discrete,Call_discrete2,N));
	Put_discrete_sum=get<0>(BlockStat(Put_discrete,Put_discrete2,N));
	Put_discrete_err=get<1>(BlockStat(Put_discrete,Put_discrete2,N));
	
	PrintError("Calldiscrete.txt",Call_discrete_sum,Call_discrete_err);
	PrintError("Putdiscrete.txt",Put_discrete_sum,Put_discrete_err);
	
	rnd.SaveSeed();
	return 0;
}