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
	
	int M=100000;
	int N=100;
	int L=M/N;
	vector<double> r;
	
	for(int i=0;i<M;i++){				//rng
		double a=rnd.Rannyu();
		r.push_back(a);
	}
	
	//average
	vector<double> ave, ave2, ave_sum, ave_error;
	for(int i=0;i<N;i++){
		double sum=0;
		for(int j=0;j<L;j++){
			int k=j+i*L;
			sum+=r[k];
		}
		sum/=L;
		ave.push_back(sum);								//block average
		ave2.push_back(sum*sum);
	}
	ave_sum=get<0>(BlockStat(ave,ave2,N));
	ave_error=get<1>(BlockStat(ave,ave2,N));
	PrintError("average.txt",ave_sum,ave_error);
	
	//variance
	vector<double> var, var2, var_sum, var_error;
	for(int i=0;i<N;i++){
		double sum=0;
		for(int j=0;j<L;j++){
			int k=j+i*L;
			sum+=pow(r[k]-0.5,2);
		}
		sum/=L;
		var.push_back(sum);								//block variance
		var2.push_back(sum*sum);
	}
	var_sum=get<0>(BlockStat(var,var2,N));
	var_error=get<1>(BlockStat(var,var2,N));
	PrintError("variance.txt",var_sum,var_error);
	
	
	//chi2 
	int m=100;
	int n=10000;
	double div=n/m;
	vector<int> chi2;
	for(int i=0;i<m;i++){
		int bin[100]={0};								//array of counters
		double x2=0;
		for(int j=0;j<n;j++){
			int x=rnd.Rannyu(0,100);					//0->1 to 0->100, increase the relative counter
			bin[x]++;
		}
		for(int k=0;k<100;k++){
			x2+=pow(bin[k]-div,2)/div;
		}
		chi2.push_back(x2);
	}
	Print("chi2.txt",chi2);
	
	return 0;
}