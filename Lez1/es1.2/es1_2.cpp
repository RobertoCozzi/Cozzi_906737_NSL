#include "random.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
using namespace std;

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
	
	int N=10000;
	int times[4]={1,2,10,100};
	
	//standard dice
	double Unif[N][4];
	ofstream outunif("Uniform.txt");
	for(int i=0;i<N;i++){
		for(int j=0;j<4;j++){
			double x=0;
			for(int k=0;k<times[j];k++){
				x+=int(rnd.Rannyu(1,7));
			}
			x/=times[j];
			Unif[i][j]=x;
		}
		outunif<<i+1<<" "<<Unif[i][0]<<" "<<Unif[i][1]<<" "<<Unif[i][2]<<" "<<Unif[i][3]<<endl;
	}
	outunif.close();
	
	//exponential dice
	double Exp[N][4];
	ofstream outexp("Exponential.txt");
	for(int i=0;i<N;i++){
		for(int j=0;j<4;j++){
			double x=0;
			for(int k=0;k<times[j];k++){
				x+=rnd.Exponential(1.);
			}
			x/=times[j];
			Exp[i][j]=x;
		}
		outexp<<i+1<<" "<<Exp[i][0]<<" "<<Exp[i][1]<<" "<<Exp[i][2]<<" "<<Exp[i][3]<<endl;
	}
	outexp.close();
	
	//Cauchy dice
	double Cau[N][4];
	ofstream outcau("Cauchy.txt");
	for(int i=0;i<N;i++){
		for(int j=0;j<4;j++){
			double x=0;
			for(int k=0;k<times[j];k++){
				x+=rnd.Cauchy(0,1);
			}
			x/=times[j];
			Cau[i][j]=x;
		}
		outcau<<i+1<<" "<<Cau[i][0]<<" "<<Cau[i][1]<<" "<<Cau[i][2]<<" "<<Cau[i][3]<<endl;
	}
	outcau.close();
	
	return 0;
}