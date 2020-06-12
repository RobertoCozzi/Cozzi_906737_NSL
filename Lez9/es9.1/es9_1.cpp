#include "random.h"
#include "individual.h"
#include "population.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
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
	
	//parameters
	double perm=0.1, block_perm=0.07, shift=0.05, inv=0.05, cross=0.65;		//mutation rates
	double radius=1, side=1;
	int cities=32, pop=200, gens=2000;
	
	cout<<"Traveling Salesman Problem using a genetic algorithm"<<endl;
	cout<<"cities: "<<cities<<", paths: "<<pop<<", generations: "<<gens<<endl;
	cout<<"mutation rates:"<<endl;
	cout<<"     permutation       "<<perm<<endl;
	cout<<"     block permutation "<<block_perm<<endl;
	cout<<"     shift             "<<shift<<endl;
	cout<<"     inversion         "<<inv<<endl;
	cout<<"     crossover         "<<cross<<endl;
	cout<<"using L2 norm"<<endl<<endl;
	
	cout<<"Cities on a circumference of radius "<<radius<<endl;
	//build cities
	vector<double> x_cities_c, y_cities_c;
	for(int c=0;c<cities;c++){
		double theta=rnd.Rannyu(0,2*M_PI);
		double x=radius*cos(theta);
		double y=radius*sin(theta);
		x_cities_c.push_back(x);
		y_cities_c.push_back(y);
	}
	//build cities
	vector<double> x_cities_s, y_cities_s;
	for(int c=0;c<cities;c++){
		double x=rnd.Rannyu(-side,side);
		double y=rnd.Rannyu(-side,side);
		x_cities_s.push_back(x);
		y_cities_s.push_back(y);
	}
	//create population
	population TSP_c(perm, block_perm, shift, inv, cross, pop, x_cities_c, y_cities_c, rnd);
	cout<<"first generation ready"<<endl;
	//evolution
	for(int g=0;g<gens;g++){
		TSP_c.next_gen();
		TSP_c.statistics("TSP_L2_cir_ave.txt");
		if(g%200==0) cout<<"generation "<<g<<endl;
	}
	TSP_c.print_best("TSP_L2_cir_best.txt");
	cout<<endl;
	
	cout<<"Cities inside a square of side "<<2*side<<endl;
	//create population
	population TSP_s(perm, block_perm, shift, inv, cross, pop, x_cities_s, y_cities_s, rnd);
	cout<<"first generation ready"<<endl;
	//evolution
	for(int g=0;g<gens;g++){
		TSP_s.next_gen();
		TSP_s.statistics("TSP_L2_sq_ave.txt");
		if(g%200==0) cout<<"generation "<<g<<endl;
	}
	TSP_s.print_best("TSP_L2_sq_best.txt");
	cout<<endl;
	
	return 0;
}