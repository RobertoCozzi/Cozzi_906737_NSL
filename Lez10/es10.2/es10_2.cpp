#include "random.h"
#include "individual.h"
#include "population.h"
#include "mpi.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

int main(int argc, char* argv[]){
	
	MPI_Init(&argc,&argv);
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	Random rnd;
	int seed[4];
	int p[size][2];
	ifstream Primes("Primes");
	if (Primes.is_open()){
		for(int i=0;i<size;i++){
			Primes >> p[i][0] >> p[i][1];
		}
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p[rank][0],p[rank][1]);
			}
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	
	//parameters
	double perm=0.15, block_perm=0.10, shift=0.09, inv=0.09, cross=0.75;		//mutation rates
	double radius=1, side=1;
	int cities=32, pop=200, gens=2000, mig=100;
	
	if(rank==0){
		cout<<"Traveling Salesman Problem using a genetic algorithm"<<endl;
		cout<<"cities: "<<cities<<", paths: "<<pop<<", generations: "<<gens<<endl;
		cout<<"mutation rates:"<<endl;
		cout<<"     permutation       "<<perm<<endl;
		cout<<"     block permutation "<<block_perm<<endl;
		cout<<"     shift             "<<shift<<endl;
		cout<<"     inversion         "<<inv<<endl;
		cout<<"     crossover         "<<cross<<endl;
		cout<<"using L2 norm"<<endl<<endl;
	}
	
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
	
	if(rank==0) cout<<"Cities on a circumference of radius "<<radius<<endl;
	//create population
	population TSP_c(perm, block_perm, shift, inv, cross, pop, x_cities_c, y_cities_c, rnd);
	//evolution
	for(int g=0;g<gens;g++){
		TSP_c.next_gen();
		if(g!=0 and g%mig==0){
			MPI_Status stat1, stat2, stat3, stat4;
			MPI_Request req1, req2;
			int root=0;
			int itag1=1, itag2=2, itag3=3, itag4=4;
			int target[3];
			if(rank==0){
				target[0]=(int)rnd.Rannyu(1,4);
				do{
					target[1]=(int)rnd.Rannyu(1,4);
				}while(target[1]==0 or target[1]==target[0]);
				do{
					target[2]=(int)rnd.Rannyu(1,4);
				}while(target[2]==0 or target[2]==target[0] or target[2]==target[1]);
			}
			MPI_Bcast(target, 3, MPI_INT, root, MPI_COMM_WORLD);
			
			vector<int> path_exc=TSP_c.get_path(0);
			vector<int> new_path(cities);
			int* path1=new int[cities];
			int* path2=new int[cities];
			for(int c=0;c<cities;c++){
				path1[c]=path_exc[c];
			}
			MPI_Barrier(MPI_COMM_WORLD);
			
			if(rank==0){
				MPI_Isend(&path1[0], cities, MPI_INTEGER, target[0], itag1, MPI_COMM_WORLD, &req1);
				MPI_Recv(&path2[0], cities, MPI_INTEGER, target[0], itag2, MPI_COMM_WORLD, &stat2);
			}
			else if(rank==target[0]){
				MPI_Send(&path1[0], cities, MPI_INTEGER, 0, itag2, MPI_COMM_WORLD);
				MPI_Recv(&path2[0], cities, MPI_INTEGER, 0, itag1, MPI_COMM_WORLD, &stat1);
			}
			else if(rank==target[1]){
				MPI_Isend(&path1[0], cities, MPI_INTEGER, target[2], itag3, MPI_COMM_WORLD, &req2);
				MPI_Recv(&path2[0], cities, MPI_INTEGER, target[2], itag4, MPI_COMM_WORLD, &stat4);
			}
			else if(rank==target[2]){
				MPI_Send(&path1[0], cities, MPI_INTEGER, target[1], itag4, MPI_COMM_WORLD);
				MPI_Recv(&path2[0], cities, MPI_INTEGER, target[1], itag3, MPI_COMM_WORLD, &stat3);
			}
			for(int c=0;c<cities;c++){
				new_path[c]=path2[c];
			}
			
			TSP_c.set_path(0, new_path);
			MPI_Barrier(MPI_COMM_WORLD);
			if(rank==0) cout<<" generation "<<g<<" migration completed"<<endl;
		}
		TSP_c.statistics(to_string(rank)+"_TSP_L2_cir_ave.txt");
	}
	TSP_c.print_best(to_string(rank)+"_TSP_L2_cir_best.txt" );
	cout<<endl;
	
	if(rank==0) cout<<"Cities inside a square of side "<<2*side<<endl;
	//create population
	population TSP_s(perm, block_perm, shift, inv, cross, pop, x_cities_s, y_cities_s, rnd);
	//evolution
	for(int g=0;g<gens;g++){
		TSP_s.next_gen();
		if(g!=0 and g%mig==0){
			MPI_Status stat1, stat2, stat3, stat4;
			MPI_Request req1, req2;
			int root=0;
			int itag1=1, itag2=2, itag3=3, itag4=4;
			int target[3];
			if(rank==0){
				target[0]=(int)rnd.Rannyu(1,4);
				do{
					target[1]=(int)rnd.Rannyu(1,4);
				}while(target[1]==0 or target[1]==target[0]);
				do{
					target[2]=(int)rnd.Rannyu(1,4);
				}while(target[2]==0 or target[2]==target[0] or target[2]==target[1]);
			}
			MPI_Bcast(target, 3, MPI_INT, root, MPI_COMM_WORLD);
			
			vector<int> path_exc=TSP_s.get_path(0);
			vector<int> new_path(cities);
			int* path1=new int[cities];
			int* path2=new int[cities];
			for(int c=0;c<cities;c++){
				path1[c]=path_exc[c];
			}
			MPI_Barrier(MPI_COMM_WORLD);
			
			if(rank==0){
				MPI_Isend(&path1[0], cities, MPI_INTEGER, target[0], itag1, MPI_COMM_WORLD, &req1);
				MPI_Recv(&path2[0], cities, MPI_INTEGER, target[0], itag2, MPI_COMM_WORLD, &stat2);
			}
			else if(rank==target[0]){
				MPI_Send(&path1[0], cities, MPI_INTEGER, 0, itag2, MPI_COMM_WORLD);
				MPI_Recv(&path2[0], cities, MPI_INTEGER, 0, itag1, MPI_COMM_WORLD, &stat1);
			}
			else if(rank==target[1]){
				MPI_Isend(&path1[0], cities, MPI_INTEGER, target[2], itag3, MPI_COMM_WORLD, &req2);
				MPI_Recv(&path2[0], cities, MPI_INTEGER, target[2], itag4, MPI_COMM_WORLD, &stat4);
			}
			else if(rank==target[2]){
				MPI_Send(&path1[0], cities, MPI_INTEGER, target[1], itag4, MPI_COMM_WORLD);
				MPI_Recv(&path2[0], cities, MPI_INTEGER, target[1], itag3, MPI_COMM_WORLD, &stat3);
			}
			for(int c=0;c<cities;c++){
				new_path[c]=path2[c];
			}
			
			TSP_s.set_path(0, new_path);
			MPI_Barrier(MPI_COMM_WORLD);
			if(rank==0) cout<<" generation "<<g<<" migration completed"<<endl;
		}
		TSP_s.statistics(to_string(rank)+"_TSP_L2_sq_ave.txt");
	}
	TSP_s.print_best(to_string(rank)+"_TSP_L2_sq_best.txt" );
	cout<<endl;
	
	
	MPI_Finalize();
	return 0;
}