#include "population.h"
#include "individual.h"
#include "random.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

population::population(double perm, double block_perm, double shift, double invert, double cross, int pop, vector<double> x, vector<double> y, Random rnd){
	//mutation rates
	m_perm=perm;
	m_block_perm=block_perm;
	m_shift=shift;
	m_invert=invert;
	m_cross=cross;
	
	m_rnd=rnd;
	m_pop=pop;
	m_ncities=x.size();
	m_x=x;
	m_y=y;
	m_gen=0;
	first_gen();
	Sort();
}

void population::first_gen(){
	for(int i=0;i<m_pop;i++){
		vector<int> path=random_path();
		individual RP(m_x, m_y, path);
		m_population.push_back(RP);
	}
}

vector<int> population::random_path(){
	vector<int> path;
	for(int i=0;i<m_ncities;i++){
		path.push_back(i);
	}
	for(int i=0;i<50;i++){
		int pos1=0, pos2=0, p=0;
		pos1=(int)m_rnd.Rannyu(1,m_ncities);
		pos2=(int)m_rnd.Rannyu(1,m_ncities);
		p=path[pos1];
		path[pos1]=path[pos2];
		path[pos2]=p;
	}
	return path;
}

void population::next_gen(){
	Sort();
	vector<individual> next;
	for(int i=0;i<m_pop;i++){
		int pos=select();
		individual chosen=m_population[pos];
		chosen.set_norm();
		next.push_back(chosen);
	}
	m_population=next;
	int k=(int)m_rnd.Rannyu(1,m_pop/5);
	for(int i=0;i<k;i++){ 
		double cross=m_rnd.Rannyu();
		if(cross<m_cross){
			crossover();
		}
	}
	mutation();
	m_gen++;
	Sort();
}

void population::mutation(){
	int k=(int)m_rnd.Rannyu(1,m_pop);
	for(int i=0;i<k;i++){
		double perm=m_rnd.Rannyu();
		if(perm<m_perm){
			int pos1=0, pos2=0, target=0;
			pos1=(int)m_rnd.Rannyu(1,m_ncities);
			pos2=(int)m_rnd.Rannyu(1,m_ncities);
			target=(int)m_rnd.Rannyu(0,m_pop);
			if(pos1==pos2) pos2+=1;
			m_population[target].permutation(pos1,pos2);
			m_population[target].set_norm();
		}
	}
	
	k=(int)m_rnd.Rannyu(1,m_pop/2);
	for(int i=0;i<k;i++){
		double block_perm=m_rnd.Rannyu();
		if(block_perm<m_block_perm){
			int start1=0, start2=0, length=0, target=0;
			start1=(int)m_rnd.Rannyu(1,m_ncities/2+1);
			start2=(int)m_rnd.Rannyu(m_ncities/2+1,m_ncities);
			length=(int)m_rnd.Rannyu(1,m_ncities/2);
			target=(int)m_rnd.Rannyu(0,m_ncities);
			m_population[target].block_permutation(start1,start2,length);
			m_population[target].set_norm();
		}
	}

	k=(int)m_rnd.Rannyu(1,m_pop/5);
	for(int i=0;i<k;i++){
		double shift=m_rnd.Rannyu();
		if(shift<m_shift){
			int start=0, length=0, target=0;
			start=(int)m_rnd.Rannyu(1,m_ncities-1);
			length=(int)m_rnd.Rannyu(1,m_ncities/2);
			target=(int)m_rnd.Rannyu(0,m_pop);
			m_population[target].shift(start,length);
			m_population[target].set_norm();
		}
	}
	
	k=(int)m_rnd.Rannyu(1,m_pop/5);
	for(int i=0;i<k;i++){
		double invert=m_rnd.Rannyu();
		if(invert<m_invert){
			int start=0, end=0, target=0;
			start=(int)m_rnd.Rannyu(1,m_ncities/2);
			end=(int)m_rnd.Rannyu(m_ncities/2,m_ncities);
			target=(int)m_rnd.Rannyu(0,m_pop);
			m_population[target].inversion(start,end);
			m_population[target].set_norm();
		}
	}
}

void population::crossover(){
	int pos1=select();
	int pos2=select();
	vector<int> parent1=m_population[pos1].get_path();
	vector<int> parent2=m_population[pos2].get_path();
	int start=int(m_rnd.Rannyu(1,m_ncities));
	vector<int> child1, child2;
	for(int i=0;i<start;i++){
		child1.push_back(parent1[i]);
		child2.push_back(parent2[i]);
	}
	for(int i=start;i<m_ncities;i++){
		for(int j=0;j<m_ncities;j++){
			if(parent1[i]==parent2[j]) child1.push_back(parent2[j]);
			if(parent2[i]==parent1[j]) child2.push_back(parent1[j]);
		}
	}
	individual Child1(m_x, m_y, child1);
	individual Child2(m_x, m_y, child2);
	Child1.set_norm();
	Child2.set_norm();
	m_population[pos1]=Child1;
	m_population[pos2]=Child2;
}

int population::select(){
	Sort();
	double r=m_rnd.Rannyu();
	return int(m_pop*pow(r,2));
}

void population::statistics(const char* filename){
	Sort();
	double ave=0, stddev=0;
	ofstream out_stat(filename,ios::app);
	for(int i=0;i<m_pop/2;i++){
		ave+=m_population[i].get_norm();
	}
	ave/=(m_pop/2);
	for(int i=0;i<m_pop/2;i++){
		stddev+=(m_population[i].get_norm()-ave)*(m_population[i].get_norm()-ave);
	}
	stddev=sqrt(stddev/(m_pop/2-1));
	out_stat<<m_gen<<" "<<ave<<" "<<stddev<<endl;
	out_stat.close();
}

void population::print_best(const char* filename){
	Sort();
	cout<<"Best ";
	m_population[0].print_path();
	vector<int> best=m_population[0].get_path();
	ofstream out_best(filename,ios::app);
	for(int i=0;i<m_ncities;i++){
		int pos=best[i];
		out_best<<i<<" "<<m_x[pos]<<" "<<m_y[pos]<<endl;
	}
	int pos=best[0];
	out_best<<m_ncities+1<<" "<<m_x[pos]<<" "<<m_y[pos]<<endl;
	out_best.close();
}