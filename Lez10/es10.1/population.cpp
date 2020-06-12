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

population::population(double perm, double block_perm, double shift, double invert, double beta, double cool, int cool_step, vector<double> x, vector<double> y, Random rnd){
	//mutation rates
	m_perm=perm;
	m_block_perm=block_perm;
	m_shift=shift;
	m_invert=invert;
	
	m_beta=beta;
	m_cool=cool;
	m_cool_step=cool_step;
	
	m_rnd=rnd;
	m_ncities=x.size();
	m_x=x;
	m_y=y;
	m_gen=0;
	m_accept=0;
	first_gen();
}

void population::first_gen(){
	vector<int> path=random_path();
	individual RP(m_x, m_y, path);
	m_path_old=RP;
	m_path_new=RP;
}

vector<int> population::random_path(){
	vector<int> path;
	for(int i=0;i<m_ncities;i++){
		path.push_back(i);
	}
	for(int i=0;i<100;i++){
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
	m_path_new=m_path_old;
	mutation();
	double p_old=Boltzmann(m_path_old.get_norm());
	double p_new=Boltzmann(m_path_new.get_norm());
	double p=p_new/p_old;
	double r=m_rnd.Rannyu();
	if(r<=p){
		m_path_old=m_path_new;
		m_accept++;
	}
	m_gen++;
	if(m_gen%m_cool_step==0) m_beta/=m_cool;
}

void population::mutation(){
	for(int i=0;i<2;i++){
		double perm=m_rnd.Rannyu();
		if(perm<m_perm){
			int pos1=0, pos2=0;
			pos1=(int)m_rnd.Rannyu(1,m_ncities);
			pos2=(int)m_rnd.Rannyu(1,m_ncities);
			if(pos1==pos2) pos2+=8;
			m_path_new.permutation(pos1,pos2);
			m_path_new.set_norm();
		}
	}

	double block_perm=m_rnd.Rannyu();
	if(block_perm<m_block_perm){
		int start1=0, start2=0, length=0;
		start1=(int)m_rnd.Rannyu(1,m_ncities/2+1);
		start2=(int)m_rnd.Rannyu(m_ncities/2+1,m_ncities);
		length=(int)m_rnd.Rannyu(1,m_ncities/2);
		m_path_new.block_permutation(start1,start2,length);
		m_path_new.set_norm();
	}

	double shift=m_rnd.Rannyu();
	if(shift<m_shift){
		int start=0, length=0;
		start=(int)m_rnd.Rannyu(1,m_ncities-1);
		length=(int)m_rnd.Rannyu(1,m_ncities/2);
		m_path_new.shift(start,length);
		m_path_new.set_norm();
	}

	double inv=m_rnd.Rannyu();
	if(inv<m_invert){
		int start=0, end=0;
		start=(int)m_rnd.Rannyu(1,m_ncities/2);
		end=(int)m_rnd.Rannyu(m_ncities/2,m_ncities);
		m_path_new.inversion(start,end);
		m_path_new.set_norm();
	}

}

double population::Boltzmann(double norm){
	return exp(-m_beta*norm);
}

void population::statistics(const char* filename){
	ofstream out_stat(filename,ios::app);
	out_stat<<m_gen<<" "<<m_path_old.get_norm()<<" "<<m_beta<<endl;
	out_stat.close();
}

void population::print_best(const char* filename){
	cout<<"Best ";
	m_path_old.print_path();
	cout<<"Accepted moves "<<m_accept<<" of "<<m_gen<<endl;
	vector<int> best=m_path_old.get_path();
	ofstream out_best(filename,ios::app);
	for(int i=0;i<m_ncities;i++){
		int pos=best[i];
		out_best<<i<<" "<<m_x[pos]<<" "<<m_y[pos]<<endl;
	}
	int pos=best[0];
	out_best<<m_ncities+1<<" "<<m_x[pos]<<" "<<m_y[pos]<<endl;
	out_best.close();
}