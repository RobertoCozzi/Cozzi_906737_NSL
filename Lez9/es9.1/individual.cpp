#include "individual.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

individual::individual(vector<double> x, vector<double> y, vector<int> path){
	m_x=x;
	m_y=y;
	m_path=path;
	m_ncities=m_path.size();
	set_norm();
}

void individual::set_norm(){
	double norm=0;
	for(int i=0;i<m_ncities;i++){
		int start=m_path[i];
		int end=m_path[Pbc(i+1)];
		double dx=m_x[start]-m_x[end];
		double dy=m_y[start]-m_y[end];
		norm+=sqrt(dx*dx+dy*dy);
	}
	
	m_norm=norm;
}

void individual::print_path() const{
	cout<<"Path: ";
	for(int i=0;i<m_ncities;i++) cout<<m_path[i]<<" ";
	cout<<"length "<<m_norm<<endl;
}

void individual::permutation(int pos1, int pos2){
	pos1=Pbc(pos1);
	pos2=Pbc(pos2);
	int p=m_path[pos1];
	m_path[pos1]=m_path[pos2];
	m_path[pos2]=p;
	check();
}

void individual::block_permutation(int start1, int start2, int length){
	for(int i=0;i<length;i++){
		permutation(Pbc(start1+i),Pbc(start2+i));
	}
	check();
}

void individual::shift(int start, int length){
	for(int i=0;i<length;i++){
		permutation(Pbc(start+i),Pbc(start+length+i));
	}
	check();
}

void individual::inversion(int start, int end){
	int length=floor((start-end)/2);
	for(int i=0;i<length;i++){
		permutation(Pbc(start+i),Pbc(end-i));
	}
	check();
}

void individual::check(){
	int sum=0;
	for(int i=0;i<m_ncities;i++){
		for(int j=i+1;j<m_ncities;j++){
			if(m_path[i]==m_path[j]){
				cout<<"city "<<m_path[i]<<" visited more than once"<<endl;
				sum++;
			}
		}
	}
	if(sum!=0) print_path();
}

int individual::Pbc(int i){
	if(i>=m_ncities) i=i-m_ncities;
	return i;
}

bool operator<(const individual& path1, const individual& path2){
	return path1.get_norm()<path2.get_norm();
}