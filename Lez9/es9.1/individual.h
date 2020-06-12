#ifndef __individual_h__
#define __individual_h__

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

class individual{
	public:
		individual(){};
		individual(vector<double> x, vector<double>y, vector<int> path);
		~individual(){};
		
		void set_norm();
		double get_norm() const {return m_norm;}
		double get_x(int pos) const {return m_x[pos];}
		double get_y(int pos) const {return m_y[pos];}
		int get_size() const {return m_ncities;}
		vector<int> get_path() const {return m_path;}
		
		void print_path() const;
		void check();
		
		//Mutations
		void permutation(int pos1, int pos2);
		void block_permutation(int start1, int start2, int length);
		void shift(int start, int length);
		void inversion(int start, int end);
	
	protected:
		int m_ncities;
		double m_norm;
		vector<int> m_path;
		vector<double> m_x, m_y;
		int Pbc(int i);
};

bool operator<(const individual& path1, const individual& path2);

#endif