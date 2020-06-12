#ifndef __population_h__
#define __population_h__

#include "individual.h"
#include "random.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

class population{
	public:
		population(){};
		population(double perm, double block_perm, double shift, double invert, double beta, double cool, int cool_step, vector<double> x, vector<double> y, Random rnd);
		~population(){};
		
		double get_perm() const {return m_perm;}
		double get_block_perm() const {return m_block_perm;}
		double get_shift() const {return m_shift;}
		double get_invert() const {return m_invert;}
		double get_beta() const {return m_beta;}
		double get_cool() const {return m_cool;}
		int get_ncities() const {return m_ncities;}
		int get_step() const {return m_cool_step;}
		Random get_rnd() const {return m_rnd;}
		
		void first_gen();
		vector<int> random_path();
		void next_gen();
		void mutation();
		double Boltzmann(double norm);
		
		void statistics(const char* filename);
		void print_best(const char* filename);
	
	protected:
		double m_perm, m_block_perm, m_shift, m_invert;
		double m_beta, m_cool;
		int m_cool_step, m_ncities, m_gen, m_accept;
		individual m_path_new, m_path_old;
		vector<double> m_x, m_y;
		Random m_rnd;
};

#endif