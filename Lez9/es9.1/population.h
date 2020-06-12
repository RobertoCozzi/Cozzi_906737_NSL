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
		population(double perm, double block_perm, double shift, double invert, double cross, int pop, vector<double> x, vector<double> y, Random rnd);
		~population(){};
		
		double get_perm() const {return m_perm;}
		double get_block_perm() const {return m_block_perm;}
		double get_shift() const {return m_shift;}
		double get_invert() const {return m_invert;}
		double get_cross() const {return m_cross;}
		int get_ncities() const {return m_ncities;}
		int get_pop() const {return m_pop;}
		
		void Sort() {sort(m_population.begin(),m_population.end());}
		
		void first_gen();
		vector<int> random_path();
		void next_gen();
		void mutation();
		void crossover();
		int select();
		
		void statistics(const char* filename);
		void print_best(const char* filename);
	
	protected:
		double m_perm, m_block_perm, m_shift, m_invert, m_cross;
		int m_pop, m_ncities, m_gen;
		vector<individual> m_population;
		vector<double> m_x, m_y;
		Random m_rnd;
};

#endif