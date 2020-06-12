#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <tuple>
using namespace std;

template <typename T> double CalcMean (const vector <T> &);
template <typename T> double CalcVariance (const vector <T> &);
template <typename T> T CalcMedian (vector <T> &);
template <typename T> tuple<vector<double>, vector<double> > BlockStat (const vector <T> &, const vector <T> &, int);
template <typename T> vector<T> ReadAll (const char*);
template <typename T> void Print (const vector <T> &);
template <typename T> void Print (const char*,const vector <T> &);
template <typename T> void PrintError (const char*, const vector <T> &, const vector <T> &);

template <typename T> vector<T> ReadAll(const char* filename){
	vector<T> v;
	ifstream fin(filename);
	assert(fin && "Input file does not exist");
	if (!fin){
		cerr<<"Cannot open file "<<filename<<endl;
		exit(11);
	} 
	else{
		while(!fin.eof()){
			T val=0;
			fin>>val;
			v.push_back(val);
		};
	}
	return v;
};

template <typename T> double CalcMean(const vector <T> & v){
	double mean=0;
	for(int k=0;k<v.size();k++){
		mean+=v[k];
	}
	mean/=double(v.size());
	return mean;
};

template <typename T> double CalcVariance(const vector <T> & v){
	double mean=CalcMean(v);
	double variance=0;
	for(int k=0;k<v.size();k++){
		variance+=((v[k]-mean)*(v[k]-mean));
	}
	variance/=double(v.size());
	return variance;
};

template <typename T> T CalcMedian(vector<T> & v){
	int a=v.size();
	sort(v.begin(),v.end());
	T median=0;
	if(a%2==0){
		median=(v[a/2-1]+v[a/2])/2;
	}
	else{
		median=v[a/2];
	}
	return median;
};

template <typename T> tuple<vector<double>, vector<double>> BlockStat(const vector <T> & mean, const vector <T> & mean2, int blocks){
	if(blocks!=mean.size() or blocks!=mean2.size()){
		cerr<<"Imput vectors of wrong dimension"<<endl;
		exit(11);
	}
	else{
		vector<double> mean_sum, mean2_sum, mean_err;
		for(int k=0;k<blocks;k++){
			double sum=0;
			double sum2=0;
			for(int j=0;j<k+1;j++){
				sum+=mean[j];
				sum2+=mean2[j];
			}
			sum/=(k+1);
			sum2/=(k+1);
			mean_sum.push_back(sum);
			mean2_sum.push_back(sum2);
			if(k==0) mean_err.push_back(0);
			else{
				double err=sqrt((sum2-sum*sum)/k);
				mean_err.push_back(err);
			}
		}
		return make_tuple(mean_sum, mean_err);
	}
};

template <typename T> void Print(const char*  filename, const vector<T> & v){		//print on file
	ofstream fout(filename);
	if(!fout){
		cerr<<"Cannot create the file "<<filename<<endl;
		return;
	}
	for(int k=0;k<v.size();k++) fout<<k+1<<" "<<v[k]<<endl;
	fout.close();
	return;
};

template <typename T> void Print(const vector<T> & v){								//print on terminal
	for(int k=0;k<v.size();k++) cout<<v[k]<<endl;
};

template <typename T> void PrintError(const char* filename, const vector<T> & v, const vector<T> & err){
	ofstream fout(filename);
	if(!fout){
		cerr<<"Cannot create the file "<<filename<<endl;
		return;
	}
	for(int k=0;k<v.size();k++) fout<<k+1<<" "<<v[k]<<" "<<err[k]<<endl;
	fout.close();
	return;
};