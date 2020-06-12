/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_NVT.h"

using namespace std;

int main(int argc, char** argv){
	if (argc<1){
		cout<<"missing input file"<<endl;
		return -1;
	}
	char* input_file = argv[1];
	
	Input(input_file); //Inizialization
	
	int i_p=instant_print;
	instant_print=0;
	for(int i=0;i<eq_step;i++){			//equilibration
		Move();
	}
	
	instant_print=i_p;
	
	int nconf = 1;
	for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
		if(iblk%10==0) cout<<"block "<<iblk<<endl;
		Reset(iblk);   //Reset block averages
		for(int istep=1; istep <= nstep; ++istep){
			Move();
			Measure();
			Accumulate(); //Update block averages
			if(istep%10 == 0 && print_config==1){
				ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
				nconf += 1;
			}
		}	
		Averages(iblk);   //Print results for current block
	}
	ConfFinal(); //Write final configuration

	return 0;
}


void Input(const char* input_file){
	ifstream ReadInput,ReadConf;

	cout << "Classic Lennard-Jones fluid        " << endl;
	cout << "Monte Carlo simulation             " << endl << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
	cout << "The program uses Lennard-Jones units " << endl;

//Read seed for random numbers
	int p1, p2;
	ifstream Primes("Primes");
	Primes >> p1 >> p2 ;
	Primes.close();

	ifstream input("seed.in");
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
	input.close();
  
//Read input informations
	ReadInput.open(input_file);

	ReadInput >> temp;
	beta = 1.0/temp;
	cout << "Temperature = " << temp << endl;

	ReadInput >> npart;
	cout << "Number of particles = " << npart << endl;

	ReadInput >> rho;
	cout << "Density of particles = " << rho << endl;
	vol = (double)npart/rho;
	box = pow(vol,1.0/3.0);
	cout << "Volume of the simulation box = " << vol << endl;
	cout << "Edge of the simulation box = " << box << endl;

	ReadInput >> rcut;
	cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
	//Tail corrections for potential energy and pressure
	vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
	ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
	cout << "Tail correction for the potential energy = " << vtail << endl;
	cout << "Tail correction for the virial           = " << ptail << endl; 

	ReadInput >> delta;
	ReadInput >> eq_step;
	ReadInput >> nblk;
	ReadInput >> nstep;
	ReadInput >> instant_print;
	ReadInput >> print_config;

	cout << "The program perform Metropolis moves with uniform translations" << endl;
	cout << "Moves parameter = " << delta << endl;
	cout << "Number of blocks = " << nblk << endl;
	cout << "Number of steps in one block = " << nstep << endl << endl;
	ReadInput.close();


//Prepare arrays for measurements
	iv = 0; //Potential energy
	iw = 1; //Virial
 
	n_props = 2; //Number of observables

//measurement of g(r)
	igofr = 2;
	nbins = 100;
	n_props = n_props + nbins;
	bin_size = (box/2.0)/(double)nbins;

//Read initial configuration
	cout << "Read initial configuration from file config.0 " << endl << endl;
	ReadConf.open("config.0");
	for (int i=0; i<npart; ++i){
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = Pbc( x[i] * box );
		y[i] = Pbc( y[i] * box );
		z[i] = Pbc( z[i] * box );
	}
	ReadConf.close();
  
//Evaluate potential energy and virial of the initial configuration
	Measure();

//Print initial values for the potential energy and virial
	cout << "Initial potential energy (with tail corrections) = " << walker[iv]/(double)npart + vtail << endl;
	cout << "Virial                   (with tail corrections) = " << walker[iw]/(double)npart + ptail << endl;
	cout << "Pressure                 (with tail corrections) = " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl << endl;
}


void Move(void){
	int o;
	double p, energy_old, energy_new;
	double xold, yold, zold, xnew, ynew, znew;

	for(int i=0; i<npart; ++i){
  //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
		o = (int)(rnd.Rannyu()*npart);

		//Old
		xold = x[o];
		yold = y[o];
		zold = z[o];
		energy_old = Boltzmann(xold,yold,zold,o);

		//New
		xnew = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
		ynew = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
		znew = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

		energy_new = Boltzmann(xnew,ynew,znew,o);

		//Metropolis test
		p = exp(beta*(energy_old-energy_new));
		if(p >= rnd.Rannyu()){  
			//Update
			x[o] = xnew;
			y[o] = ynew;
			z[o] = znew;
    
			accepted = accepted + 1.0;
		}
		attempted = attempted + 1.0;
	}
}

double Boltzmann(double xx, double yy, double zz, int ip){
	double ene=0.0;
	double dx, dy, dz, dr;

	for (int i=0; i<npart; ++i){
		if(i != ip){
// distance ip-i in pbc
			dx = Pbc(xx - x[i]);
			dy = Pbc(yy - y[i]);
			dz = Pbc(zz - z[i]);

			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);

			if(dr < rcut){
				ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
			}
		}
	}

	return 4.0*ene;
}

void Measure(){
	int bin;
	double v = 0.0, w = 0.0;
	double vij, wij;
	double dx, dy, dz, dr;
	
//reset the hystogram of g(r)
	for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

//cycle over pairs of particles
	for (int i=0; i<npart-1; ++i){
		for (int j=i+1; j<npart; ++j){

			dx = Pbc(x[i] - x[j]);
			dy = Pbc(y[i] - y[j]);
			dz = Pbc(z[i] - z[j]);

			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);

//g(r)
			for(bin=0;bin<nbins;bin++){
				if(dr>=(bin*bin_size) && dr<((bin+1)*bin_size)) walker[igofr+bin]+=2;
			}

			if(dr < rcut){
				vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
				wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
// contribution to energy and virial
				v += vij;
				w += wij;
			}
		}          
	}

	walker[iv] = 4.0 * v;
	walker[iw] = 48.0 * w / 3.0;
}


void Reset(int iblk){ //Reset block averages
   if(iblk == 1){
		for(int i=0; i<n_props; ++i){
			glob_av[i] = 0;
			glob_av2[i] = 0;
		}
	}

	for(int i=0; i<n_props; ++i){
		blk_av[i] = 0;
	}
	blk_norm = 0;
	attempted = 0;
	accepted = 0;
}


void Accumulate(void){ //Update block averages

	for(int i=0; i<n_props; ++i){
		blk_av[i] = blk_av[i] + walker[i];
	}
	blk_norm = blk_norm + 1.0;
}


void Averages(int iblk){ //Print results for current block
	double r, gdir;
	ofstream Gofr, Gave, Epot, Pres;
	const int wd=12;
    
    //cout << "Block number " << iblk << endl;
	//cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Epot.open("output.epot.0",ios::app);
    Pres.open("output.pres.0",ios::app);
    Gofr.open("output.gofr.0",ios::app);
    Gave.open("output.gave.0",ios::app);
    
	
    stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_press=Error(glob_av[iw],glob_av2[iw],iblk);

//Potential energy per particle
    Epot<<iblk<<" "<<stima_pot<<" "<<glob_av[iv]/(double)iblk<<" "<<err_pot<<endl;
//Pressure 
    Pres<<iblk<<" "<<stima_pres<<" "<<glob_av[iw]/(double)iblk<<" "<<err_press<<endl;
	
//g(r)
	double stima_gofr[nbins], err_gofr[nbins];
	r=bin_size;
	for(int ibin=0;ibin<nbins;ibin++){
		double dv=4*M_PI/3*(pow(r+bin_size,3)-pow(r,3));
		stima_gofr[ibin]=blk_av[ibin+igofr]/blk_norm/(rho*npart*dv);
		glob_av[ibin+igofr]+=stima_gofr[ibin];
		glob_av2[ibin+igofr]+=stima_gofr[ibin]*stima_gofr[ibin];
		err_gofr[ibin]=Error(glob_av[ibin+igofr],glob_av2[ibin+igofr],iblk);
		Gofr<<iblk<<" "<<ibin<<" "<<stima_gofr[ibin]<<endl;
		r+=bin_size;
	}
	
	if(iblk==nblk){
		r=bin_size;
		for(int bin=0;bin<nbins;bin++){
			Gave<<bin<<" "<<r<<" "<<stima_gofr[bin]<<" "<<glob_av[bin+igofr]/(double)iblk<<" "<<err_gofr[bin]<<endl;
			r+=bin_size;
		}
	}
	
    //cout << "----------------------------" << endl << endl;

    Epot.close();
    Pres.close();
    Gofr.close();
	Gave.close();
}


void ConfFinal(void){
	ofstream WriteConf;

	cout << "Print final configuration to file config.final " << endl << endl;
	WriteConf.open("config.final");
	for (int i=0; i<npart; ++i){
		WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
	}	
	WriteConf.close();

	rnd.SaveSeed();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
	ofstream WriteXYZ;

	WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
	WriteXYZ << npart << endl;
	WriteXYZ << "This is only a comment!" << endl;
	for (int i=0; i<npart; ++i){
		WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
	}
	WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk){
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
