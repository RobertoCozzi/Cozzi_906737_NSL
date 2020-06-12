/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"
#include "random.h"

using namespace std;

int main(int argc, char** argv ){
	if (argc<1){
		cout<<"missing input file"<<endl;
		return -1;
	}
	char* input_file = argv[1];            //Inizialization
	Input(input_file);
	
	//equilibration
	print_eq=1;
	cout<<"Equilibration"<<endl<<endl;
	for(int i=0;i<eqstep;i++){
		Move();
		if(i%10==0) Measure();
		if(i==500) ConfFinal("old.0");
	}
	ConfFinal("old.final");
	print_eq=0;
	
	int nconf = 1;
	for(int iblk=1;iblk<=nblock;iblk++){
		if(iblk%10==0) cout<<"block "<<iblk<<endl;
		Reset(iblk);
		for(int istep=0;istep<nstep;istep++){
			Move();
			if(istep%10==0){
				Measure();
				if(print_config==1) ConfXYZ(nconf);
				nconf++;
				Accumulate();
			}
		}
		Averages(iblk);
	}
	
	ConfFinal("config.final");         //Write final configuration to restart

	return 0;
}


void Input(const char* input_file){ //Prepare all stuff for the simulation
	ifstream ReadInput,ReadConf;
	double ep, ek, pr, et, vir;
	
	cout << "Classic Lennard-Jones fluid        " << endl;
	cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "The program uses Lennard-Jones units " << endl;
	
	//Read seed for random numbers
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

	ReadInput.open(input_file); //Read input
	ReadInput >> temp;
	ReadInput >> npart;
	cout << "Number of particles = " << npart << endl;
	ReadInput >> rho;
	cout << "Density of particles = " << rho << endl;
	vol = (double)npart/rho;
	cout << "Volume of the simulation box = " << vol << endl;
	box = pow(vol,1.0/3.0);
	cout << "Edge of the simulation box = " << box << endl;

	ReadInput >> rcut;
	ReadInput >> dt;
	ReadInput >> nstep;
	ReadInput >> nblock;
	ReadInput >> eqstep;
	ReadInput >> old;
	ReadInput >> print_config;
	ReadInput >> instant_print;

	cout << "The program integrates Newton equations with the Verlet method " << endl;
	cout << "Time step = " << dt << endl;
	cout << "Number of blocks = "<< nblock <<endl; 
	cout << "Number of steps = " << nstep << endl << endl;
	ReadInput.close();

	//Prepare array for measurements
	iv = 0; //Potential energy
	ik = 1; //Kinetic energy
	ie = 2; //Total energy
	it = 3; //Temperature
	iw = 4; //Pressure
	n_props = 5; //Number of observables
	
	//measurement of g(r)
	igofr = 5;
	nbins = 100;
	n_props = n_props + nbins;
	bin_size = (box/2.0)/(double)nbins;
	
	//Tail corrections for potential energy and pressure
	vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
	ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
	cout << "Tail correction for the potential energy = " << vtail << endl;
	cout << "Tail correction for the virial           = " << ptail << endl;

	//Read initial configuration
	cout << "Read initial configuration from file config.0 " << endl << endl;
	ReadConf.open("config.0");
	for (int i=0; i<npart; ++i){
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = x[i] * box;
		y[i] = y[i] * box;
		z[i] = z[i] * box;
	}
	ReadConf.close();
	
	//Prepare initial velocities
	if(old==1){
		cout<<"Read the initial configuration from file old.0"<<endl<<endl;
		ifstream ReadConfOld("old.0");
		for(int i=0;i<npart;i++){
			ReadConfOld>>xold[i]>>yold[i]>>zold[i];
			xold[i]=xold[i]*box;
			yold[i]=yold[i]*box;
			zold[i]=zold[i]*box;
		}
		ReadConfOld.close();
		double sumv2=0, fs, xnew, ynew, znew;
		Move();
		for(int i=0;i<npart;i++) sumv2+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
		sumv2 /= (double)npart;
		fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
		for (int i=0; i<npart; ++i){
			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;

			xnew = Pbc(x[i] - vx[i] * dt*2);
			ynew = Pbc(y[i] - vy[i] * dt*2);
			znew = Pbc(z[i] - vz[i] * dt*2);
			
			xold[i]=x[i];
			yold[i]=y[i];
			zold[i]=z[i];
			
			x[i]=xnew;
			y[i]=ynew;
			z[i]=znew;
		}
	}

	else{
		cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
		double sumv[3] = {0.0, 0.0, 0.0};
		for (int i=0; i<npart; ++i){
			vx[i] = rnd.Rannyu() - 0.5;
			vy[i] = rnd.Rannyu() - 0.5;
			vz[i] = rnd.Rannyu() - 0.5;

			sumv[0] += vx[i];
			sumv[1] += vy[i];
			sumv[2] += vz[i];
		}
		for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
		double sumv2 = 0.0, fs;
		for (int i=0; i<npart; ++i){
			vx[i] = vx[i] - sumv[0];
			vy[i] = vy[i] - sumv[1];
			vz[i] = vz[i] - sumv[2];

			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		}
		sumv2 /= (double)npart;

		fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
		for (int i=0; i<npart; ++i){
			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;

			xold[i] = Pbc(x[i] - vx[i] * dt);
			yold[i] = Pbc(y[i] - vy[i] * dt);
			zold[i] = Pbc(z[i] - vz[i] * dt);
		}
	}
	return;
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
}

void Accumulate(void){ //Update block averages

	for(int i=0; i<n_props; ++i){
		blk_av[i] = blk_av[i] + walker[i];
	}
	blk_norm = blk_norm + 1.0;
}

void Move(void){ //Move particles with Verlet algorithm
	double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

	for(int i=0; i<npart; ++i){ //Force acting on particle i
		fx[i] = Force(i,0);
		fy[i] = Force(i,1);
		fz[i] = Force(i,2);
	}

	for(int i=0; i<npart; ++i){ //Verlet integration scheme

		xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(dt,2) );
		ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(dt,2) );
		znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(dt,2) );

		vx[i] = Pbc(xnew - xold[i])/(2.0 * dt);
		vy[i] = Pbc(ynew - yold[i])/(2.0 * dt);
		vz[i] = Pbc(znew - zold[i])/(2.0 * dt);

		xold[i] = x[i];
		yold[i] = y[i];
		zold[i] = z[i];

		x[i] = xnew;
		y[i] = ynew;
		z[i] = znew;
	}
	return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
	double f=0.0;
	double dvec[3], dr;

	for (int i=0; i<npart; ++i){
		if(i != ip){
			dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
			dvec[1] = Pbc( y[ip] - y[i] );
			dvec[2] = Pbc( z[ip] - z[i] );
	
			dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
			dr = sqrt(dr);

			if(dr < rcut){
				f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
			}
		}
	}
  
	return f;
}

void Measure(){ //Properties measurement
	int bin;
	double v, k, w, vij, wij;
	double dx, dy, dz, dr;
	

	v = 0.0; //reset observables
	k = 0.0;
	w = 0.0;
	for(int i=0;i<nbins;i++){
		walker[i+igofr]=0;
	}

//cycle over pairs of particles
	for (int i=0; i<npart-1; ++i){
		for (int j=i+1; j<npart; ++j){

			dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
			dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
			dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);
			
			for(bin=0;bin<nbins;bin++){
				if(dr>=(bin*bin_size) && dr<((bin+1)*bin_size)) walker[bin+igofr]+=2;
			}
//g(r)
			if(dr < rcut){
				vij = 1./pow(dr,12) - 1./pow(dr,6);
				wij = 1./pow(dr,12) - 0.5/pow(dr,6);
//Potential energy
				v+=vij;
				w+=wij;
			}
		}          
	}

//Kinetic energy
	for (int i=0; i<npart; ++i) k += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
	
	walker[iv] = 4.0*v; 					//Potential energy
	walker[ik] = k; 						//Kinetic energy
	walker[it] = (2.0 / 3.0) * k; 			//Temperature
	walker[ie] = walker[iv]+walker[ik];		//Total energy
	walker[iw] = 48.0 * w / 3.0; 			//Virial
						
	if(instant_print==1){
		ofstream IEpot, IEkin, IEtot, ITemp, IPres;

		IEpot.open("instant.epot.0",ios::app);
		IEkin.open("instant.ekin.0",ios::app);
		ITemp.open("instant.temp.0",ios::app);
		IEtot.open("instant.etot.0",ios::app);
		IPres.open("instant.pres.0",ios::app);
		
		IEpot<<walker[iv]/npart+vtail<< endl;
		IEkin<<walker[ik]/npart<<endl;
		ITemp<<walker[it]/npart<<endl;
		IEtot<<walker[ie]/npart<<endl;
		IPres<<rho*walker[it]+(walker[iw]+ptail*(double)npart)/vol<<endl;

		IEpot.close();
		IEkin.close();
		ITemp.close();
		IEtot.close();
	}
	
	if(print_eq==1){
		ofstream Tempeq;
		Tempeq.open("temp_eq.0", ios::app);
		Tempeq<<walker[it]/(double)npart<<endl;
		Tempeq.close();
	}
    return;
}

void Averages(int iblk){ //Print results for current block
	double r, gdir;
	ofstream Gofr, Gave, Epot, Ekin, Etot, Temp, Pres;
	const int wd=12;
    
    //cout << "Block number " << iblk << endl;
	//cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Epot.open("output.epot.0",ios::app);
	Ekin.open("output.ekin.0",ios::app);
	Etot.open("output.etot.0",ios::app);
	Temp.open("output.temp.0",ios::app);
	Pres.open("output.pres.0",ios::app);
    Gofr.open("output.gofr.0",ios::app);
    Gave.open("output.gave.0",ios::app);
    
    stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);

    stima_kin = blk_av[ik]/blk_norm/(double)npart; //Kinetic energy
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_kin=Error(glob_av[ik],glob_av2[ik],iblk);

    stima_etot = blk_av[ie]/blk_norm/(double)npart + vtail ;//Total energy
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);

    stima_temp = blk_av[it]/blk_norm/(double)npart; //Temperature
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);
	
	stima_pres = rho * stima_temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_pres=Error(glob_av[iw],glob_av2[iw],iblk);

    Epot<<iblk<<" "<<stima_pot<<" "<<glob_av[iv]/(double)iblk<<" "<<err_pot<<endl;
	Ekin<<iblk<<" "<<stima_kin<<" "<<glob_av[ik]/(double)iblk<<" "<<err_kin<<endl;
	Etot<<iblk<<" "<<stima_etot<<" "<<glob_av[ie]/(double)iblk<<" "<<err_etot<<endl;
	Temp<<iblk<<" "<<stima_temp<<" "<<glob_av[it]/(double)iblk<<" "<<err_temp<<endl;
    Pres<<iblk<<" "<<stima_pres<<" "<<glob_av[iw]/(double)iblk<<" "<<err_pres<<endl;
	
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
	
	if(iblk==nblock){
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

void ConfFinal(const char* filename){ //Write final configuration
	ofstream WriteConf;

	cout << "Print configuration to file " <<filename << endl;
	WriteConf.open(filename);

	for (int i=0; i<npart; ++i){
		WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
	}
	WriteConf.close();
	return;
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
	if (iblk==1) return 0;
	else return sqrt((sum2/(double)iblk - sum/(double)iblk*sum/(double)iblk)/(double)iblk);
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