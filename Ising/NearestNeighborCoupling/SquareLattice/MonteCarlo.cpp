// Monte Carlo code by Addison Richards


#include <fstream>
#include <math.h>
#include <random>

// Lattice size
const int L = 64;
const int N = L*L;

// Temperature grid
const double TN = 64;
const double Tmax = 3.0;
const double Tmin = 0.1;

// Mersenne twister RNG
std::random_device rd;
const int seed = rd();
std::mt19937 gen(seed);
std::uniform_real_distribution<double> dis(0.0, 1.0);

// Monte Carlo sweeps
const int eqmcs = 10000; // equilibration
const int mcs = 20000; // measurement
const int runs = 10; // number of independent runs

// output
const int binSize = 10000; // output sample average size
std::string DATADIR = "./data/";
std::string OUTFILE = DATADIR + "L=" + std::to_string(L) + ".txt";

void out(double &m,double &m2,double &u,double &u2,int &samples, double T)
{
	m2 /= N; u2 /= N;
	m /= (N*binSize);  m2 /= (N*binSize);
	u /= (N*binSize);  u2 /= (N*binSize);

	std::ofstream binOut;
	binOut.open(OUTFILE, std::ios_base::app);
	binOut << T << ',' << m << ',' << m2 << ',' << u << ',' << u2 << '\n';
	binOut.close();

	m = 0;  m2 = 0;
	u = 0;  u2 = 0;
	samples = 0;
}

void MCSweep(double &E,double &M, double S[L][L], double T)
{
	for(int n=0; n<N; n++)
	{
		int rx = int(L*dis(gen));
		int ry = int(L*dis(gen));
		
		double delE = 2*S[rx][ry]*( S[(rx+1)%L][ry] + S[rx][(ry+1)%L] + S[(rx+L-1)%L][ry] + S[rx][(ry+L-1)%L]);

		if(delE < 0)
		{
			S[rx][ry] = -S[rx][ry]; 
			M += 2.0*S[rx][ry]; 
			E += delE;
		}
		else if(dis(gen) < exp(-delE/T))
		{
			S[rx][ry] = -S[rx][ry]; 
			M += 2.0*S[rx][ry]; 
			E += delE;
		}
	}
}

void Metropolis(double T)
{
	//------ Initializing lattice ------//
	double S[L][L];
	for(int x=0; x<L; x++){for(int y=0; y<L; y++){S[x][y] = 1;}}
	double E = -2*N;
	double M = N;

	//------ Equilibration ------//
	for(int i=0; i<eqmcs; i++){MCSweep(E,M,S,T);}

	//------ Measurement ------//
	double m=0, m2=0, u=0, u2=0;
	int samples = 0;

	for(int i=0; i<mcs; i++)
	{
		MCSweep(E,M,S,T);

		m += fabs(M);  	m2 += M*M;
		u += E;  		u2 += E*E;
		samples++;
		
		if(samples % binSize == 0){out(m, m2, u, u2, samples, T);}
	}
}

int main()
{
	for(int run=0; run<runs; run++){
		for(int Ti=0; Ti<TN; Ti++)
		{
			double T = Tmax - Ti*(Tmax-Tmin)/(TN-1);
			Metropolis(T);
		}	
	}
}