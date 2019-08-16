#include <fstream>
#include <math.h>
#include <random>

// Lattice size
const int L = 64;
const int N = L*L;

// Coupling
const double Jx = 1.0;
const double Jy = 1.0;

// Temperature grid
const double TN = 64;
const double Tmax = 0.5;
const double Tmin = 0.01;

// Mersenne twister RNG
std::random_device rd;
const int seed = rd();
std::mt19937 gen(seed);
std::uniform_real_distribution<double> dis(0.0, 1.0);

// Monte Carlo sweeps
const int eqmcs = 100000; // equilibration
const int mcs = 100000; // measurement
const int runs = 100; // number of independent runs

// output
const int binSize = 100000; // output sample average size
std::string DATADIR = "./data/";
std::string OUTFILE = DATADIR + "L=" + std::to_string(L) + ".txt";

void out(double &q,double &q2,double &u,double &u2,int &samples, double T)
{
	q2 /= N; u2 /= N;
	q /= (N*binSize);  q2 /= (N*binSize);
	u /= (N*binSize);  u2 /= (N*binSize);

	std::ofstream binOut;
	binOut.open(OUTFILE, std::ios_base::app);
	binOut << T << ',' << q << ',' << q2 << ',' << u << ',' << u2 << ',' << 0 << '\n';
	binOut.close();

	q = 0;  q2 = 0;
	u = 0;  u2 = 0;
	samples = 0;
}

void MCSweep(double &E,double &Q, double S[L][L][2], double T)
{
	for(int n=0; n<N; n++)
	{
		int rx = int(L*dis(gen));
		int ry = int(L*dis(gen));

		double Sxo = S[rx][ry][0];
		double Syo = S[rx][ry][1];
		
		double th = 2*M_PI*dis(gen);
		double Sxn = cos(th);
		double Syn = sin(th);

		double delE = -Jx*(Sxn-Sxo)*(S[(rx+1)%L][ry][0] + S[(rx+L-1)%L][ry][0])
					  -Jy*(Syn-Syo)*(S[rx][(ry+1)%L][1] + S[rx][(ry+L-1)%L][1]);

		if(delE < 0)
		{
			S[rx][ry][0] = Sxn; 
			S[rx][ry][1] = Syn;
			Q += (Sxn*Sxn-Syn*Syn) - (Sxo*Sxo-Syo*Syo);
			E += delE;
		}
		else if(dis(gen) < exp(-delE/T))
		{
			S[rx][ry][0] = Sxn; 
			S[rx][ry][1] = Syn;
			Q += (Sxn*Sxn-Syn*Syn) - (Sxo*Sxo-Syo*Syo); 
			E += delE;
		}
	}

}

void Metropolis(double T)
{
	//------ Initializing lattice ------//
	double S[L][L][2];
	for(int x=0; x<L; x++){for(int y=0; y<L; y++){S[x][y][0] = 1; S[x][y][1] = 0;}}
	double E = -N;
	double Q = N;

	//------ Equilibration ------//
	for(int i=0; i<eqmcs; i++){MCSweep(E,Q,S,T);}

	//------ Measurement ------//
	double q=0, q2=0, u=0, u2=0;
	int samples = 0;

	for(int i=0; i<mcs; i++)
	{
		MCSweep(E,Q,S,T);

		q += fabs(Q);  	q2 += Q*Q;
		u += E;  		u2 += E*E;

		samples++;
		
		if(samples % binSize == 0){out(q, q2, u, u2, samples, T);}
	}
}

int main()
{
	for(int run=0; run<runs; run++){
		for(int Ti=0; Ti<TN; Ti++)
		{
			double T = Tmin + Ti*(Tmax-Tmin)/(TN-1);
			Metropolis(T);
		}	
	}
}
