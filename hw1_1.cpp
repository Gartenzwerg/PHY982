#include<iostream>
#include<cmath>
#include<fstream>

using namespace std;

double PI = 3.14159265;

long double f(int l, long double r, long double E);
void SecondOrderODESolver(long double* u, long double* u_derivative, double step, long N, int l, double E); 	//RK2
long double PhaseShift(long double R, double E, int l, double matching_radius); 	//return a phaseshift value in radians between PI/2 and PI/2

int main(){
	long N = 100001;
	long double h = 0.01;		//step size: 0.01fm
	long double *u_0, *u_1;		//u_0 is the original wave function; u_1 is the first derivative of the wave function.
	u_0 = new long double[N];
	u_1 = new long double[N];
	int l = 0;		//orbital angular momentum
	double E = 0.1;		//Energy in MeV

	//initialization
	u_0[0] = 0.0;
	u_1[0] = 1.0;

	//wave function for l = 0, E = 0.1MeV
	/*
	SecondOrderODESolver(u_0, u_1, h, N, l, E);
	for(int i = 0; i < N; i++) {
		cout<<i*h<<"	"<<u_0[i]<<endl;
	}
	cout<<PhaseShift(u_0[100000]/u_1[100000], E, l, 1000)<<endl;
	*/

	//output a data file of energies and corresponding phase shifts
	ofstream fout;
	fout.open("ps_0.dat");		
	for(int i = 0; i < 391; i++) {		//scan the energy interval between 0.1 MeV and 4.0 MeV
		SecondOrderODESolver(u_0, u_1, h, N, l, E);
		fout<<E<<"	"<<PhaseShift(u_0[100000]/u_1[100000], E, l, 1000)<<endl;
		E += 0.01;
	}
	fout.close();

}

long double f(int l, long double r, long double E){
	if(r < 0.001) {
		r = 0.0001;
	}
	return l*(l+1.0)/(r*r)+0.0478450*(-E-61.1/(1+exp((r-2.5853)/0.65)));
}

void SecondOrderODESolver(long double* u, long double* u_derivative, double step, long N, int l, double E){
	long double k0, k1;
	long double u0, u1; //intermediate values
	for(int i = 0; i < N; i++){
		k1 = step*f(l, i*step, E)*u[i];
		k0 = step*u_derivative[i];
		u1 = u_derivative[i] + k1/2;
		u0 = u[i] + k0/2;
		u_derivative[i+1] = u_derivative[i]+step*f(l, (i+0.5)*step, E)*u0;
		u[i+1] = u[i] + step*u1;
	}
}

//Note that the returned value is between -PI/2 and PI/2; to remove discontinuities, one may need to add or subtract PI by hand.
long double PhaseShift(long double R, double E, int l, double matching_radius){
	long double wave_number = sqrt(0.0478450*E);
	long double shift = atan(wave_number*R)+l*PI/2-wave_number*matching_radius;
	long double res = fmod(shift, PI);
	if(res < -PI/2) {
		return res + PI;
	}else if(res > PI/2) {
		return res - PI;
	}else{
		return res;
	}
}