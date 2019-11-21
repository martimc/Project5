#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include "lib.cpp"
#include <sstream>
#include <string>
#include <random>
#include <armadillo>

using namespace std;
ofstream ofile;

void Output(int n, double* u);
//random_device rd;
//mt19937_64 gen(rd());
//uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
//uniform_int_distribution<int> RandomPosition(0, size-1);
//looping over all spins

/*void Forward_Euler(int n, int tsteps, double dx, double dt, double alpha) {
	int size = n;
	double u[size + 1], u_new[size + 1];
	u[n] = 1; u_new[n] = 0;
	for (int i = 0; i < n; i++) {
		u[i] = 0; u_new[i] = 0;
	}
	for (int t = 1; t <= tsteps; t++) {
		for (int i = 1; i < n; i++) {
			u_new[i] = alpha * u[i - 1] + (1 - 2 * alpha) * u[i] + alpha * u[i + 1];
		}
		for (int i = 1; i < n; i++) {
			u[i] = u_new[i];
		}
	}
	Output(n, u);
}*/

void Output(int n, double* u) {
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	for (int i = 0; i < n; i++) {
		ofile << setw(8) << setprecision(8) << " " << u[i];
	}
	ofile << setw(8) << setprecision(8) << " " << u[n] << endl;
}

void read_input(int& n, int& tsteps, double& dx, double& dt) {
	cout << "size of array: ";
	cin >> n;
	cout << "number of time steps: ";
	cin >> tsteps;
	cout << "dx: ";
	cin >> dx;
	cout << "dt: ";
	cin >> dt;
}

int main(int argc, char* argv[]) {
	char* outfilename;
	int n, tsteps;
	double dx, dt, alpha;
	
	if (argc <= 1) {
		cout << "Bad Usage: " << argv[0] << " read also output file on same line" << endl;
		exit(1);
	}
	else {
		outfilename = argv[1];
	}

	ofile.open(outfilename);

	read_input(n, tsteps, dx, dt);
	alpha = dt / dx / dx;

	vec u(n + 1);

	for (int i = 0; i < n; i++) {
		u[i] = 0;
	}
	u[n] = 1;

	Output(n, u);
	//Forward_Euler(n, tsteps, dx, dt, alpha);
}