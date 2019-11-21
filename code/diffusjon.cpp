#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include "lib.cpp"
#include <sstream>
#include <string>
#include <random>
//#include <armadillo>

using namespace std;
ofstream ofile;
void Output(int n, double* u, int tsteps, double dt);

void Forward_Euler(int n, int tsteps, double dx, double dt, double alpha) {
	double *u, *u_new;
	u = new double[n+1]; u_new = new double[n + 1];

	u[n] = 1; u_new[n] = 0; //boundary conditions and initial conidtions
	for (int i = 0; i < n; i++) {
		u[i] = 0; u_new[i] = 0;
	}
	for (int t = 1; t <= tsteps; t++) { //Forward Euler algorithm
		for (int i = 1; i < n; i++) {
			u_new[i] = alpha * u[i - 1] + (1 - 2 * alpha) * u[i] + alpha * u[i + 1];
		}
		//Output(n, u);
		for (int i = 1; i < n; i++) {
			u[i] = u_new[i];
		}
	}
	Output(n, u, tsteps, dt);
}

void Output(int n, double* u, int tsteps, double dt) {
	//function for writing the results in a file
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	//ofile << setprecision(8) << "t = " << dt * tsteps << endl;
	for (int i = 0; i < n; i++) {
		ofile << setw(4) << setprecision(8) << " " << u[i];
	}
	ofile << setw(4) << setprecision(8) << " " << u[n] << endl;
}

void read_input(int& tsteps, double& dx, double& dt) {
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

	read_input(tsteps, dx, dt);
	n = 1 / dx;
	alpha = dt / dx / dx;
	
	Forward_Euler(n, tsteps, dx, dt, alpha);
}