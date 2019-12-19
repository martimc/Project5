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

void tridiag(double a, double b, double c, double* y, double* &u, int n) {
	double denom;
	double *c_new;
	c_new = new double[n];
	if (b == 0) throw("error 1 in tridiag");
	u[0] = y[0]; u[n] = y[n];

	denom = b;
	u[1] = y[1] / denom;
	c_new[1] = c / denom;
	for (int i = 2; i < n; i++) {
		denom = b - a * c_new[i - 1];
		if (denom == 0) throw("error 2 in tridiag");
		u[i] = (y[i] - a * u[i - 1]) / denom;
		c_new[i] = c / denom;
	}
	for (int j = (n - 1); j > 0; j--) {
		u[j] -= c_new[j] * u[j + 1];
	}
	//Output(n, u, 1, 0.005);
}

void Forward_Euler(int n, int tsteps, double dt, double alpha) {
	double *u, *y;
	u = new double[n+1]; y = new double[n + 1];

	u[n] = 1; y[n] = 0; //boundary conditions and initial conidtions
	for (int i = 0; i < n; i++) {
		u[i] = 0;
		y[i] = 0;
	}
	for (int t = 1; t <= tsteps; t++) { //Forward Euler algorithm, y = Au
		for (int i = 1; i < n; i++) {
			y[i] = alpha * u[i - 1] + (1 - 2 * alpha) * u[i] + alpha * u[i + 1];
		}
		//Output(n, u);
		for (int i = 1; i < n; i++) {
			u[i] = y[i];
		}
	}
	Output(n, u, tsteps, dt);
}

void Backward_Euler(int n, int tsteps, double dt, double alpha){
	double a, b, c;
	double *u, *y;
	u = new double[n + 1]; y = new double[n + 1];
	a = -alpha; c = -alpha;
	b = 1 + 2 * alpha;

	for(int i=0; i<=n; i++){
		u[i]=0;
		y[i]=0;
	}
	y[n] = 1;
	for(int t=1; t<=tsteps; t++){
		tridiag(a, b, c, y, u, n);
		for (int j = 0; j <= n; j++) {
			y[j] = u[j];
		}
	}
	Output(n, u, tsteps, dt);
}

void Crank_Nicolson(int n, int tsteps, double dt, double alpha){
	double a, b, c;
	double *u, *y;
	u = new double[n + 1]; y = new double[n + 1];
	a = -alpha; c = -alpha;
	b = 2 + 2 * alpha;
	for (int i = 0; i < n; i++) {
		u[i] = 0;
	}
	u[n] = 1;

	for (int t = 1; t <= tsteps; t++) {
		for (int i = 1; i < n; i++) { // using forward euler to find the vector y in the equation Au = y
			y[i] = alpha * u[i - 1] + (2 - 2 * alpha) * u[i] + alpha * u[i + 1];
		}
		y[0] = 0; y[n] = 1;
		// Tridiagonal solver for finding u in the equaiton Au = y
		tridiag(a, b, c, y, u, n);
	}
	Output(n, u, tsteps, dt);
}

void Output(int n, double* u, int tsteps, double dt) {
	//function for writing the results in a file
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	for (int i = 0; i < n; i++) {
		ofile << setprecision(8) << u[i] << " ";
	}
	ofile << setprecision(8) << u[n] << endl;
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

	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << setprecision(8) << dt * tsteps << " ";
	ofile << setprecision(1) << n << endl;

	Forward_Euler(n, tsteps, dt, alpha);
	Backward_Euler(n, tsteps, dt, alpha);
	Crank_Nicolson(n, tsteps, dt, alpha);
}
