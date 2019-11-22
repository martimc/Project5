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
	ofile << setprecision(8) << dt * tsteps << " ";
	ofile << setprecision(1) << n << endl;
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

void LU_decomp(double a, double b, double c, vector<double> &L, vector<double> &U){

	U[0]=b;

	for (int k = 1; k < L.size(); k++){

		L[k] = a/U[k-1];
		U[k] = b - L[k]*c;
		cout << L[k] << endl;
		cout << U[k] << endl;
	}
}

void solve_Ly_f(vector<double> &f, vector<double> &y, vector<double> &L){
	y[0] = f[0];

	for (int k=1; k< y.size(); k++){
		y[k] = f[k] - L[k]*y[k-1];
		cout << y[k] << endl;
	}

}

void solve_Uu_f(vector<double> &y, vector<double> &U, vector<double> &u, double c){
	u[u.size()] = y[u.size()]/U[u.size()];
	for (int k = y.size()-1; k != 0; k--){
		u[k] = y[k]/U[k] - c * u[k+1] / U[k];
		cout << u[k] << endl;
	}
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

	// read_input(tsteps, dx, dt);
	n = 1 / dx;
	alpha = dt / dx / dx;


	int length = 10;
	int a, b, c;
	a=1;
	b=2;
	c=1;
	vector<double> L(length,0);
	vector<double> U(length,0);
	LU_decomp(a,b,c, L, U);
	vector<double> y(length,0);
	vector<double> f(length,2);
	solve_Ly_f(f,y,L);
	vector<double> u(length,0);
	solve_Uu_f(y,U,u,c);
	// Forward_Euler(n, tsteps, dx, dt, alpha);
}
