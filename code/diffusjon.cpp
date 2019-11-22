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

struct Variables{
	double *L;
	double *U;
	double *y;
	double *f;
	double *u;
	int length;
	int a;
	int b;
	int c;
};

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

void LU_decomp(Variables& var){

	var.U[0]=var.b;

	for (int k = 1; k < var.length; k++){

		var.L[k] = var.a/var.U[k-1];
		var.U[k] = var.b - var.L[k]*var.c;
	}
}

void solve_Ly_f(Variables& var){
	var.y[0] = var.f[0];

	for (int k=1; k< var.length; k++){
		var.y[k] = var.f[k] - var.L[k]*var.y[k-1];
	}

}

void solve_Uu_f(Variables& var){
	var.u[var.length-1] = var.y[var.length-1]/var.U[var.length-1];
	for (int k = var.length-2; k != -1; k--){
		var.u[k] = var.y[k]/var.U[k] - var.c * var.u[k+1] / var.U[k];
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



	Variables sol;
	sol.length = 3;
	sol.L = new double[sol.length];
	sol.U = new double[sol.length];
	sol.y = new double[sol.length];
	sol.f = new double[sol.length];
	sol.u = new double[sol.length];
	sol.a=1;
	sol.b=2;
	sol.c=1;
	for (int i=0; i<sol.length; i++){
		sol.f[i] = 2;
	}

	LU_decomp(sol);
	solve_Ly_f(sol);
	solve_Uu_f(sol);


	// Forward_Euler(n, tsteps, dx, dt, alpha);
}
