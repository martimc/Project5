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
	double *u_new;
	int length;
	int a;
	int b;
	int c;
	int n;
	int tsteps;
	double dx;
	double dt;
	double alpha;
};

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

void Forward_Euler(int n, int tsteps, double dt, double alpha) {
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

void Backward_Euler(Variables& var){
	var.u[var.length-1]=1; // boundary conditions and initial conditions
	var.u_new[var.length-1]=0;
	for(int i=1; i<var.length; i++){
		var.u[i]=0;
		var.u_new[i]=0;
	}
	for(int t=1; t<var.tsteps; t++){
		LU_decomp(var);
		solve_Ly_f(var);
		solve_Uu_f(var);
		var.u[0]=0;
		var.u[var.length-1]=1; //boundary conditions again
		for(int i=0; i<var.length; i++){
			var.y[i]=var.u[i];
		}
	}
	// Output(n, u, tsteps, dt);
}

void tridiag(double a, double b, double c, double* y, double* &u, int n) {
	double denom;
	double *c_new;
	c_new = new double[n];
	if (b == 0) throw("error 1 in tridiag");
	denom = b;

	u[0] = y[0] / denom;
	for (int i = 1; i <= n; i++) {
		c_new[i-1] = c / denom;
		denom = b - a * c_new[i-1];
		if (denom == 0) throw("error 2 in tridiag");
		u[i] = (y[i] - a * u[i - 1]) / denom;
	}
	for (int j = (n - 1); j >= 0; j--) {
		u[j] -= c_new[j] * u[j + 1];
	}
	//Output(n, u, 1, 0.005);
}

void Crank_Nicholsen(int n, int tsteps, double dt, double alpha){
	double a, b, c;
	double *u, *y;
	u = new double[n + 1]; y = new double[n + 1];
	a = -alpha; c = -alpha;
	b = 2 + 2 * alpha;
	for (int i = 0; i <= n; i++) {
		u[i] = 0;
	}

	for (int t = 1; t <= tsteps; t++) {
		for (int i = 1; i < n; i++) { // using forward euler to find the vector y in the equation Au = y
			y[i] = alpha * u[i - 1] + (2 - 2 * alpha) * u[i] + alpha * u[i + 1];
		}
		y[0] = 0; y[n] = 1;
		// Tridiagonal solver for finding u in the equaiton Au = y
		tridiag(a, b, c, y, u, n);
		u[0] = 0; u[n] = 1;
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



int main(int argc, char* argv[]) {
	char* outfilename;
	int n, tsteps;
	double dx, dt, alpha;

	/*Variables sol;
	sol.length = 3;
	sol.L = new double[sol.length];
	sol.U = new double[sol.length];
	sol.y = new double[sol.length];
	sol.f = new double[sol.length];
	sol.u = new double[sol.length];
	sol.u_new = new double[sol.length];
	sol.y = new double[sol.length];
	sol.a=1;
	sol.b=2;
	sol.c=1;
	for (int i=0; i<sol.length; i++){
		sol.f[i] = 2;
		}*/

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




	//LU_decomp(sol);
	//solve_Ly_f(sol);
	//solve_Uu_f(sol);

	//Backward_Euler(sol);
	
	Forward_Euler(n, tsteps, dt, alpha);
	Crank_Nicholsen(n, tsteps, dt, alpha);
}
