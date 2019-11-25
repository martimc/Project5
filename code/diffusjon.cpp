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
	double *u;
	double *u_new;
	int length;
	int tsteps;
	double dx;
	double dt;
	double alpha;
	int a;
	int b;
	int c;
};

void LU_decomp(Variables& var){

	var.U[0]=var.b;

	for (int k = 1; k < var.length; k++){

		var.L[k] = var.a/var.U[k-1];
		var.U[k] = var.b - var.L[k]*var.c;
		// cout << var.L[k] << endl;
		// cout << var.U[k] << endl;
	}
}

void solve_Ly_u(Variables& var){
	var.y[0] = var.u[0];
	// cout << var.u_new[0] << endl;
	// cout << var.y[0] << endl;
	for (int k=1; k< var.length; k++){
		var.y[k] = var.u[k] - var.L[k]*var.y[k-1];
		// cout << var.y[k] << endl;
	}

}

void solve_UuNew_u(Variables& var){
	var.u_new[var.length-1] = var.y[var.length-1]/var.U[var.length-1];
	// cout << var.u_new[var.length-1];
	for (int k = var.length-2; k != -1; k--){
		var.u_new[k] = var.y[k]/var.U[k] - var.c * var.u_new[k+1] / var.U[k];
		// cout << var.u_new[k] << endl;
	}
}

void thomas_algorithm(Variables& var){
	LU_decomp(var);
	solve_Ly_u(var);
	solve_UuNew_u(var);
}

void Forward_Euler(int n, int tsteps, double dx, double dt, double alpha) {
	double *u;
	double *u_new;
	u = new double[n+1]; u_new = new double[n + 1];

	u[n] = 1; u_new[n] = 0; //boundary conditions and initial conidtions
	for (int i = 0; i < n; i++) {
		u[i] = 0;
		u_new[i] = 0;
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
	var.u[var.length]=1; // boundary conditions and initial conditions

	for(int i=0; i<var.length; i++){
		var.u[i]=0;
		var.u_new[i]=0;
	}
	for(int t=1; t<var.tsteps; t++){
		LU_decomp(var);
		solve_Ly_u(var);
		solve_UuNew_u(var);
		for(int i=1; i<var.length; i++){
			var.u[i]=var.u_new[i];
		}
	}
	Output(var.length, var.u, var.tsteps, var.dt);
}

// void Crank_Nicholsen(Variables& var){
// 	var.u[var.length-1]=1;
// 	var.u_new[var.length-1]=0
// }

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

	read_input(tsteps, dx, dt);

	Variables sol;
	sol.dt = dt;
	sol.dx = dx;
	sol.tsteps = tsteps;
	sol.alpha = sol.dt/sol.dx/sol.dx;
	sol.length = (int)(1/sol.dx);
	sol.L = new double[sol.length+1];
	sol.U = new double[sol.length+1];
	sol.y = new double[sol.length+1];
	sol.u = new double[sol.length+1];
	sol.u_new = new double[sol.length+1];

	sol.a = -sol.alpha;
	sol.c = -sol.alpha;
	sol.b = 1 + 2*sol.alpha;






	if (argc <= 1) {
		cout << "Bad Usage: " << argv[0] << " read also output file on same line" << endl;
		exit(1);
	}
	else {
		outfilename = argv[1];
	}

	ofile.open(outfilename);

	Backward_Euler(sol);
	Forward_Euler(n, tsteps, dx, dt, alpha);
}
