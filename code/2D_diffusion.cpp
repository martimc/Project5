#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include "lib.cpp"
#include <sstream>
#include <string>

using namespace std;
ofstream ofile;
void Output(int, double**, int, double);

void Forward_Euler(int n, int tsteps, double dt, double alpha) {
	double **u;
	double **y;
	double u_xx, u_yy;
	u = new double* [n + 1]; y = new double* [n + 1];
	for (int i = 0; i <= n; i++) {
		u[i] = new double[n + 1];
		y[i] = new double[n + 1];
	}
	
	//boundary conditions and initial conditions
	for (int i = 0; i<=n; i++){
		for (int j = 0; j<=n; j++){
			u[i][j]=0;
			y[i][j]=0;
		}
	}
	
	for (int i = 0; i<=n; i++){
		u[i][n]=1;
	}

	for (int t = 1; t <= tsteps; t++) { //Forward Euler algorithm, y = Au
		for (int i = 1; i < n; i++) {
			for (int j = 1; j < n; j++) {
				u_xx = (u[i+1][j] - 2*u[i][j] + u[i-1][j]);
				u_yy = (u[i][j+1] - 2*u[i][j] + u[j][j-1]);
				y[i][j] = u[i][j] + alpha*(u_xx + u_yy);
				//y[i][j] = alpha * u[i + 1][j] + alpha * u[i][j + 1] + (1 - 4 * alpha) * u[i][j] + alpha * u[i - 1][j] + alpha * u[i][j - 1];
			}
		}
		for (int i = 1; i < n; i++) {
			for (int j = 1; j < n; j++){
				u[i][j]=y[i][j];
			}
		}
	}
	Output(n, u, tsteps, dt);
}

void read_input(int& tsteps, double& dx, double& dt) {
	cout << "number of time steps: ";
	cin >> tsteps;
	cout << "dx and dy: ";
	cin >> dx;
	cout << "dt: ";
	cin >> dt;
}

void Output(int n, double** u, int tsteps, double dt) {
	//function for writing the results in a file
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << setprecision(8) << dt * tsteps << " ";
	ofile << setprecision(8) << n << endl;
	printf("start");
	for (int i = 0; i <= n; i++) {
		for (int j = 0; j < n; j++) {
			ofile << setprecision(8) << u[i][j] << " ";
		}
		ofile << setprecision(8) << u[i][n] << endl;
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

	read_input(tsteps, dx, dt);
	n = 1 / dx;
	alpha = dt / dx / dx;

	ofile << setiosflags(ios::showpoint | ios::uppercase);
	//ofile << setprecision(8) << alpha << endl;

	Forward_Euler(n, tsteps, dt, alpha);
}
