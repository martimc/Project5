#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include "lib.cpp"
#include <sstream>
#include <string>

void Forward_Euler(int n, int tsteps, double dt, double alpha) {
	double **u=0;
	double **y=0;
	int dx;
	int dy;
	dx = 1/n;
	dy=dx;
	u = new double*[n+1];
	y = new double*[n+1];
	double u_xx, u_yy;
	//initialize matrix
	for (int i = 0; i<=n+1; i++){
		u[i] = new double[n+1];
		y[i] = new double[n+1];
	}
	//boundary conditions and initial conditions
	for (int i = 0; i<=n; i++){
		for (int j = 0; j<=n; j++){
			u[i][j]=0;
			y[i][j]=0;
		}
	}
	for (int i = 0; i<=n; i++){
		u[i,n]=(double)1;
	}

	for (int t = 1; t <= tsteps; t++) { //Forward Euler algorithm, y = Au
		for (int i = 1; i < n; i++) {
			for (int j = 1; j <n; j++){
				u_xx = (u[i+1][j] - 2*u[i][j] + u[i-1][j])/dx/dx;
				u_yy = (u[i][j+1] - 2*u[i][j] + u[j][j-1])/dy/dy;
				y[i][j] = u[i][j] + dt*(u_xx + u_yy);
			}
		}
		for (int i = 1; i < n; i++) {
			for (int j = 1; i < n; i++){
				u[i,j]=y[i,j];
			}
		}
	}
}

void read_input(int& tsteps, double& dx, double& dt) {
	cout << "number of time steps: ";
	cin >> tsteps;
	cout << "dx and dy: ";
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

	read_input(tsteps, dx, dt);
	Forward_Euler(n, tsteps, dt, alpha);

}
