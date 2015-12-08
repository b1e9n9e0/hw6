/*
4th order runge kutta
b1e9n9e0
03.12.15


iterate from x= 0 to x = 100
dx from 0.1 to 0.01
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>


using namespace std;

// global variabels
static double k_1[3];
static double k_2[3];
static double k_3[3];
static double k_4[3];

static double a = 10;
static double b = 28;
static double c = 8./3;

void func(const double* const y, double* const f) {
	f[0] = a * (y[1]-y[0]);
	f[1] = y[0] * (b-y[2]) - y[1];
	f[2] = y[0]*y[1] - c*y[2];
}

void fill_k (const double* const y, const double dt) {

	double temp[3];

	// calculate the first k
	func(y, k_1);

	// fill array with new values and calculate the next k
	temp[0] = y[0] + 0.5*dt*k_1[0];
	temp[1] = y[1] + 0.5*dt*k_1[1];
	temp[2] = y[2] + 0.5*dt*k_1[2];

	func(temp, k_2);

	// fill array with new values and calculate the next k
	temp[0] = y[0] + 0.5*dt*k_2[0];
	temp[1] = y[1] + 0.5*dt*k_2[1];
	temp[2] = y[2] + 0.5*dt*k_2[2];

	func(temp, k_3);

	// fill array with new values and calculate the next k
	temp[0] = y[0] + dt*k_3[0];
	temp[1] = y[1] + dt*k_3[1];
	temp[2] = y[2] + dt*k_3[2];

	func(temp, k_4);
}

int main(int argcount, char** argvector){


	// defaults
	const double t_min = 0;
	const double t_max = 100;
	double y[3] = {1.,1.,1.};

	double dt;

	// variables to build the filename
	const string prefix = "stepsize_";
	const string suffix = ".dat";

	// streams
	stringstream stream;
	ofstream out;

	// check if the right amount of arguments is given
	if (argcount != 2) {
		// print error message, return to caller
		cout << "usage: " << argvector[0] << " stepsize" << endl;
		return -1;
	}

	// get the value of dt and implement into double value
	dt = atof(argvector[1]);

	// build the filename and open the file
	stream << prefix << dt << suffix;
	out.open(stream.str().c_str());

	// output first line to file
	out << t_min << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << endl;

	// iterate over all timesteps
	for (double t=t_min; t<=t_max; t+=dt) {

		// get the k values
		fill_k(y, dt);

		// calculate the next step of y
		y[0] += dt/6. * (k_1[0] + 2*k_2[0] + 2*k_3[0] + k_4[0]);
		y[1] += dt/6. * (k_1[1] + 2*k_2[1] + 2*k_3[1] + k_4[1]);
		y[2] += dt/6. * (k_1[2] + 2*k_2[2] + 2*k_3[2] + k_4[2]);

		// output to file
		out << t << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << endl;
	}

	// close output stream
	out.close();

return 0;
}