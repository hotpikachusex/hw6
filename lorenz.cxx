#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace std;



int main(){

 
 
  int T = 10;
  double a = 8.0, b = 28, c = 8/3;
  double dt = 1/1000.0;
  int N = T/dt;
  

 double* x = new double[N];
 double* y = new double[N];
 double* z = new double[N];
 x[0] = 1;
 y[0] = 1;
 z[0] = 1;
 double k1x = 0.0, k2x = 0.0,k3x = 0.0,k4x = 0.0;
 double k1y = 0.0, k2y = 0.0,k3y = 0.0,k4y = 0.0;
 double k1z = 0.0, k2z = 0.0,k3z = 0.0,k4z = 0.0;

 ofstream out("lorenz.txt");

 for(int i = 0; i<N-1;i++){

		  
	k1x = a*(y[i]-x[i]);
	k2x = a*(y[i]-(x[i]+dt/2 * k1x));
	k3x = a*(y[i]-(x[i]+dt/2 * k2x));
	k4x = a*(y[i]-(x[i]+dt*k3x));
	x[i+1] = x[i] + dt/6 * (k1x+2*k2x+2*k3x+k4x);

	k1y = x[i] * (b-z[i]) - y[i];
	k2y = x[i] * (b-z[i]) - (y[i]+dt/2 * k1y);
	k3y = x[i] * (b-z[i]) - (y[i]+dt/2 * k2y);
	k4y = x[i] * (b-z[i]) - (y[i]+dt * k3y);
	y[i+1] = y[i] + dt/6 * (k1y+2*k2y+2*k3y+k4y);

	k1z = x[i] * y[i] - c * z[i];
	k2z = x[i] * y[i] - c * (z[i]+dt/2 * k1z);
	k3z = x[i] * y[i] - c * (z[i]+dt/2 * k2z);
	k4z = x[i] * y[i] - c * (z[i]+dt * k3z);
	z[i+1] = z[i] + dt/6 * (k1z+2*k2z+2*k3z+k4z);

	
	
	out << i*dt << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << endl;

 }


 out.close();

 delete[] x;
 delete[] y;
 delete[] z;

 return 0;
}