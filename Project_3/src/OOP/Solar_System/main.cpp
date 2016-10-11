#include <iostream>
#include <celestial.h>
#include <math.h>
using namespace std;

int main()
{
    double M = 1.0;
    double x_0=1.0;
    double y_0=0.0;
    double Vx_0=0.0;
    double Vy_0=2.0*M_PI;
    //Create object
    celestial Earth(M, x_0, y_0, Vx_0, Vy_0);

    double t_max = 1.0;
    int N = 100;
    char outfile_verlet = 'verlet';
    char outfile_euler = 'euler';
    //Calls Verlet and Euler methods
    Earth.Verlet(N, t_max, outfile_verlet);
    Earth.Euler(N, t_max, outfile_euler);
    //Earth.print();
    return 0;
}

