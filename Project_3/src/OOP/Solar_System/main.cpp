#include <iostream>
#include <celestial.h>

using namespace std;

int main()
{
    double M = 56.0;
    double x_0=66.0;
    double y_0=67.0;
    double Vx_0=68.0;
    double Vy_0=69.0;
    //Create object
    celestial Earth(M, x_0, y_0, Vx_0, Vy_0);

    double t_max = 3.0;
    int N = 100;
    //char outfile_verlet = 'verlet';
    char outfile_euler = 'e';
    //Calls Verlet and Euler methods
    //Earth.Verlet(N, t_max, outfile_verlet);
    Earth.Euler(N, t_max, outfile_euler);
    Earth.print();
    return 0;
}

