#include "celestial.h"
#include <iostream>
using namespace std;

celestial::celestial(double M, double x0, double y0, double Vx0, double Vy0)
{
    mass=M;
    startx=x0;
    starty=y0;
    startvx0=Vx0;
    startvy0=Vy0;
}

int celestial::Euler(int N, double t_finish, char output_filename )
{
    cout << "N = " << N << ", t_finish = " << t_finish << ", output = " << output_filename << endl;
}

int celestial::print()
{
    cout << " X_0 === " << startx << endl;
}
