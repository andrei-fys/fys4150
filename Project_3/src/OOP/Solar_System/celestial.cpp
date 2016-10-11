#include "celestial.h"
#include <iostream>
#include <time.h>
#include <math.h>
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
    double *grid_points = new double[N];
    double t_start = 0.0;
    double h=((double) (t_finish - t_start)/N);

    for (int i=0; i<N; i++) {
        grid_points[i] =  t_start + i*h;
    }

    double *x = new double[N];
    double *y = new double[N];
    double *Vx = new double[N];
    double *Vy = new double[N];
    double a = 4.0*M_PI*M_PI*h;
    clock_t g_start, g_finish;
    x[0]=startx;
    y[0]=starty;
    Vx[0]=startvx0;
    Vy[0]=startvy0;

    g_start = clock();
    // EULER START
    for (int i=0;i<N;i++){
        x[i+1]=x[i] + Vx[i]*h;
        Vx[i+1]=Vx[i] - a*x[i];
        y[i+1]=y[i] + Vy[i]*h;
        Vy[i+1]=Vy[i] - a*y[i];
    }
    // EULER END
    g_finish=clock();

    for (int i=0; i<N; i++) {
            cout << grid_points[i] << "," << x[i] << "," << y[i] << endl;
        }

    delete [] grid_points;
    delete [] x;
    delete [] y;
    delete [] Vx;
    delete [] Vy;

    cout << "Time of Euler " << ((double) (t_finish-t_start)/CLOCKS_PER_SEC) << endl;

}

int celestial::Verlet(int N, double t_finish, char output_filename )
{
    double *grid_points = new double[N];
    double t_start = 0.0;
    double h=((double) (t_finish - t_start)/N);

    for (int i=0; i<N; i++) {
        grid_points[i] = t_start + i*h;
    }

    double *x = new double[N];
    double *y = new double[N];
    double *Vx = new double[N];
    double *Vy = new double[N];
    double a = 2.0*M_PI*M_PI*h;
    double b = 1.0 - a*h;
    clock_t g_start, g_finish;
    x[0]=startx;
    y[0]=starty;
    Vx[0]=startvx0;
    Vy[0]=startvy0;

    g_start = clock();
    // VERLET START
    for (int i=0;i<N;i++){
        x[i+1]=x[i]*b + Vx[i]*h;
        Vx[i+1]=Vx[i] - a*(x[i+1]+x[i]);
        y[i+1]=y[i]*b + Vy[i]*h;
        Vy[i+1]=Vy[i] - a*(y[i+1]+y[i]);
    }
    // VERLET END
    g_finish=clock();

    for (int i=0; i<N; i++) {
            cout << grid_points[i] << "," << x[i] << "," << y[i] << endl;
        }

    delete [] grid_points;
    delete [] x;
    delete [] y;
    delete [] Vx;
    delete [] Vy;

    cout << "Time of Verlet(velocity) " << ((double) (t_finish-t_start)/CLOCKS_PER_SEC) << endl;
}


int celestial::print()
{
    cout << " X_0 === " << startx << endl;
}
