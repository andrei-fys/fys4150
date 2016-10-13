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

    //cout << "Time of Verlet(velocity) " << ((double) (t_finish-t_start)/CLOCKS_PER_SEC) << endl;
}

int celestial::VerletTwoBody(int N, double t_finish, char output_filename )
{
    double *grid_points = new double[N];
    double t_start = 0.0;
    double h=((double) (t_finish - t_start)/N);

    for (int i=0; i<N; i++) {
        grid_points[i] = t_start + i*h;
    }
    double *r = new double[N];
    double *r1 = new double[N];
    double *x1 = new double[N];
    double *y1 = new double[N];
    double *x2 = new double[N];
    double *y2 = new double[N];
    double *Vx1 = new double[N];
    double *Vy1 = new double[N];
    double *Vx2 = new double[N];
    double *Vy2 = new double[N];
    double h1 =h/2;
    double h2 =h*h/2;
    double m1=1.9/2.0*0.001; //mass jupiter devided by mass sun
    double m2=3.0*(1e-6); //mass earth  devided by mass sun
    double a = 4.0*M_PI*M_PI;
    double r2= 5.2*5.2*5.2; // distace from sun to Juipiter quibed
    clock_t g_start, g_finish;
    x1[0]=startx;
    y1[0]=starty;
    x2[0]=startx;
    y2[0]=starty;
    Vx1[0]=startvx0;
    Vy1[0]=startvy0;
    Vx2[0]=startvx0;
    Vy2[0]=startvy0;

    g_start = clock();
    // VERLET TWO BODY START
    for (int i=0;i<N;i++){

        r[i] = sqrt((x1[i] - x2[i])*(x1[i] - x2[i]) +
                    (y1[i] - y2[i])*(y1[i] - y2[i]));
        r1[i] = r[i]*r[i]*r[i];
        x1[i+1] = x1[i] + h1*Vx1[i] +
                h2*(-a*(1+m1/r1[i])*x1[i] + a*x2[i]/r1[i]);
        x2[i+1] = x2[i] + h1*Vx2[i] +
                h2*(-a*(1/r2 + m2/r1[i])*x2[i] + a*x1[i]/r1[i]);
        r[i+1] = sqrt((x1[i+1] - x2[i+1])*(x1[i+1] - x2[i+1]) +
                (y1[i+1] - y2[i+1])*(y1[i+1] - y2[i+1]));
        r1[i+1] = r[i+1]*r[i+1]*r[i+1];

        Vx1[i+1]=Vx1[i] + h1*(-a*(1+m1/r1[i+1])*x1[i+1]+a*x2[i+1]/r1[i+1] -
                a*(1+m1/r1[i])*x1[i]+a*x2[i]/r1[i]);

        Vx2[i+1]=Vx2[i] + h1*(-a*(1/r2 + m2/r1[i+1])*x2[i+1] + a*x1[i+1]/r1[i+1] -
                a*(1/r2 + m2/r1[i])*x2[i] + a*x1[i]/r1[i]);

        y1[i+1] = y1[i] + h1*Vy1[i] +
                h2*(-a*(1+m1/r1[i])*y1[i] + a*y2[i]/r1[i]);
        y2[i+1] = y2[i] + h1*Vy2[i] +
                h2*(-a*(1/r2 + m2/r1[i])*y2[i] + a*y1[i]/r1[i]);

        Vy1[i+1]=Vy1[i] + h1*(-a*(1+m1/r1[i+1])*y1[i+1]+a*y2[i+1]/r1[i+1] -
                a*(1+m1/r1[i])*y1[i]+a*y2[i]/r1[i]);

        Vy2[i+1]=Vy2[i] + h1*(-a*(1/r2 + m2/r1[i+1])*y2[i+1] + a*y1[i+1]/r1[i+1] -
                a*(1/r2 + m2/r1[i])*y2[i] + a*y1[i]/r1[i]);


    }
    // VERLET TWO BODY END
    g_finish=clock();

    for (int i=0; i<N; i++) {
            cout << grid_points[i] << "," <<'Position for Erth' << x1[i] << "," << y1[i] << endl;
            cout << grid_points[i] << "," <<'Position for Jupiter' << x2[i] << "," << y2[i] << endl;
        }

    delete [] grid_points;
    delete [] x1;
    delete [] y1;
    delete [] Vx1;
    delete [] Vy1;
    delete [] x2;
    delete [] y2;
    delete [] Vx2;
    delete [] Vy2;

    cout << "Time of VerletTwoBody(velocity) " << ((double) (t_finish-t_start)/CLOCKS_PER_SEC) << endl;
}


int celestial::print()
{
    cout << "M = " << mass << "X0 = " << startx <<"Y0 = "<< starty <<"Vx,Vy = "<< startvx0 <<","<< startvy0 << endl;
}
