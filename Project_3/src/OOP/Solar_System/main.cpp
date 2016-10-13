#include <iostream>
#include <vector>
#include "planet.h"
#include <math.h>
using namespace std;

//int main()
//{
//    double M = 1.0;
//    double x_0=1.0;
//    double y_0=0.0;
//    double Vx_0=0.0;
//    double Vy_0=2.0*M_PI;
//    //Create object
//    celestial Earth(M, x_0, y_0, Vx_0, Vy_0);
//    //Create object
//    celestial Juipiter(M, x_0, y_0, Vx_0, Vy_0);


//    double t_max = 1.0;
//    int N = 100;
//    char outfile_verlet = 'verlet';
//    char outfile_euler = 'euler';
//    char outfile_verlet_two_body = 'verlet_two_body';
//    //Calls Verlet and Euler methods
//    //Earth.Verlet(N, t_max, outfile_verlet);
//    //Earth.Euler(N, t_max, outfile_euler);
//    //Earth.print();
//    Juipiter.VerletTwoBody(N, t_max, outfile_verlet_two_body);
//   // TwoBody.print();
//    return 0;
//}

void computeForces(vector<Planet*> bodies) {
    // Reset all the forces
    for(Planet *planet : bodies) {
        planet->f[0] = 0.0;
        planet->f[1] = 0.0;
        planet->f[2] = 0.0;
    }

    for(int i=0; i<bodies.size(); i++) {
        Planet *planet1 = bodies[i];
        for(int j=i+1; j<bodies.size(); j++) {
            Planet *planet2 = bodies[j];

            double dx = planet2->r[0] - planet1->r[0];
            double dy = planet2->r[1] - planet1->r[1];
            double dz = planet2->r[2] - planet1->r[2];
            double dr2 = dx*dx + dy*dy + dz*dz;
            double r = sqrt(dr2);
            double G = 4*M_PI*M_PI;
            double F = G * planet1->mass * planet2->mass / (r*r*r);
            double fx = F*dx;
            double fy = F*dy;
            double fz = F*dz;
            planet1->f[0] += fx;
            planet1->f[1] += fy;
            planet1->f[2] += fz;

            planet2->f[0] -= fx;
            planet2->f[1] -= fy;
            planet2->f[2] -= fz;
        }
    }
}

void integrateEuler(vector<Planet*> bodies, double dt) {
    computeForces(bodies);

    for(Planet *planet : bodies) {
        for(int i=0; i<3; i++) {
            planet->v[i] += planet->f[i] / planet->mass * dt;
            planet->r[i] += planet->v[i] * dt;
        }
    }
}

int main() {
    Planet *earth = new Planet(1, 0, 0, 0, 2*M_PI, 0, 0, 0, 0, 3e-6);
    Planet *sun = new Planet(0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0);
    Planet *jupiter = new Planet (5.2, 0, 0, 0, 0.88*M_PI, 0, 0, 0, 0, 0.95e-3);

    vector<Planet*> bodies;

    bodies.push_back(earth);
    bodies.push_back(sun);

    int N = 1000;
    double T = 1.0;
    double dt = T / N;

    for(int i=0; i<N; i++) {
        integrateEuler(bodies, dt);
    }

}
