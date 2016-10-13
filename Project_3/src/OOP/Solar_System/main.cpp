#include <iostream>
#include <vector>
#include "planet.h"
#include <math.h>
#include <fstream>

using namespace std;

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
            double G = 4.0*M_PI*M_PI;
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
            planet->r[i] += planet->v[i] * dt;
            planet->v[i] += (planet->f[i] / planet->mass) * dt;
        }
    }
}

void integrateVerlet(vector<Planet*> bodies, double dt) {

    computeForces(bodies);

    for(Planet *planet : bodies) {
        for(int i=0; i<3; i++) {
            planet->old_a[i] = planet->f[i] / planet->mass;
            planet->r[i] += planet->v[i] * dt + ((dt * dt) / 2)* planet->f[i]/planet->mass;
        }
    }

    computeForces(bodies);

    for(Planet *planet : bodies) {
        for(int i=0; i<3; i++) {
            planet->v[i] += (((planet->f[i] / planet->mass) + planet->old_a[i]) * dt/2.0);
        }
    }
}


int main() {
    Planet *sun = new Planet(0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0);
    Planet *earth = new Planet(1, 0, 0, 0, 2*M_PI, 0, 0, 0, 0, 3e-6);
    Planet *jupiter = new Planet (5.2, 0, 0, 0, 0.88*M_PI, 0, 0, 0, 0, 0.95e-3);

    vector<Planet*> bodies;

    bodies.push_back(sun);
    bodies.push_back(earth);
    bodies.push_back(jupiter);

    int N = 10000;
    double T = 30.0;
    double dt = (double) (T / N);
    for(int i=0; i<N; i++) {
        integrateEuler(bodies, dt);
        //integrateVerlet(bodies, dt);
        ofstream ofile;
        ofile.open("earth", std::ios::app);
        ofile << earth->r[0] << ',' << earth->r[1] << ',' << earth->r[2] << endl;
        ofile.close();

        ofile.open("sun", std::ios::app);
        ofile << sun->r[0] << ',' << sun->r[1] << ',' << sun->r[2] << endl;
        ofile.close();

        ofile.open("jupiter", std::ios::app);
        ofile << jupiter->r[0] << ',' << jupiter->r[1] << ',' << jupiter->r[2] << endl;
        ofile.close();

    }

}
