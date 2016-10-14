#include <iostream>
#include <vector>
#include "planet.h"
#include <math.h>
#include <fstream>
#include <string>
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

void writeToFile(vector<Planet*> bodies, string filename)
{
    ofstream m_file;
    m_file.open(filename, std::ios::app);
    m_file << bodies.size() << "\n";
    m_file << "some comment" << "\n";
    for(Planet *planet : bodies) {

        m_file << planet->body_radius << " " << planet->body_name << " " << planet->r[0] << " " << planet->r[1] << " " << planet->r[2] << "\n";
    }
}

int main() {
    double earth_radius_au = 4.3e-5;
    Planet *sun = new Planet(0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, "Sun", 109.3*earth_radius_au);
    Planet *earth = new Planet(1, 0, 0, 0, 2*M_PI, 0, 0, 0, 0, 3e-6, "Earth", 1.0*earth_radius_au);
    Planet *jupiter = new Planet (5.2, 0, 0, 0, 0.88*M_PI, 0, 0, 0, 0, 0.95e-3, "Jupiter", 10.97*earth_radius_au );

    Planet *venus = new Planet (0.72, 0, 0, 0, 2.3*M_PI, 0, 0, 0, 0, 0.4e-6, "Venus", 0.9499*earth_radius_au);

    vector<Planet*> bodies;

    bodies.push_back(sun);
    bodies.push_back(earth);
    bodies.push_back(jupiter);
    bodies.push_back(venus);

    int N = 10000;
    double T = 30.0;
    double dt = (double) (T / N);
    string my_file = "system";
    for(int i=0; i<N; i++) {
        //integrateEuler(bodies, dt);
        integrateVerlet(bodies, dt);
        writeToFile(bodies, my_file);
    }

}
