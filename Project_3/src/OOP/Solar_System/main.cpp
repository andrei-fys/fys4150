#include <iostream>
#include <vector>
#include "celestial.h"
#include <math.h>
#include <fstream>
#include <string>
using namespace std;

void computeForces(vector<Celestial*> bodies) {
    // Reset all the forces
    for(Celestial *celestial : bodies) {
        celestial->f[0] = 0.0;
        celestial->f[1] = 0.0;
        celestial->f[2] = 0.0;
    }

    for(int i=0; i<bodies.size(); i++) {
        Celestial *celestial1 = bodies[i];
        for(int j=i+1; j<bodies.size(); j++) {
            Celestial *celestial2 = bodies[j];

            double dx = celestial2->r[0] - celestial1->r[0];
            double dy = celestial2->r[1] - celestial1->r[1];
            double dz = celestial2->r[2] - celestial1->r[2];
            double dr2 = dx*dx + dy*dy + dz*dz;
            double r = sqrt(dr2);
            double G = 4.0*M_PI*M_PI;
            double F = G * celestial1->mass * celestial2->mass / (r*r*r);
            double fx = F*dx;
            double fy = F*dy;
            double fz = F*dz;
            celestial1->f[0] += fx;
            celestial1->f[1] += fy;
            celestial1->f[2] += fz;

            celestial2->f[0] -= fx;
            celestial2->f[1] -= fy;
            celestial2->f[2] -= fz;
        }
    }
}

void integrateEuler(vector<Celestial*> bodies, double dt) {

    computeForces(bodies);
    for(Celestial *celestial : bodies) {
        for(int i=0; i<3; i++) {
            celestial->r[i] += celestial->v[i] * dt;
            celestial->v[i] += (celestial->f[i] / celestial->mass) * dt;
        }
    }
}

void integrateVerlet(vector<Celestial*> bodies, double dt) {

    computeForces(bodies);

    for(Celestial *celestial : bodies) {
        for(int i=0; i<3; i++) {
            celestial->old_a[i] = celestial->f[i] / celestial->mass;
            celestial->r[i] += celestial->v[i] * dt + ((dt * dt) / 2)* celestial->f[i]/celestial->mass;
        }
    }

    computeForces(bodies);

    for(Celestial *celestial : bodies) {
        for(int i=0; i<3; i++) {
            celestial->v[i] += (((celestial->f[i] / celestial->mass) + celestial->old_a[i]) * dt/2.0);
        }
    }
}

void writeToFile(vector<Celestial*> bodies, string filename)
{
    ofstream m_file;
    m_file.open(filename, std::ios::app);
    m_file << bodies.size() << "\n";
    m_file << "some comment" << "\n";
    for(Celestial *celestial : bodies) {

        m_file << celestial->body_radius << " " << celestial->body_name << " " << celestial->r[0] << " " << celestial->r[1] << " " << celestial->r[2] << "\n";
        celestial->writeMyCoordinates();
    }
}

int main() {
    double earth_radius_au = 4.3e-5;
    Celestial *sun = new Celestial(0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, "Sun", 109.3*earth_radius_au);
    Celestial *earth = new Celestial(1, 0, 0, 0, 2*M_PI, 0, 0, 0, 0, 3e-6, "Earth", 1.0*earth_radius_au);
    Celestial *jupiter = new Celestial (5.2, 0, 0, 0, 0.88*M_PI, 0, 0, 0, 0, 0.95e-3, "Jupiter", 10.97*earth_radius_au );
    //Celestial *venus = new Celestial (0.72, 0, 0, 0, 2.3*M_PI, 0, 0, 0, 0, 0.4e-6, "Venus", 0.9499*earth_radius_au);

    vector<Celestial*> bodies;

    bodies.push_back(sun);
    bodies.push_back(earth);
    bodies.push_back(jupiter);
    //bodies.push_back(venus);

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
