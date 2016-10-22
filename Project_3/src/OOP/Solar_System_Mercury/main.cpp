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

    for(int i=0; (unsigned)i<bodies.size(); i++) {
        Celestial *celestial1 = bodies[i];
        for(int j=i+1; (unsigned)j<bodies.size(); j++) {
            Celestial *celestial2 = bodies[j];

            double dx = celestial2->r[0] - celestial1->r[0];
            double dy = celestial2->r[1] - celestial1->r[1];
            double dz = celestial2->r[2] - celestial1->r[2];
            double dr2 = dx*dx + dy*dy + dz*dz;
            double r = sqrt(dr2);
            double G = 4.0*M_PI*M_PI;
            double c = 63145; //c, in AU/year
            double l = pow((celestial2->r[1]*celestial2->v[2] - celestial2->r[2]*celestial2->v[1]),2) -
                       pow((celestial2->r[0]*celestial2->v[2] - celestial2->r[2]*celestial2->v[0]),2) +
                       pow((celestial2->r[0]*celestial2->v[1] - celestial2->r[1]*celestial2->v[0]),2);
            double F = (G * celestial1->mass * celestial2->mass / (r*r*r)) * (1.0 + (3.0*l)/(r*r*c*c));
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

void integrateVerlet(vector<Celestial*> bodies, double dt) {

    computeForces(bodies);

    for(Celestial *celestial : bodies) {  // Computing axeleration to use on next time step(ald_a attribute)
        for(int i=0; i<3; i++) {
            celestial->old_a[i] = celestial->f[i] / celestial->mass;
            celestial->r[i] += celestial->v[i] * dt + ((dt * dt) / 2)* celestial->f[i]/celestial->mass;
        }
    }

    computeForces(bodies);                // Getting actual axeleration

    for(Celestial *celestial : bodies) {  // Compute velocity using current and previous step axxelerations
        for(int i=0; i<3; i++) {
            celestial->v[i] += (((celestial->f[i] / celestial->mass) + celestial->old_a[i]) * dt/2.0);
        }
    }
}

int main() {
    double earth_radius_au = 4.3e-5; //Used just to for visualization in Ovito

    Celestial *sun = new Celestial(0, 0, 0,
                                   0, 0, 0,
                                   0, 0, 0, 1.0, "Sun", 109.3*earth_radius_au, -1);
    // Mercury starts at x-axis with 12.44 AU/year moving conterclockwise
    Celestial *mercury = new Celestial (0.3075,  0,     0,
                                        0,       12.44, 0,
                                        0,       0,     0, 1.2e-7, "Mercury", 0.3829*earth_radius_au, 1);

    vector<Celestial*> bodies;
    bodies.push_back(sun);
    bodies.push_back(mercury);

    double T = 100.0;
    double dt = 1e-8;
    long N = T/dt;

    // Set some helper variables before we start the time integration.
    double thetaPrevious 	= 0;	   // The perihelion angle of the previous time step.
    double thetaCurrent;

    double rPreviousPrevious 	= 0;	// Mercury-Sun-distance two times steps ago.
    double rPrevious   	 	    = 0;    // Mercury-Sun-distance of the previous time step.

    ofstream perihelion_file;

    for(long i=0; i<N; i++) {
        integrateVerlet(bodies, dt);
        double x = mercury->r[0] - sun->r[0];
        double y = mercury->r[1] - sun->r[1];
        thetaCurrent = atan2( y, x );
        double rCurrent = sqrt(pow((mercury->r[0] - sun->r[0]),2) + pow((mercury->r[1] - sun->r[1]),2));
        if ( rCurrent > rPrevious && rPrevious < rPreviousPrevious ) {
            // If we are perihelion, print angle (in radians) to terminal and write to file
            cout << "Perihelion angle: " << thetaPrevious << endl;
            perihelion_file.open("perihelion", std::ios::app);
            perihelion_file << i*dt << "," << thetaPrevious << endl;
            perihelion_file.close();
        }
        rPreviousPrevious 	= rPrevious;
        rPrevious           = rCurrent;
        thetaPrevious		= thetaCurrent;
        perihelion_file.close();
    }

}
