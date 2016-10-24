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
            double P_energ = G * celestial1->mass * celestial2->mass / (r);
            double fx = F*dx;
            double fy = F*dy;
            double fz = F*dz;
            celestial1->f[0] += fx;
            celestial1->f[1] += fy;
            celestial1->f[2] += fz;
            celestial1->P = P_energ;

            celestial2->f[0] -= fx;
            celestial2->f[1] -= fy;
            celestial2->f[2] -= fz;
            celestial2->P = P_energ;
        }
    }
}


void computeEnergies(vector<Celestial*> bodies, double timePoint) {
    for(Celestial *celestial : bodies) {
        celestial->K = celestial->mass*0.5*(pow(celestial->v[0],2) + pow(celestial->v[1],2) + pow(celestial->v[2],2));
    }
    double K_sum, P_sum, E_tot;
    for(Celestial *celestial : bodies) {
        K_sum += celestial->K;
        P_sum += celestial->P;
        E_tot += K_sum + P_sum;
    }
    ofstream e_file;
    e_file.open("Energy", std::ios::app);
    e_file << timePoint <<","<< K_sum << "," << P_sum << "," << E_tot << "\n";
    e_file.close();

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
        celestial->writeMyCoordinates(); /* Every time step it writes additional
                                          * file with celestial as file name only
                                          * coordinates of chosen celestial
                                          * Usefull for python static plots.
                                          */
    }
}

int main() {

    double earth_radius_au = 4.3e-5; // not used in calulations. Fancy to have for Ovito vizualization

    // 11. October 2016, data provided by NASA

    Celestial *sun = new Celestial(3.573492533603374E-03, 3.382190685467812E-03, -1.599607179196645E-04,
                                  -1.961253767828543E-06,  6.848684614020077E-06,  3.982709603458882E-08,
                                   0, 0, 0, 1.0, "Sun", 109.3*earth_radius_au, -1);
    Celestial *earth = new Celestial(9.536270134019510E-01,  3.098418103295414E-01, -1.776444588972545E-04,
                                    -5.568672077931286E-03,  1.631182923688181E-02,  2.223611254304490E-08,
                                     0, 0, 0, 3e-6, "Earth", 1.0*earth_radius_au, 1);
    Celestial *jupiter = new Celestial (-5.430639955496683E+00, -4.249175327527568E-01,  1.232159332227928E-01,
                                        5.011810831526324E-04, -7.166398278883576E-03,  1.856716507962747E-05,
                                        0, 0, 0, 0.95e-2, "Jupiter", 10.97*earth_radius_au, 1);


    vector<Celestial*> bodies;

    bodies.push_back(sun);
    bodies.push_back(earth);
    bodies.push_back(jupiter);



    // Transfering NASA velocities to AU/year
    for(Celestial *celestial : bodies) {
        for(int i=0; i<3; i++) {
            celestial->v[i] *= 365.0 ;
        }
    }

    double dt = 0.00001;
    //int N = 1000000;
    double T = 10.0;
    int N = (int) (T / dt);
    string my_file = "solar_system.xyz";
    for(int i=0; i<N; i++) {
        //integrateEuler(bodies, dt);
        integrateVerlet(bodies, dt);
        //computeEnergies(bodies, i*dt); //UNIT tests here
        //computeAngularMomentum(bodies, i*dt); // not saved because of sun, futher research
        writeToFile(bodies, my_file); // animation file
    }

}
