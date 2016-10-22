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
            double c = 63145; //c, in AU/year
            double l = pow((celestial2->r[1]*celestial2->v[2] - celestial2->r[2]*celestial2->v[1]),2) -
                       pow((celestial2->r[0]*celestial2->v[2] - celestial2->r[2]*celestial2->v[0]),2) +
                       pow((celestial2->r[0]*celestial2->v[1] - celestial2->r[1]*celestial2->v[0]),2);
            //double l = pow((dy*celestial2->v[2] - dz*celestial2->v[1]),2) -
            //           pow((dx*celestial2->v[2] - dz*celestial2->v[0]),2) +
            //           pow((dx*celestial2->v[1] - dy*celestial2->v[0]),2);
            //cout << l << endl;
            double F = (G * celestial1->mass * celestial2->mass / (r*r*r)) * (1.0 + (3.0*l)/(r*r*c*c));
            //double F = G * celestial1->mass * celestial2->mass / (r*r*r);
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

void computeAngularMomentum(vector<Celestial*> bodies, double timePoint)
{
    for(Celestial *celestial : bodies) {
            double V = sqrt(celestial->v[0]*celestial->v[0] + celestial->v[1]*celestial->v[1] + celestial->v[2]*celestial->v[2]);
            double R = sqrt(celestial->r[0]*celestial->r[0] + celestial->r[1]*celestial->r[1] + celestial->r[2]*celestial->r[2]);
            celestial->L = V*R*celestial->mass*celestial->spin;
            ofstream aaa_file;
            aaa_file.open("AWFUL", std::ios::app);
            aaa_file << celestial->body_name <<","<<timePoint<<","<<celestial->L<<","<<V<<","<<R<<"\n";
            aaa_file.close();
    }
    double L_sum;
    for(Celestial *celestial : bodies) {
        L_sum += celestial->L;
    }
    ofstream l_file;
    l_file.open("Angular", std::ios::app);
    l_file << timePoint <<","<< L_sum << "\n";
    l_file.close();

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

double periheliumDataSave(Celestial* body)
{
    double r_Mercury = sqrt(body->r[0]*body->r[0] + body->r[1]*body->r[1] + body->r[2]*body->r[2]);
    return r_Mercury;
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


    Celestial *sun = new Celestial(0, 0, 0,
                                   0, 0, 0,
                                   0, 0, 0, 1.0, "Sun", 109.3*earth_radius_au, -1);

    Celestial *mercury = new Celestial (0.3075,  0,     0,
                                        0,       12.44, 0,
                                        0,       0,     0, 1.2e-7, "Mercury", 0.3829*earth_radius_au, 1);

    vector<Celestial*> bodies;

    bodies.push_back(sun);
    bodies.push_back(mercury);

    // 1e12
    double T = 100.0;
    double dt = 1e-7;
    long N = T/dt;
/*
    double ** merc_coord = new double*[N];
    for (int i=0;i<N;i++){
        merc_coord[i] = new double[4];
    }

    string my_file = "system";
    for(long i=0; i<N; i++) {
        //integrateEuler(bodies, dt);
        integrateVerlet(bodies, dt);
        //computeEnergies(bodies, i*dt);
        //computeAngularMomentum(bodies, i*dt);
        writeToFile(bodies, my_file);
        merc_coord[i][0] = periheliumDataSave(mercury);
        merc_coord[i][1] = mercury->r[0];
        merc_coord[i][2] = mercury->r[1];
        merc_coord[i][3] = mercury->r[2];
    }

    vector<int> merc_per;
    for(int i=1; i<N-1; i++) {
        if ( merc_coord[i][0] < merc_coord[i-1][0] && merc_coord[i][0] < merc_coord[i+1][0] ){
            merc_per.push_back(i);
            cout << "INDEX  "<< i << endl;
        }
    }

    for (int i : merc_per) {
        cout << "ATAMAN    " << atan(merc_coord[i][2]/merc_coord[i][1]) << endl;
    }
    cout << merc_per.size() << endl;
    /*for (vector<int>::iterator it = merc_per.begin(); it != merc_per.end(); ++it) {
        cout << it << endl;
    }*/

    // Set some helper variables before we start the time integration.
    double thetaPrevious 	= 0;	// The perihelion angle of the previous time step.
    double theta 		= 0;	// The perihelion angle of the current time step.
    double thetaCurrent;

    double rPreviousPrevious 	= 0;	// Mercury-Sun-distance two times steps ago.
    double rPrevious   	 	= 0;	// Mercury-Sun-distance of the previous time step.
    double r 		 	= 0;	// Mercury-Sun-distance of the current time step.

    ofstream perihelion_file;

    for(long i=0; i<N; i++) {
        integrateVerlet(bodies, dt);

        double x = mercury->r[0] - sun->r[0];
        double y = mercury->r[1] - sun->r[1];

        thetaCurrent = atan2( y, x );

        double rCurrent = sqrt(pow((mercury->r[0] - sun->r[0]),2) + pow((mercury->r[1] - sun->r[1]),2));

        if ( rCurrent > rPrevious && rPrevious < rPreviousPrevious ) {

            // If we are perihelion, print angle (in radians) to terminal.
            cout << "Perihelion angle: " << thetaPrevious << endl;
            perihelion_file.open("perihelion", std::ios::app);
            perihelion_file << i*dt << "," << thetaPrevious << endl;
            perihelion_file.close();

            // Here you should also probably write it to file for later plotting or something.
        }
        rPreviousPrevious 	= rPrevious;
        rPrevious		= rCurrent;
        thetaPrevious		= thetaCurrent;
        perihelion_file.close();
    }

}
