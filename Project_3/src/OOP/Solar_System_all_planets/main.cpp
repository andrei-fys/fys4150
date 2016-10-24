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
			/* Iterates over all bodies and calculates fores by coupling them */
			/* distance between selestials */
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
	/* first step is to compute axeleration components needed on nex step */    
	/* at the same time we compute new coordinates */ 
    for(Celestial *celestial : bodies) {
        for(int i=0; i<3; i++) {
            celestial->old_a[i] = celestial->f[i] / celestial->mass;
            celestial->r[i] += celestial->v[i] * dt + ((dt * dt) / 2)* celestial->f[i]/celestial->mass;
        }
    }
	/* new accelerations */
    computeForces(bodies);
	/* finally velocity */
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
    double earth_radius_au = 4.3e-5;

    /*
    * Yahooo! Escape velocity test is working. Looks like we implemented velocity Verlet correct
    *
    */
    //Celestial *sun = new Celestial(0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0, "Sun", 109.3*earth_radius_au, -1);
    //Celestial *earth = new Celestial(1, 0, 0, 0, 0.75*sqrt(2.0)*2*M_PI, 0, 0, 0, 0, 3e-6, "Earth", 1.0*earth_radius_au, 1);
	
	// 11. October 2016, data provided by NASA
    
	Celestial *sun = new Celestial(3.573492533603374E-03, 3.382190685467812E-03, -1.599607179196645E-04,
                                  -1.961253767828543E-06,  6.848684614020077E-06,  3.982709603458882E-08,
                                   0, 0, 0, 1.0, "Sun", 109.3*earth_radius_au, -1);
    Celestial *earth = new Celestial(9.536270134019510E-01,  3.098418103295414E-01, -1.776444588972545E-04,
                                    -5.568672077931286E-03,  1.631182923688181E-02,  2.223611254304490E-08,
                                     0, 0, 0, 3e-6, "Earth", 1.0*earth_radius_au, 1);
    Celestial *jupiter = new Celestial (-5.430639955496683E+00, -4.249175327527568E-01,  1.232159332227928E-01,
                                        5.011810831526324E-04, -7.166398278883576E-03,  1.856716507962747E-05,
                                        0, 0, 0, 0.95e-3, "Jupiter", 10.97*earth_radius_au, 1);
    Celestial *venus = new Celestial (1.225258696979552E-01, -7.140769027371650E-01, -1.686184874263122E-02,
                                      1.981659775036702E-02,  3.244564766878831E-03, -1.099221302378845E-03,
                                      0, 0, 0, 0.4e-6, "Venus", 0.9499*earth_radius_au,1);
    Celestial *saturn = new Celestial (-2.287558104465964E+00, -9.769680590273008E+00,  2.609098499564027E-01,
                                       5.125440738199292E-03, -1.288971396833438E-03, -1.817263069662650E-04,
                                       0, 0, 0, 2.75e-4, "Saturn", 9.14*earth_radius_au, 1);
    Celestial *uranus = new Celestial (1.846929716791574E+01,  7.547597393735306E+00, -2.112414878691660E-01,
                                       -1.516589563922049E-03,  3.457480606816019E-03,  3.250985990748335E-05,
                                       0, 0, 0, 4.4e-5, "Uranus", 3.981*earth_radius_au, 1);
    Celestial *mars = new Celestial (1.128390641775210E+00, -8.012991904292689E-01, -4.462796684198677E-02,
                                       8.671623920281516E-03,  1.258717477587190E-02,  5.080016699199911E-05,
                                       0, 0, 0, 3.3e-7, "Mars", 0.532*earth_radius_au, 1);
    Celestial *mercury = new Celestial (-2.939616533863970E-01,  1.796773137244946E-01,  4.154186207623273E-02,
                                        -2.011072661591013E-02, -2.301806775661566E-02, -3.655500038055793E-05,
                                       0, 0, 0, 1.2e-7, "Mercury", 0.3829*earth_radius_au, 1);
    Celestial *neptune = new Celestial (2.825685429079379E+01, -9.934214960318053E+00, -4.466312684684883E-01,
                                        1.020709773953188E-03,  2.980178063753966E-03, -8.515600087600944E-05,
                                       0, 0, 0, 0.15e-4, "Neptune", 3.865*earth_radius_au,1);

    Celestial *pluto = new Celestial (9.411553037487414E+00, -3.181920205260207E+01,  6.824607170720862E-01,
                                      3.066338725295883E-03,  2.292371257231461E-04, -9.168773634294448E-04,
                                       0, 0, 0, 0.655e-8, "Pluto", 0.186*earth_radius_au, 1);



    vector<Celestial*> bodies;

    bodies.push_back(sun);
    bodies.push_back(earth);
    bodies.push_back(jupiter);
    bodies.push_back(venus);
    bodies.push_back(saturn);
    bodies.push_back(uranus);
    bodies.push_back(mars);
    bodies.push_back(mercury);
    bodies.push_back(neptune);
    bodies.push_back(pluto);

    for(Celestial *celestial : bodies) {
        for(int i=0; i<3; i++) {
            celestial->v[i] *= 365.0 ;
        }
    }

    int N = 10000;
    double T = 1.0;
    double dt = (double) (T / N);
    string my_file = "system";
    for(int i=0; i<N; i++) {
        //integrateEuler(bodies, dt);
        integrateVerlet(bodies, dt);
        computeEnergies(bodies, i*dt);
        computeAngularMomentum(bodies, i*dt);
        writeToFile(bodies, my_file);
        //for(int i=0; i<3; i++) {
        //    cout << bodies[0]->v[0] << "," << bodies[0]->v[1] << "," << bodies[0]->v[2] << endl;
        //}
    }

}
