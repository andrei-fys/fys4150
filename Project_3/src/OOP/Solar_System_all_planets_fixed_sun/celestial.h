#ifndef CELESTIAL_H
#define CELESTIAL_H
#include <string>
#include <fstream>

struct Celestial
{
    Celestial(double x0, double y0, double z0,
           double vx0, double vy0, double vz0,
           double fx0, double fy0, double fz0,
           double mass,
           std::string body_name,
           double radius,
           double c_spin);
    void writeMyCoordinates();

    double r[3];
    double v[3];
    double f[3];
    double old_a[3];
    double mass;
    std::string body_name;
    double body_radius;
    std::ofstream my_file;
    double K;
    double P;
    double L;
    double spin;

};

#endif // CELESTIAL_H
