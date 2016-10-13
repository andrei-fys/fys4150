#ifndef PLANET_H
#define PLANET_H

struct Planet
{
    Planet(double x0, double y0, double z0,
           double vx0, double vy0, double vz0,
           double fx0, double fy0, double fz0,
           double mass);


    double r[3];
    double v[3];
    double f[3];
    double mass;

};

#endif // PLANET_H
