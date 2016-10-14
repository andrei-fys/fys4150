#include "celestial.h"
#include <string>
Celestial::Celestial(double x0, double y0, double z0,
               double vx0, double vy0, double vz0,
               double fx0, double fy0, double fz0,
               double mass0,
               std::string p_name,
               double radius)
{
    r[0] = x0;
    r[1] = y0;
    r[2] = z0;

    v[0] = vx0;
    v[1] = vy0;
    v[2] = vz0;

    f[0] = fx0;
    f[1] = fy0;
    f[2] = fz0;

    old_a[0] = 0.0;
    old_a[1] = 0.0;
    old_a[2] = 0.0;

    mass = mass0;
    body_name = p_name;
    body_radius = radius;
    my_file.open(body_name);
}

void Celestial::writeMyCoordinates()
{
    my_file << r[0] << "," << r[1] << std::endl;
}



