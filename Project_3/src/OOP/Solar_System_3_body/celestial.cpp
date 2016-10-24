#include "celestial.h"
#include <string>
Celestial::Celestial(double x0, double y0, double z0,
               double vx0, double vy0, double vz0,
               double fx0, double fy0, double fz0,
               double mass0,
               std::string p_name,
               double radius,
               double c_spin)
{
    r[0] = x0;  /* Coordinates */
    r[1] = y0;
    r[2] = z0;

    v[0] = vx0;  /* Velocities */
    v[1] = vy0;
    v[2] = vz0;

    f[0] = fx0;  /* Forces */
    f[1] = fy0;
    f[2] = fz0;

    old_a[0] = 0.0;  /* Verlet solver axeleration components on previous step */
    old_a[1] = 0.0;
    old_a[2] = 0.0;

    mass = mass0;
    body_name = p_name;
    body_radius = radius;  // used just to play with Ovito vizualisation
    my_file.open(body_name);  // celestial body name is file name for coordinates
    spin = c_spin; // clockwise/counterclockwise
}

void Celestial::writeMyCoordinates()
{
    my_file << r[0] << "," << r[1] << std::endl;
}



