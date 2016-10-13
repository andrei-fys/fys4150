#ifndef CELESTIAL_H
#define CELESTIAL_H

class celestial
{
private:

public:
    double mass, startx, starty, startvx0, startvy0;
    celestial(double M, double x0, double y0, double Vx0, double Vy0);
    int Euler(int, double, char);
    int Verlet(int, double, char);
    int VerletTwoBody(int, double, char);
    int print();

//signals:

//public slots:
};

#endif // CELESTIAL_H
