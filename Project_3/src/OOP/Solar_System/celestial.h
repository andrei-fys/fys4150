#ifndef CELESTIAL_H
#define CELESTIAL_H


class celestial
{
private:

public:
    double mass, startx, starty, startvx0, startvy0;
    celestial(double M, double x0, double y0, double Vx0, double Vy0);

    //void setX(double x) { components[0] = x; }
    //void setY(double y) { components[1] = y; }
    //void setZ(double z) { components[2] = z; }


    int Euler(int, double, char);
    int print();
//signals:

//public slots:
};

#endif // CELESTIAL_H
