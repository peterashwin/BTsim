
#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>

//using namespace std;
class Particle
{
public:
    Particle();
    ~Particle();
    Particle(const Particle &);
    Particle & operator=(const  Particle &);
    int index; // which is unique for each particle and distinguishes with other particles
    int type;
    std::vector<int> location; // the location in the bundle
    bool labeled;

    void moveTo(const std::vector<int> & newLocation, const int newType);// move to newLocation and be in newType
};

#endif // PARTICLE_H
