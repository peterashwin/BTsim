/*
 *particle objects are assigned to units.
 *If a unit is occupied, one or more particle objects will be associated with, depending on the occupied number in units
 */
#include "particle.h"

Particle::Particle()
{
    labeled=false;
    type=0;
}
Particle::~Particle(){}

Particle::Particle(const Particle & copyparticle){
    index=copyparticle.index;
    type=copyparticle.type;
    location=copyparticle.location;
    labeled=copyparticle.labeled;
}

Particle & Particle::operator =(const Particle & otherparticle){
    index=otherparticle.index;
    type=otherparticle.type;
    location=otherparticle.location;
    labeled=otherparticle.labeled;
}

void Particle::moveTo(const std::vector<int> & newLocation, const int newType){
    type=newType;
    location=newLocation;
}
