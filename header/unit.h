#ifndef UNIT_H
#define UNIT_H

#include <vector>
#include <iostream>
#include "particle.h"


/*
 This class defines a unit in a lattice. It is a single site in geometry and is used to accommodate particles
 */

class Unit{
public:
    Unit(const std::vector<int> &);//default constructor
    Unit(const Unit &);//copy
    Unit & operator=(const Unit &);//assign equal
    Unit();
    ~Unit();

    std::vector<int> occupiedNum;//occupied number for each type in the unit
    std::vector<Particle> occupiedParticle;// a list of particles in the unit
    int typeNum;// number of particle types
    int getTotalOccupiedNum();
    bool unblocked;  //true - able to move into; false - unable to move into

    //Unit(const int);
    //void unitRemove();
    //void setOccupiedNum();
   // void getTypeNum();

};

#endif

