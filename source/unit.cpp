/*
 *this unit class is the smallest entity for a bundle object. A unit object is represented by a vector of non negative intergers:
 *zero being empty or occupied by a postive number of particle objects
 */
#include <vector>
#include <numeric>
#include "unit.h"

Unit::Unit(const std::vector<int> & Ini_occupiedNum)
    : occupiedNum(Ini_occupiedNum){unblocked=true;typeNum=occupiedNum.size();
                                  }

Unit::Unit(const Unit & copyUnit){
    occupiedNum=copyUnit.occupiedNum;
    typeNum=copyUnit.typeNum;
    unblocked=copyUnit.unblocked;
}

Unit & Unit::operator =(const Unit & otherUnit){
    occupiedNum=otherUnit.occupiedNum;
    typeNum=otherUnit.typeNum;
    unblocked=otherUnit.unblocked;
}
Unit::Unit(){}
Unit::~Unit(){}


int Unit::getTotalOccupiedNum(){
    return std::accumulate(occupiedNum.begin(),occupiedNum.end(),0);
}

/* this is annother constructor
 Unit::Unit(const int Ini_typeNum):typeNum(Ini_typeNum){
    unblocked=true;setOccupiedNum();
}*/


/*void Unit::unitRemove(){
    for (int i=0;i<typeNum;i++) occupiedNum[i]=0;
    unblocked=false;
}*/

/*void Unit::getTypeNum(){
    typeNum=occupiedNum.size();
}*/
/*
void Unit::setOccupiedNum(){
    occupiedNum.resize(typeNum);
    for (int i=0;i<typeNum;i++) occupiedNum[i]=0;
}
*/


