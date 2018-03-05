/*
 *a track object is a matrix (of size width*length) of units in geometry. the xCoordinate sets the location for the first units in the track.
 *each unit in the track has an orientation depending on the minus and plus locations of the track.
 *near the plus end, the track can be assocated with inhomogeneous transition rates to be determined in the bundle object
 *
 */
#include <vector>
#include "track.h"
#include "parameter.h"
#include <algorithm>

using namespace std;
Track::Track(const int & Ini_length, const int & Ini_width)
    : length(Ini_length), width(Ini_width){
    initializeTrack(); openBoundary=flagopenboundary;xCoordinate=0;initializePlusMinusLocations();initializeInhomoPara();//orientation=0;
}

Track::Track(const Track &copyTrack){
    istate=copyTrack.istate;
    length=copyTrack.length;
    width=copyTrack.width;
    //orientation=copyTrack.orientation;
    xCoordinate=copyTrack.xCoordinate;
    openBoundary=copyTrack.openBoundary;
    minusEndLocations=copyTrack.minusEndLocations;
    plusEndLocations=copyTrack.plusEndLocations;
    dyneinNum=copyTrack.dyneinNum;
    inhomoDistances=copyTrack.inhomoDistances;
}

Track::Track(){

}

Track & Track::operator =(const Track & otherTrack){
    istate=otherTrack.istate;
    length=otherTrack.length;
    width=otherTrack.width;
    //orientation=otherTrack.orientation;
    xCoordinate=otherTrack.xCoordinate;
    openBoundary=otherTrack.openBoundary;
    minusEndLocations=otherTrack.minusEndLocations;
    plusEndLocations=otherTrack.plusEndLocations;
    dyneinNum=otherTrack.dyneinNum;
    inhomoDistances=otherTrack.inhomoDistances;
}
Track::~Track(){}

void Track::getLength(){
    getXCoordinate();
    length=max(minusEndLocations.back(),plusEndLocations.back())-xCoordinate;
}

void Track::getXCoordinate(){
    std::vector<int> temp;
    temp.insert(temp.end(),minusEndLocations.begin(),minusEndLocations.end());
    temp.insert(temp.end(),plusEndLocations.begin(),plusEndLocations.end());
    xCoordinate=*min_element(temp.begin(),temp.end());
    //xCoordinate=min(minusEndLocations[0],plusEndLocations[0]);
}

int Track::getOrientation(const int site){//site can be -1 or maxLength for boundary

      if (site<minusEndLocations[0]) return 1;
      else if (site>=minusEndLocations.back()) return 0;
      else return 2;// more than one minus ends.  not included yet

}

Unit* Track::getUnit(const int lane, const int site){
    return &(istate[lane][site]);
}


void Track::initializeInhomoPara(){
    dyneinNum.clear();
    inhomoDistances.clear();
    for (std::vector<int>::size_type i=0;i<plusEndLocations.size();i++){
        dyneinNum.insert(dyneinNum.end(),num_of_dynein);
        inhomoDistances.insert(inhomoDistances.end(),int(distancePlus/hs));
    }
}

void Track::initializePlusMinusLocations(){
    minusEndLocations.clear();
    minusEndLocations.insert(minusEndLocations.begin(),0);
    plusEndLocations.clear();
    plusEndLocations.insert(plusEndLocations.begin(),length);
}

void Track::initializeTrack(){
    vector<int> temp(num_of_type,0);
    Unit empty_unit(temp);
    vector<Unit> temp_lane(length,empty_unit);
    istate.resize(width,temp_lane);
    temp.clear();
    temp_lane.clear();
}




void Track::showAllTypes(){
    for (int i=0;i<width;i++){
        for (int j=0;j<length+xCoordinate;j++){
            if (j<xCoordinate){
                cout<<" (";
                for (vector<int>::size_type x=0;x<istate[0][0].occupiedNum.size();x++) cout<<"  ";
                cout<<")";
            }
            else {
                cout<<" (";
                for (vector<int>::size_type x=0;x<istate[0][0].occupiedNum.size();x++)
                    cout<<istate[i][j-xCoordinate].occupiedNum[x]<<" ";
                cout<<")";
            }
        }
        cout<<endl;
    }
}


void Track::showSingleType(int type){
    if (type<num_of_type){
        cout<<"the state of "<<type+1<<"th type on this track is shown as below:"<<endl;
        for (int i=0;i<width;i++){
            for (int j=0;j<length+xCoordinate;j++){
                if (j<xCoordinate)  cout<<" ";
                else
                    cout<<istate[i][j-xCoordinate].occupiedNum[type]<<" ";
            }
            cout<<endl;
        }
    }
    else {
        cout<<"wrong input of type index, which should be less than the number of types"<<endl;
    }
}





