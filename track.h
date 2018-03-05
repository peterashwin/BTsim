#ifndef TRACK_H
#define TRACK_H

#include <vector>
#include "unit.h"


class Track{
public:
    Track(const int &, const int &);
    Track(const Track &);
    Track();
    Track & operator=(const Track &);
    std::vector< std::vector<Unit> > istate;
    int length;
    int width;
    //int orientation;  //0- minus end on the left hand side; 1 - minus end on the right side
    int xCoordinate;  //the location of the first site alongth the length
    bool openBoundary;  // allowing injection/exit or not
    std::vector<int> minusEndLocations;
    std::vector<int> plusEndLocations;


    // the following are parameters for inhomogeneous turning rates
     std::vector<int> dyneinNum;// for inhomogeneity turning rate near plus ends
     std::vector<int> inhomoDistances; //unit um; for inhomogeneity turning rate regions
     std::vector<std::vector<int> > unaccessibleRegion;

    ~Track();
    void initializeInhomoPara();
    void initializeTrack();
    void showSingleType(int type);// show state of track for a single given type
    void showAllTypes();// show state of the track
    void getLength();// get track length from plus and minus locations
    int getOrientation(const int site);
    Unit* getUnit(const int lane , const int site);
    void getXCoordinate();// get xCoordinate from plus and minus locations
private:
    void initializePlusMinusLocations();
};

#endif
