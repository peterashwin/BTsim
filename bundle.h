
#ifndef BUNDLE_H
#define BUNDLE_H

#include <vector>
#include "track.h"
#include "particle.h"
#include "parameter.h"
//using namespace std;// dangerous to use this in head files
class Bundle{

private:
    double rMin;
    void initializeBundle();

    //check a particle of one type (third input, currently only type=1 is defined) at leaving site (first input) of the track(second input) meet a minus or not
    bool meetMinus(const int site, const int track, const int type);

    bool meetPlus(const int site, const int track, const int type);
    // get a vector of neigh locations for a given site(first input) and a fixe radius (second input),
    void getNeigh(const int site,const int radius, std::vector<std::vector<int> > & neigh);

    // add all possible rates occur among the region centered at the given site and a radius one
    void addRate(const int site, const int track);
    double getAvailableLane(const int track, const int site);
    int getForwardSite(const int orientation,const int type);
public:
    Bundle(const int & Ini_tracknumber,const int & Ini_lanenumber, const int & Ini_sitenumber);
    ~Bundle();

    bool flagRunLengthEnd;
    bool flagRunLengthBundle;
    bool flagDen;
    bool flagCurrent;
    bool flagInhomo;
    bool flagTipSize;
    bool unaccessibility;//set unaccessible region as blockages
    bool flagColocalizeTip;
    bool flagPause;
    bool flagDirPro;
    bool flagEEDis;
    bool flagEEDistriScore;

    double labelPer;
    std::vector<int> labelRegion;
    int currentCounter;
    int maxLength;
    int maxWidth;
    int particleNum;
    int typeNum;// only set typeNum=2
    int trackNum;
    int laneNum;//only for initialize
    int siteNum;// only for initialize
    int blockedNum;
    int tipRegion;
    double pLoad;
    double stepSize;
    int bundleRegion_1,bundleRegion_2;
    //physical quantities
    std::vector<double> runLengthUnbundled;
    std::vector<double>runLengthBundled;
    std::vector<double> velocity;

    TransitionRates transitionRates;
    std::vector<Track> bundleState;
    std::vector<TransitionEvent> rate;
    std::vector<std::vector<int> > colocalization;
    std::vector<std::vector<std::vector<std::vector<int> > > >  nstate; //4D in [type, track, lane, site]
    std::vector<std::vector<std::vector<double> > >denSum; //3D in [type, lane, site]
    std::vector<std::vector<std::vector<double> > >currentSum;  //3D in[type lane site]
    std::vector<std::vector<double> > runLength;// vector of [particle index, location, time]
    std::vector<int> runLenghValues;// store run length in unit of sites of individual particles
    std::vector<std::vector<int> > vkymograph; //vecotr of coarsened state
    std::vector<std::vector<double> >tipSize; //vector of [steps, tipsize]
    std::vector<int> particleIndexes; // available index for new arrival particles; only used for open boundary
    std::vector<Particle> listOfParticles;// list of particles in the bundle in order
    std::vector<std::vector<double> > pauses;// list of pauses
    void initializeBundleVar();
    void addOneTrack(const Track & track);
    void removePause(const Particle & particle);
    void recordPause(const Particle & particle);
    void labelParticles(const int type, const int start, const int end);
    void clearLable();
    // get the 'bundled region' of the whole bundle: the region apart from the unbundled region at two ends;
    void getBundleRegion(std::vector<int> & region);
    Unit* getUnit(const int track, const int lane, const int site);
    Track* getTrack(const int track);
    bool isMinus(const int site, const int track);
    bool isPlus(const int site, const int track);
    bool isBundled(const int site);
    void bundledPosition(std::vector<int> & pp,const int site_1,const int site_2,const int track_1,const int track_2);// return 0 notbundle, 1 bundle boundary and 2 interior of bundle
    void getDimension();//set maxLength and maxWidth
    int getXCoordinate(const int track);
    int getTotalSites();
    void getKPercentage(std::vector<double>& per);//calculate the percentage of type 0 particles in the bundle
    void getDirectMove(std::vector<int> & num, const std::vector<int> & region);
    //add a new line in the vkymograph variable by the current state if ture, otherwise, add a line as the same as the front row in vkymograph
    void kymographRecord(const int time_size, bool current);

    // return the update index in the transition events
    int gillespieUpdate(long *idum, double & t);


    /*    particle related below */

    //get total number of particles in the bundle
    int getParticleNum();

    // reset particle label
    void setParticleLabel(Particle & partilec, bool & labeled);

    // give a new particle with given particle parameters
    void assignParticle(Particle & newParticle, int *particleParameters);

    //get an index for a new particle
    int getParticleIndex();

    int getHeight(const int track, const int lane);
    /*    state related  below       */

    // set the state in the bundle randomly with fixed number of particles
    void randomInitializeState(long *idum);
    void tipendInitializeState();// initial alll particles at the tip


    // set the state in the bundle with particles at the right hand side
    void oneendInitializeState();

    //set the state in the bundle with particles at the first lane of the first track, starting from the right side
    void defaultInitializeState();
    void sumColocalize();
    // add the current state to denSum variable;
    void sumState(const int n);//n rate index to occur
    void sumStateParticles();
    void sumCurrent(const int n);
    void updateRunLength(const std::vector<int> & loc_1, const std::vector<int> & loc_2, Particle & particle, int ty_to_move);
    void addRunLengthVale(Particle & particle);
    void removeRunLength(Particle & particle);
    void addRunLength(Particle & particle);

    void labelBundle();
    void labelTip(const int type, const double per);
    // Copy the state in the bundle into the member nstate for easy calculation in density
    void copyState();

    // set n blockages randomly along the bundle
    void setRandomBlockages(const int n, long *idum);
    void addRandomBlockages(const int n, long *idum);
    void setUnaccessibleRegion();
    void removeUnaccessibleRegion();
    void setUnaccessibleBlockage();
    void setRandomBundleGeometry(std::vector<std::vector<int> > & randBundle, long *idum,int n);


    void generatePMends(std::vector<int> & minus, std::vector<int> & plus, long * idum, int n);
    /*     rate related      */

    //calculate the smallest possible rate value rMin;
   // void minRate();

    //initialize possible transition events with default rates
    void initializeTransitionRates();

    //get the turning rate of type 0 on given site and given track
    double inhomoWk(const int site, const int track);

    //obtain the rate given the type and location before and after the transition
    double getRateValue(const int type_to_move, const int type_to_be, const std::vector<int> & loc_1, const std::vector<int> & loc_2);

    //update the rates for transitions occuring around the input site with input radius
    void updateRate(const int ind, const int radius);

    // initialize the rate variable
    void initializeRateVector();

    /*                 calculations below             */

    //calculate tip size in the current state
    int calculateTipSize();

    // store tipsize into a vector
    void recordTipSize();

    int getTipTrack();

    /*    display state and rates     */
    void showSingleType(const int type);// display the state of one type
    void showAllTypes();//display the state of all types
    void showAllRates();//display all the possible rates with corresponding transtion-- the rate variable
    void showSingleRate(const int);//display the specific rate in rate variable with a given index
};

#endif
