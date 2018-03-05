#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <vector>
#define flaginhomo true //false true--> introduce linear turning rate of kinesin-carried (change to be dynein carried) near the plus end
#define flagden true // calculation of density
#define flagcurrent false // calculation of current
#define flagrunlengthend false //calculation of run length from tip end
#define flagrunlengthbundle false //calculation of run length from middle 
#define flagtipsize false //calculation of tip size
#define flagunaccessible true //true : consdier bundle geometry restrictions which setts some lanes unaccesible for particles to tranpsort
#define flagopenboundary false // boundary condition: false - close boundary condition
#define flagcolocalizetip false // calcualtion for colocalization
#define flagpause false  // calculation of pause time
#define flagdirpro false //calcualte density according to moving direction
#define flageedis false //calcualte distance between neiboring particles
#define flageeprojection false //calcualte the projected state of particles along the bundle axis
#define flageedistriscore false //compare EE spatial along cell axis to an even distribution at instantaneous time

#define flagdefault false // default initialization
#define flagrandom true //for random initialization
#define flagtipend false //initial particles at bundle tip end


// k--kinesin carried (type-0); d--dynein carried(type-1)
const int num_of_run=1;
const int num_of_type=2;
const int num_of_track=2;//default
const int num_of_lane=13;//lanes of each track
const double ttran=1000;//initial simulation time before averaging;
const double tmax=ttran+1000;//maximum time for simulation
const double ratio=0.025;//
const double fMAP=0.0;// MAP obstacles concentrationn
const int LL=100; //length of hyphae in unit of um
const double hs=0.1;//step space in unit of um  
const int L=int(double(LL)/hs);  // num of sites of a track
const int num_of_particle=53;//53;//int(theta*double(L*num_of_lane*num_of_track)); // number of particles
const int B1=4;//start of bundle region for switching or inhomogneous turning
const int B2=96;//end of bundle region for switching or inhomogneous turning

//rates
const double v_k=2.2;//velocity for kinesin carried
const double v_d=2.2;//velocity for dynein carried
const double M_k=39;//free run length of kinesin carried
const double M_d=62.37;//free run length of dyneinn carried
const double delta=5.1;// dynein dancing rate

const double w_k=v_k/M_k;//direction-change rate
const double w_d=v_d/M_d;
const double lane_d=delta*v_d/2;//p_k*delta*hs/2;//delta*v_d/2;
const double lane_k=0;//p_k/2;//lane-change rate
const double q=0.7;//
const double se=0.3*0.9*v_d/hs;//se track switching rates for all allowable switching
const double p_mturn1=0.1*0.3*v_d/hs;//minus on its own track 
const double p_mturn2=p_mturn1;//minus on another track
const double p_pturn2=w_d;// rate of turning when running into a plus end of another MT
const double si=0.002;//track swith within bundle
const double sek=1.0*v_k/hs;//track switch for type-0 at ends
const double sed=v_d/hs;//track switch for type-1 at ends

const double rbreak=0;//track break rate
const double rheal=0;//track heal rate

//information of inhomogeneous density of dynein

const double distancePlus=hs;//in the same unit of L
const int num_of_dynein=8;//number of dynein accumulated at plus ends
const double p_load=0.06;//loading of EEs at plus ends

const double pausemin=0.25;//minimum time spent at same site to be considered as a pause
//function declare
double ran2(long *);
int gillespieIndex(const std::vector<double >& , const double & );

struct TransitionEvent{
    int typeToMove; //the type before a transition
    std::vector<int> loc1; // the location before a transition //location is a vector of size 3 (track, lane, site)
    std::vector<int> loc2; // the location after a transition
    int typeToBe; // the type after a transition
    double rateValue;  // the corresponding transition rate
};//site x_coordinate-1 and x_coordinate+length indicate the boundary rates


struct TransitionRates{
     std::vector<double> forward;
     std::vector<double> forwardEnd;// reduced forward rate when trackswitching is possible
     std::vector<double> turning;
     std::vector<double> laneChange;
     std::vector<double> trackSwitchAtEnd;// into; out of
     std::vector<double> trackSwitchAtInner;
     std::vector<std::vector<double> >boundaryIn;  //track, type
     std::vector<std::vector<double> >boundaryOut;
};

struct DialogPara{
    int trackindex;
    int width;
    std::vector<int> minus;
    std::vector<int> plus;
    std::vector<int> distance;
    std::vector<int> dynein;
};

int maxAllowance(const int &);

bool exclusion(const std::vector<int> &, const std::vector<int> &, const std::vector<int> &,const int);

int sum(std::vector<std::vector<std::vector<std::vector<int> > > > &);

bool checkBundleConnect(const std::vector<std::vector<int> > &);


bool checkBundleConnect(const std::vector<std::vector<int> > &);
void writeToFile1(const std::vector<int> & vectors, char * pathname, char *filename);
void writeToFile2(const std::vector<std::vector<double> > & vectors, char * pathname, char *filename);
void writeToFile3(const std::vector<std::vector<std::vector<double> > >& vectors, char * pathname, char *filename);
void evendistribution(std::vector<double> & distribution,const int n, const int N);
#endif
