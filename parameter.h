#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <vector>
#define flaginhomo true //false true--> introduce linear turning rate of kinesin-carried (change to be dynein carried) near the plus end
#define flagden true
#define flagcurrent false //true
#define flagrunlengthend false //true //false //true //false
#define flagrunlengthbundle false
#define flagtipsize false //true
#define flagunaccessible true //false //true
#define flagopenboundary false
#define flagcolocalizetip false
#define flagpause false
#define flagdirpro false //true
#define flageedis false //true // ee distance
#define flageeprojection false //true
#define flageedistriscore false //true

#define flagdefault false // particles at end of first track
#define flagrandom true //false //true //false //
#define flagtipend false //true //false //true


// k--kinesin carried (type-0); d--dynein carried(type-1)
const int num_of_run=1;
const int num_of_type=2;
const int num_of_track=2;
const int num_of_lane=13;//lanes of each track
const double ttran=0.0;//3000;//300.0;//10.0;
const double tmax=ttran+10;//0;//.00;//3010.0;
const double ratio=0.025;//0.025;//0.5;//1.0;
const double fMAP=0.0;
const int LL=100; //length of hyphae in unit of um
const double hs=0.1;//0.1;//0.1;//1; //step space in unit of um  require hs<1/5.1 to keep p_d>0
const int L=int(double(LL)/hs);  // num of sites of a track
const int num_of_particle=53;//53;//int(theta*double(L*num_of_lane*num_of_track)); // number of particles
const int B1=4;//0;//4;//bundle region for switching or inhomogneous turning
const int B2=96;//100;//96;
//rates
const double v_k=2.2;//2.2;//PNAS 1.7;///10.0;//velocity
const double v_d=2.2;//v_k;
const double M_k=39;//10000;//39;//21.5;///10.0;
const double M_d=62.37;//30;//62.37;//EMBO/10.0;scaled with mean 92um of hyphal length gives 67.8 run length
const double delta=5.1;///10.0; times/um

const double w_k=v_k/M_k;//direction-change rate
const double w_d=v_d/M_d;
const double lane_d=delta*v_d/2;//p_k*delta*hs/2;//delta*v_d/2;
const double lane_k=0;//p_k/2;//lane-change rate
const double q=0.7;//1;//0.7;//0.5*0.83;
const double se=0.3*0.9*v_d/hs;//0.5*0.17*v_d/hs;// (1-q)/q=13/66 gives q=0.83// se track switching rates for all allowable switching
const double p_mturn1=0.1*0.3*v_d/hs;//minus on its own track //1.6*se;//w_d;//p_mturn/se=21/13 gives p-mtrun=
const double p_mturn2=p_mturn1;//minus on another track
const double p_pturn2=w_d;// rate of turning when running into a plus end of another MT
const double si=0.002;//0.8;//0.02;// track swith within bundle
const double sek=1.0*v_k/hs;//track switch for type-0 at ends
const double sed=v_d/hs;//sek;//se+v_d/hs*q;//track switch for type-1 at ends

const double rbreak=0;//0.001; //per um per s
const double rheal=0;//0.1; // per s

//information of inhomogeneous density of dynein

const double distancePlus=hs;//in the same unit of L
const int num_of_dynein=8;//55*2+8*5+3*0//55;
// for point inhomogeneity, i.e., distancePlus=hs, turn rate= 2N_EE * pload -w_d
const double p_load=0.06;//0.06;//0.4;//2.7;//p_load hs/v~~0.05 PNAS

//const double fluxin=w_k*v_k/(hs*p_load);
const double pausemin=0.25;//hs*2.0/v_d;
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
