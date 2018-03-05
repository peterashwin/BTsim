#ifndef WORKER_H
#define WORKER_H

#include "bundle.h"
using namespace std;
class Worker
{
public:
    Worker();
    ~Worker();
    Bundle *bundle;
    double speed;// speed in showing kymograph
    int previousStep;//the step number in the previous kymograph recording
    double densityCounter;// for density calculation
    bool threeDEffect;//show 3D effect on the bundle
    bool blocking;//show blocked units in different colors
    bool ending;//show plus/minus ends in view
    bool running;// running the gillespie update
    bool saving; // for kymograph data saving
    bool recordBundle;// save images of bundles
    bool recordMT;//save images of MT
    bool recordKymograph;//save images of kymograph
    bool flagKymograph;

    std::vector<double> accum_dis;
    std::vector<std::vector<double> > score_t;

    vector<vector<vector<double> > > density;
    vector<vector<vector<double> > > current;

public:
    double calculateEvendistributionSore(const std::vector<int> state);
    void calculateEEdistribution(vector<double> & distr, const vector<int> state);
    void projectEE(std::vector<int> & state);
    void calculateEEDistance(std::vector<int> & dis, const std::vector<int> & region);
    void clearCalculateVar();
    void initialzieWorkerVar();
    void initializeCalculate();
    void inputParameterSetup(char *filename);
    double inputGeometrySetUp(char *filename);
    void exportGeometry(char *filename);
    void inputCompleteSetUp(char * filename);
    void initializeFourTrack();
    void initializeSimpleTwoTrack();
    void initializeTwoTrack();
    void initializeOneTrack();
    void oneStep();
    void nsteps(int);
    void returnState();
    void initialzePauseCalculation();
    void initializeColocalizeCalculation();
    void initializeDensityCalculation();
    void densityCalculate();
    void sumDenCurr(std::vector<std::vector<double> > & den, const std::vector<std::vector<std::vector<double> > > & vectors,int ind);//ind=0 sum in types; ind=1 sum in direction
    void initializeCurrentCalculation();
    void currentCalculate();
    void binDensity(std::vector<std::vector<double> > & bindensity, const int binsize);
    int setRandomBundle(std::vector<std::vector<int> > &randBundle, long *idum, int mdouble);
    void generateRandomBundles();
    void setOutBundleBlockage();
};


#endif // WORKER_H
