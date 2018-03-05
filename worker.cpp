/*
 *this class does the main work running in a bundle object.
 *Jobs include initialize the bundle, update loops, calculations of density and current
 */

#include "worker.h"
#include <time.h>
#include <algorithm>
#include <functional>   // std::plus
#include <fstream>
#include <algorithm>// for find_first_of in vector
#include <numeric>// for partial_sum and accumulate in vector
#include <math.h>       /* fabs */
extern long* pdum;
extern double t;
extern int step;
extern int kymSize;
extern int tStep;
bool compi (int a,int b) {return (a>=b);}
double integral_L1(double a,double b) {return (fabs(a)+fabs(b));}
double diff(double a, double b) {return fabs(a-b);}
using namespace std;

Worker::Worker(){
     initialzieWorkerVar();
     initializeSimpleTwoTrack();
     //initializeFourTrack();
}

Worker::~Worker(){

}

void Worker::binDensity(std::vector<std::vector<double> >  & bindensity,const int binsize){

}

void Worker::clearCalculateVar(){
    bundle->particleIndexes.clear();
    bundle->listOfParticles.clear();
    bundle->runLength.clear();
    bundle->colocalization.clear();
    bundle->pauses.clear();
    bundle->rate.clear();
    bundle->nstate.clear(); //4D in [type, track, lane, site]
    bundle->denSum.clear(); //3D in [type, lane, site]
    bundle->currentSum.clear();  //3D in[type lane site]
    bundle->runLength.clear();// vector of [particle index, location, time]
    bundle->runLenghValues.clear();// store run length in unit of sites of individual particles
    bundle->vkymograph.clear(); //vecotr of coarsened state
    bundle->tipSize.clear(); //vector of [steps, tipsize]
    density.clear();
    current.clear();
    accum_dis.clear();
    score_t.clear();
}

void Worker::initialzieWorkerVar(){
    running=false;
    saving=false;
    recordBundle=false;
    recordMT=false;
    recordKymograph=false;
    flagKymograph=false;
    densityCounter=0;
    blocking=false;
    ending=false;
    speed=ratio;
    previousStep=0;
    threeDEffect=false;
}

void Worker::initializeCalculate(){
    if (flagden || bundle->flagDen) {
        bundle->flagDen=true;
        initializeDensityCalculation();
    }
    if (flagrunlengthend || bundle->flagRunLengthEnd) {
        bundle->flagRunLengthEnd=true;
        int a=bundle->maxLength;
        bundle->labelRegion.clear();
        bundle->labelRegion.push_back(a-1);
        bundle->labelRegion.push_back(a);
        std::vector<int> plus,minus;
        int tr=bundle->getTipTrack();
        plus=bundle->getTrack(tr)->plusEndLocations;
        minus=bundle->getTrack(tr)->minusEndLocations;
        if (plus[plus.size()-1]==bundle->maxLength)
            bundle->labelParticles(1,a-1,a);
        else
            bundle->labelParticles(0,a-1,a);

    }
    if (flagtipsize || bundle->flagTipSize){
        bundle->flagTipSize=true;
    }
    if (flagcurrent || bundle->flagCurrent){
       bundle->flagCurrent=true;
       initializeCurrentCalculation();
    }
    if(flageedistriscore){
        bundle->flagEEDistriScore=true;
        evendistribution(accum_dis,bundle->maxLength,num_of_particle);
    }
}

void Worker::calculateEEDistance(vector<int> &dis, const vector<int> &region){
   dis.clear();
   vector<int> sitelist;
   for (vector<Particle>::size_type p=0;p<bundle->listOfParticles.size();p++){
       sitelist.insert(sitelist.end(),bundle->listOfParticles[p].location[2]);
   }
   sort(sitelist.begin(),sitelist.end());
   vector<int>::iterator it;
   it=find_first_of(sitelist.begin(), sitelist.end(),region.begin(),region.begin()+1, compi);
   sitelist.erase(sitelist.begin(),it);
   it=find_first_of(sitelist.begin(), sitelist.end(),region.end()-1,region.end(), compi);
   sitelist.erase(it,sitelist.end());
   dis.resize(sitelist.size()-1);
   transform(sitelist.begin()+1,sitelist.end(),sitelist.begin(),dis.begin(),std::minus<int>());
}

void Worker::calculateEEdistribution(vector<double> &distr, const vector<int> state){
    distr.clear();
    distr.resize(state.size());
    partial_sum(state.begin(),state.end(),distr.begin());
}

double Worker::calculateEvendistributionSore(const std::vector<int> state){
    std::vector<double>dis;
    calculateEEdistribution(dis,state);
    evendistribution(accum_dis,bundle->maxLength,num_of_particle);
    std::transform(dis.begin(),dis.end(),accum_dis.begin(), dis.begin(), diff);
    double score=*std::max_element(dis.begin(),dis.end());
   // for (int i=0;i<static_cast<int>(state.size());i++) qDebug()<<dis[i]<<" "<<accum_dis[i];
   // double score=std::accumulate(dis.begin(),dis.end(),0,integral_L1);
    //score=score*worker->bundle->stepSize;
    return score;
}

void Worker::inputParameterSetup(char *filename){// boundar rate are homogeneus in the set up at the moment
    bundle->velocity.clear();
    bundle->runLengthBundled.clear();
    bundle->runLengthUnbundled.clear();
    bundle->transitionRates.forward.clear();
    bundle->transitionRates.forwardEnd.clear();
    bundle->transitionRates.laneChange.clear();
    bundle->transitionRates.trackSwitchAtEnd.clear();
    bundle->transitionRates.trackSwitchAtInner.clear();
    bundle->transitionRates.turning.clear();
    bundle->transitionRates.boundaryIn.clear();
    bundle->transitionRates.boundaryOut.clear();

    ifstream myBundle;
    myBundle.open(filename);
    //char output[256];
    if (myBundle.is_open()){
    double rateIndex;
    myBundle >> bundle->pLoad;

    for (int i=0;i<2;i++){
        myBundle >> rateIndex;
        bundle->velocity.insert(bundle->velocity.end(),rateIndex);
    }
    for (int i=0;i<2;i++){
        myBundle >> rateIndex;
        bundle->runLengthUnbundled.insert(bundle->runLengthUnbundled.end(),rateIndex);
    }
    for (int i=0;i<2;i++){
        myBundle >> rateIndex;
        bundle->runLengthBundled.insert(bundle->runLengthBundled.end(),rateIndex);
    }
    for (int i=0;i<3;i++){// have included the cross minus end on the same track (which is set to be 0)
        myBundle >> rateIndex;
        bundle->transitionRates.forward.insert(bundle->transitionRates.forward.end(),rateIndex);
    }
    for (int i=0;i<5;i++){
        myBundle >> rateIndex;
        bundle->transitionRates.turning.insert(bundle->transitionRates.turning.end(),rateIndex);
    }
    for (int i=0;i<4;i++){
        myBundle >> rateIndex;
        bundle->transitionRates.laneChange.insert(bundle->transitionRates.laneChange.end(),rateIndex);
    }
    for (int i=0;i<12;i++){
        myBundle >> rateIndex;
        bundle->transitionRates.trackSwitchAtEnd.insert(bundle->transitionRates.trackSwitchAtEnd.end(),rateIndex);
    }
    for (int i=0;i<4;i++){
        myBundle >> rateIndex;
        bundle->transitionRates.trackSwitchAtInner.insert(bundle->transitionRates.trackSwitchAtInner.end(),rateIndex);
    }
    for (int i=0;i<4;i++){
        myBundle >> rateIndex;
        bundle->transitionRates.forwardEnd.insert(bundle->transitionRates.forwardEnd.end(),rateIndex);
    }

    myBundle >> rateIndex;
    vector<double> alpha_k(bundle->trackNum,rateIndex);
    myBundle >> rateIndex;
    vector<double> alpha_d(bundle->trackNum,rateIndex);
    myBundle >> rateIndex;
    vector<double> beta_k(bundle->trackNum,rateIndex);
    myBundle >> rateIndex;
    vector<double> beta_d(bundle->trackNum,rateIndex);
    bundle->transitionRates.boundaryIn.insert(bundle->transitionRates.boundaryIn.end(),alpha_k);
    bundle->transitionRates.boundaryIn.insert(bundle->transitionRates.boundaryIn.end(),alpha_d);
    bundle->transitionRates.boundaryOut.insert(bundle->transitionRates.boundaryOut.end(),beta_k);
    bundle->transitionRates.boundaryOut.insert(bundle->transitionRates.boundaryOut.end(),beta_d);


    //bundle->minRate();
    bundle->initializeRateVector();
    }
}

void Worker::exportGeometry(char *filename){
    FILE * fgeometry;
    fgeometry=fopen(filename,"w");

    //bundle structures
    fprintf(fgeometry,"%d \n",bundle->trackNum);

    for (int tr=0;tr<bundle->trackNum;tr++){
        fprintf(fgeometry,"%d %d\n",bundle->getTrack(tr)->minusEndLocations.size(),bundle->getTrack(tr)->plusEndLocations.size());
        fprintf(fgeometry,"%d ",bundle->getTrack(tr)->width);
        for (vector<double>::size_type m=0;m<bundle->getTrack(tr)->minusEndLocations.size();m++)
                fprintf(fgeometry,"%d ",bundle->getTrack(tr)->minusEndLocations[m]);
        for (vector<double>::size_type p=0;p<bundle->getTrack(tr)->plusEndLocations.size();p++)
                fprintf(fgeometry,"%d ",bundle->getTrack(tr)->plusEndLocations[p]);
        for (vector<double>::size_type p=0;p<bundle->getTrack(tr)->inhomoDistances.size();p++)
                fprintf(fgeometry,"%d ",bundle->getTrack(tr)->inhomoDistances[p]);
        for (vector<double>::size_type p=0;p<bundle->getTrack(tr)->dyneinNum.size();p++)
                fprintf(fgeometry,"%d ",bundle->getTrack(tr)->dyneinNum[p]);
        fprintf(fgeometry,"\n\n\n");//    s.append("\n\n\n");
    }
    fclose(fgeometry);
}

double Worker::inputGeometrySetUp(char *filename){
    initialzieWorkerVar();
    bundle->initializeBundleVar();
    clearCalculateVar();

    ifstream myBundle;
    myBundle.open(filename);
    //char output[256];
    if (myBundle.is_open()){
           myBundle >> bundle->trackNum;
           bundle->bundleState.clear();
           int tr=0;
            while(tr<bundle->trackNum){
                int m1,m2;
                myBundle >> m1;//number of minus ends;
                myBundle >> m2;//          plus ends

                int list2[m2*3+m1+1];
                int counter=0;
                while (counter < m2*3+m1+1) {
                        myBundle >> list2[counter++];
                        cout << list2[counter - 1] << endl;
                    }

                int len=*(max_element(list2+1,list2+m1+m2+1))-*(min_element(list2+1,list2+1+m1+m2));
                Track temp(len,list2[0]);//tempary track with length L and width list2[0].toInt()
                bundle->bundleState.insert(bundle->bundleState.end(),temp);

                bundle->getTrack(tr)->minusEndLocations.clear();
                bundle->getTrack(tr)->plusEndLocations.clear();
                bundle->getTrack(tr)->inhomoDistances.clear();
                bundle->getTrack(tr)->dyneinNum.clear();
                for (int i=0;i<m1;i++)
                    bundle->getTrack(tr)->minusEndLocations.insert(bundle->getTrack(tr)->minusEndLocations.end(),list2[i+1]);

                for (int j=0;j<m2;j++){
                    bundle->getTrack(tr)->plusEndLocations.insert(bundle->getTrack(tr)->plusEndLocations.end(),list2[j+m1+1]);
                    bundle->getTrack(tr)->inhomoDistances.insert(bundle->getTrack(tr)->inhomoDistances.end(),list2[j+m1+m2+1]);
                    bundle->getTrack(tr)->dyneinNum.insert(bundle->getTrack(tr)->dyneinNum.end(),list2[j+m2*2+m1+1]);
                }
                bundle->getTrack(tr)->getLength();
                bundle->getTrack(tr)->getXCoordinate();
                tr++;
             }
            bundle->getDimension();            
            bundle->tipRegion=int(0.85*double(bundle->maxLength));
    }
    int nMAP=int(fMAP*double(bundle->getTotalSites()));
   // cout<<nMAP;
    bundle->addRandomBlockages(nMAP,pdum);
    if(bundle->unaccessibility){
        bundle->setUnaccessibleRegion();
        bundle->setUnaccessibleBlockage();
    }
    if (flagdefault)
        bundle->defaultInitializeState();
    else if (flagrandom)
        bundle->randomInitializeState(pdum);
    else if (flagtipend)
        bundle->tipendInitializeState();
    else
        bundle->oneendInitializeState();

    //
    bundle->vkymograph.clear();
    bundle->kymographRecord(kymSize,true);
    //bundle->initializeRateVector();
    return double(LL)/double(bundle->maxLength);
}

void Worker::inputCompleteSetUp(char *filename){
    initialzieWorkerVar();
    bundle->initializeBundleVar();
    clearCalculateVar();

    ifstream myBundle;
    myBundle.open(filename);
    //char output[256];
    if (myBundle.is_open()){
            myBundle >> bundle->trackNum;
            bundle->bundleState.clear();
           int tr=0;
            while(tr<bundle->trackNum){
                int m1,m2;
                myBundle >> m1;//number of minus ends;
                myBundle >> m2;//          plus ends
                int list2[m2*3+m1+1];
                int counter=0;
                while (counter < m2*3+m1+1) {
                        myBundle >> list2[counter++];
                        cout << list2[counter - 1] << endl;
                    }

                int len=*(max_element(list2+1,list2+m1+m2+1))-*(min_element(list2+1,list2+1+m1+m2));
                Track temp(len,list2[0]);
                bundle->bundleState.insert(bundle->bundleState.end(),temp);

                bundle->getTrack(tr)->minusEndLocations.clear();
                bundle->getTrack(tr)->plusEndLocations.clear();
                bundle->getTrack(tr)->inhomoDistances.clear();
                bundle->getTrack(tr)->dyneinNum.clear();
                for (int i=0;i<m1;i++)
                    bundle->getTrack(tr)->minusEndLocations.insert(bundle->getTrack(tr)->minusEndLocations.end(),list2[i+1]);

                for (int j=0;j<m2;j++){
                    bundle->getTrack(tr)->plusEndLocations.insert(bundle->getTrack(tr)->plusEndLocations.end(),list2[j+m1+1]);
                    bundle->getTrack(tr)->inhomoDistances.insert(bundle->getTrack(tr)->inhomoDistances.end(),list2[j+m1+m2+1]);
                    bundle->getTrack(tr)->dyneinNum.insert(bundle->getTrack(tr)->dyneinNum.end(),list2[j+m2*2+m1+1]);
                }
                bundle->getTrack(tr)->getLength();
                bundle->getTrack(tr)->getXCoordinate();
                tr++;
             }
            bundle->getDimension();
            bundle->tipRegion=int(0.85*double(bundle->maxLength));
           //  cout<<bundle->maxLength<<" "<<bundle->maxWidth;

             // bundle state
             int parNum=0;
             int maxL=bundle->maxLength;
             for (int tr=0;tr<bundle->trackNum;tr++){
                 int xco=bundle->getTrack(tr)->xCoordinate;
                 int len=bundle->getTrack(tr)->length;
                 for (int la=0;la<bundle->getTrack(tr)->width;la++){
                     int list0[maxL];
                     int counter=0;
                     while (counter < maxL) {
                             myBundle >> list0[counter++];
                     }
                     int list1[maxL];
                     counter=0;
                     while (counter < maxL) {
                         myBundle >> list1[counter++];
                     }
                     for (int x=0;x<len;x++){
                             int occu0=list0[x+xco];
                             int occu1=list1[x+xco];
                             if (occu0==1 && occu1==0){
                                 bundle->bundleState[tr].istate[la][x].occupiedNum[0]=1;
                                 bundle->bundleState[tr].istate[la][x].occupiedNum[1]=0;
                                 Particle particle;
                                 bool labeled=false;
                                 int a[6]={parNum,labeled,tr,la,x+xco,0};
                                 bundle->assignParticle(particle,a);                                 
                                 bundle->listOfParticles.insert(bundle->listOfParticles.end(),particle);
                                 bundle->bundleState[tr].istate[la][x].occupiedParticle.clear();
                                 bundle->bundleState[tr].istate[la][x].occupiedParticle.push_back(particle);
                                 parNum++;
                             }
                             else if(occu0==0 && occu1==1){
                                 bundle->bundleState[tr].istate[la][x].occupiedNum[0]=0;
                                 bundle->bundleState[tr].istate[la][x].occupiedNum[1]=1;
                                 Particle particle;
                                 bool labeled=false;
                                 int a[6]={parNum,labeled,tr,la,x+xco,1};
                                 bundle->assignParticle(particle,a);                                 
                                 bundle->listOfParticles.insert(bundle->listOfParticles.end(),particle);
                                 bundle->bundleState[tr].istate[la][x].occupiedParticle.clear();
                                 bundle->bundleState[tr].istate[la][x].occupiedParticle.push_back(particle);
                                 parNum++;
                             }
                             else if(occu0==0 && occu1==0){
                                 bundle->bundleState[tr].istate[la][x].occupiedNum[0]=0;
                                 bundle->bundleState[tr].istate[la][x].occupiedNum[1]=0;
                                 bundle->bundleState[tr].istate[la][x].occupiedParticle.clear();
                             }
                             else if (occu0==1 && occu1==1){
                                 cout<<"bundle state does not satisfy the exclusion principle at site !";

                             }
                             else if(occu0>1 || occu1>1){//blocked
                                 bundle->getTrack(tr)->istate[la][x].unblocked=false;
                                 bundle->getTrack(tr)->istate[la][x].occupiedNum.clear();
                                 bundle->getTrack(tr)->istate[la][x].occupiedParticle.clear();
                             }
                       }
                 }
             }
             bundle->copyState();
             bundle->particleNum=parNum;

             // transition rates
            double rateIndex;
            myBundle >> bundle->pLoad;

            for (int i=0;i<2;i++){
                myBundle >> rateIndex;
                bundle->velocity.insert(bundle->velocity.end(),rateIndex);
            }
            for (int i=0;i<2;i++){
                myBundle >> rateIndex;
                bundle->runLengthUnbundled.insert(bundle->runLengthUnbundled.end(),rateIndex);
            }
            for (int i=0;i<2;i++){
                myBundle >> rateIndex;
                bundle->runLengthBundled.insert(bundle->runLengthBundled.end(),rateIndex);
            }
            for (int i=0;i<3;i++){// have included the cross minus end on the same track (which is set to be 0)
                myBundle >> rateIndex;
                bundle->transitionRates.forward.insert(bundle->transitionRates.forward.end(),rateIndex);
            }
            for (int i=0;i<5;i++){
                myBundle >> rateIndex;
                bundle->transitionRates.turning.insert(bundle->transitionRates.turning.end(),rateIndex);
            }
            for (int i=0;i<4;i++){
                myBundle >> rateIndex;
                bundle->transitionRates.laneChange.insert(bundle->transitionRates.laneChange.end(),rateIndex);
            }
            for (int i=0;i<12;i++){
                myBundle >> rateIndex;
                bundle->transitionRates.trackSwitchAtEnd.insert(bundle->transitionRates.trackSwitchAtEnd.end(),rateIndex);
            }
            for (int i=0;i<4;i++){
                myBundle >> rateIndex;
                bundle->transitionRates.trackSwitchAtInner.insert(bundle->transitionRates.trackSwitchAtInner.end(),rateIndex);
            }
            for (int i=0;i<4;i++){
                myBundle >> rateIndex;
                bundle->transitionRates.forwardEnd.insert(bundle->transitionRates.forwardEnd.end(),rateIndex);
            }

            myBundle >> rateIndex;
            vector<double> alpha_k(bundle->trackNum,rateIndex);
            myBundle >> rateIndex;
            vector<double> alpha_d(bundle->trackNum,rateIndex);
            myBundle >> rateIndex;
            vector<double> beta_k(bundle->trackNum,rateIndex);
            myBundle >> rateIndex;
            vector<double> beta_d(bundle->trackNum,rateIndex);

            bundle->transitionRates.boundaryIn.insert(bundle->transitionRates.boundaryIn.end(),alpha_k);
            bundle->transitionRates.boundaryIn.insert(bundle->transitionRates.boundaryIn.end(),alpha_d);
            bundle->transitionRates.boundaryOut.insert(bundle->transitionRates.boundaryOut.end(),beta_k);
            bundle->transitionRates.boundaryOut.insert(bundle->transitionRates.boundaryOut.end(),beta_d);

            //bundle->minRate();
            bundle->initializeRateVector();
    }
    bundle->vkymograph.clear();
    bundle->kymographRecord(kymSize,true);

}

void Worker::initializeOneTrack(){
    bundle=new Bundle(1,num_of_lane,L);
    bundle->getTrack(0)->openBoundary=true;
    bundle->getDimension();
    if(bundle->unaccessibility){
        bundle->setUnaccessibleRegion();
        bundle->setUnaccessibleBlockage();
    }
    bundle->defaultInitializeState();
    bundle->setRandomBlockages(bundle->blockedNum,pdum);
    bundle->kymographRecord(kymSize,true);
    bundle->initializeRateVector();
}

void Worker::initializeTwoTrack(){
    bundle=new Bundle(num_of_track,num_of_lane,9*L/10);
    bundle->bundleState[0].minusEndLocations.clear();
    bundle->bundleState[0].plusEndLocations.clear();
    bundle->bundleState[0].minusEndLocations.push_back(3*L/4);
    bundle->bundleState[0].plusEndLocations.push_back(L/10);
    bundle->bundleState[0].plusEndLocations.push_back(L);
    bundle->bundleState[1].minusEndLocations.clear();
    bundle->bundleState[1].plusEndLocations.clear();
    bundle->bundleState[1].minusEndLocations.push_back(L/4);
    bundle->bundleState[1].plusEndLocations.push_back(0);
    bundle->bundleState[1].plusEndLocations.push_back(9*L/10);
    sort(bundle->getTrack(0)->plusEndLocations.begin(),bundle->getTrack(0)->plusEndLocations.end());
    sort(bundle->getTrack(1)->plusEndLocations.begin(),bundle->getTrack(1)->plusEndLocations.end());
    // make sure each track is of length 9*L/10;

    for (int i=0;i<bundle->trackNum;i++) {
        bundle->getXCoordinate(i);
        bundle->bundleState[i].getLength();
        bundle->getTrack(i)->initializeInhomoPara();
    }
    bundle->getDimension();
    int nMAP=int(fMAP*double(bundle->getTotalSites()));
    bundle->addRandomBlockages(nMAP,pdum);
    if(bundle->unaccessibility){
        bundle->setUnaccessibleRegion();
        bundle->setUnaccessibleBlockage();
    }
    bundle->defaultInitializeState();

    bundle->kymographRecord(kymSize,true);
    bundle->initializeRateVector();
    vector<double> line(bundle->maxLength,0);
    vector<vector<double> > matrix(bundle->maxWidth,line);
    density.resize(bundle->typeNum,matrix);
    density.resize(bundle->typeNum);
}

void Worker::initializeSimpleTwoTrack(){
    int len=int(9*L/10);
    bundle=new Bundle(num_of_track,num_of_lane,len);
    bundle->bundleState[0].minusEndLocations.clear();
    bundle->bundleState[0].plusEndLocations.clear();
    bundle->bundleState[0].minusEndLocations.push_back(L-len);
    bundle->bundleState[0].plusEndLocations.push_back(L);
    bundle->bundleState[1].minusEndLocations.clear();
    bundle->bundleState[1].plusEndLocations.clear();
    bundle->bundleState[1].plusEndLocations.push_back(0);
    bundle->bundleState[1].minusEndLocations.push_back(len);

    sort(bundle->getTrack(0)->plusEndLocations.begin(),bundle->getTrack(0)->plusEndLocations.end());
    sort(bundle->getTrack(1)->plusEndLocations.begin(),bundle->getTrack(1)->plusEndLocations.end());
    for (int i=0;i<bundle->trackNum;i++){
        bundle->getXCoordinate(i);
        bundle->bundleState[i].getLength();
        bundle->getTrack(i)->initializeInhomoPara();
    }
    bundle->getDimension();
    bundle->tipRegion=int(0.85*double(bundle->maxLength));
    int nMAP=int(fMAP*double(bundle->getTotalSites()));
   // cout<<nMAP;
    bundle->addRandomBlockages(nMAP,pdum);
    if(bundle->unaccessibility){
        bundle->setUnaccessibleRegion();
        bundle->setUnaccessibleBlockage();
    }
    if (flagdefault)   bundle->defaultInitializeState();
    else      bundle->randomInitializeState(pdum);
    bundle->kymographRecord(kymSize,true);
    bundle->initializeRateVector();
    vector<double> line(bundle->maxLength,0);
    vector<vector<double> > matrix(bundle->maxWidth,line);
    density.resize(bundle->typeNum,matrix);
    density.resize(bundle->typeNum);
}

void Worker::initializeFourTrack(){

    int len=int(9*L/10);
    bundle=new Bundle(2,num_of_lane,len);
    bundle->bundleState[0].minusEndLocations.clear();
    bundle->bundleState[0].plusEndLocations.clear();
    bundle->bundleState[0].minusEndLocations.push_back(L-len);
   // bundle->bundleState[0].plusEndLocations.push_back(L/10);
    bundle->bundleState[0].plusEndLocations.push_back(L);
    bundle->bundleState[1].minusEndLocations.clear();
    bundle->bundleState[1].plusEndLocations.clear();
   // bundle->bundleState[1].minusEndLocations.push_back(L/4);
    bundle->bundleState[1].plusEndLocations.push_back(0);
    bundle->bundleState[1].minusEndLocations.push_back(len);


    Track temp2(int(0.3*double(L)),num_of_lane);
    temp2.minusEndLocations.clear();
    temp2.plusEndLocations.clear();
    temp2.minusEndLocations.push_back(int(0.35*double(L)));
    temp2.plusEndLocations.push_back(int(0.35*double(L))-temp2.length);
     bundle->bundleState.insert(bundle->bundleState.end(),temp2);
    //for (int i=0;i<4;i++){
    Track temp3(int(0.3*double(L)),num_of_lane);
    temp3.minusEndLocations.clear();
    temp3.plusEndLocations.clear();
    temp3.minusEndLocations.push_back(int(0.65*double(L)));
    temp3.plusEndLocations.push_back(int(0.65*double(L))+temp3.length);

    bundle->bundleState.insert(bundle->bundleState.end(),temp3);

    //}
    bundle->trackNum=bundle->bundleState.size();
    /*bundle->bundleState[2].minusEndLocations.clear();
    bundle->bundleState[2].plusEndLocations.clear();
    //bundle->bundleState[2].minusEndLocations.push_back(L/4);
    bundle->bundleState[2].minusEndLocations.push_back(int(0.35*double(L)));
    bundle->bundleState[2].plusEndLocations.push_back(int(0.05*double(L)));

    bundle->getTrack(3)->plusEndLocations.clear();
    bundle->getTrack(3)->minusEndLocations.clear();
    bundle->getTrack(3)->minusEndLocations.push_back(int(0.65*double(L)));
    bundle->getTrack(3)->plusEndLocations.push_back(int(0.95*double(L)));*/
    // make sure each track is of length 9*L/10;
   // for (int i=0;i<4;i++)
    //sort(bundle->getTrack(i)->plusEndLocations.begin(),bundle->getTrack(i)->plusEndLocations.end());
    //sort(bundle->getTrack(1)->plusEndLocations.begin(),bundle->getTrack(1)->plusEndLocations.end());
    for (int i=0;i<bundle->trackNum;i++){
        sort(bundle->getTrack(i)->plusEndLocations.begin(),bundle->getTrack(i)->plusEndLocations.end());
        bundle->getXCoordinate(i);
        bundle->bundleState[i].getLength();
        bundle->getTrack(i)->initializeInhomoPara();
    }


    bundle->getDimension();
    int nMAP=int(fMAP*double(bundle->getTotalSites()));
   // cout<<nMAP;
    bundle->addRandomBlockages(nMAP,pdum);

    if(bundle->unaccessibility){
        bundle->setUnaccessibleRegion();
        bundle->setUnaccessibleBlockage();
    }
    //bundle->randomInitializeState(pdum);
    bundle->defaultInitializeState();
    bundle->kymographRecord(kymSize,true);
    bundle->initializeRateVector();
    vector<double> line(bundle->maxLength,0);
    vector<vector<double> > matrix(bundle->maxWidth,line);
    density.resize(bundle->typeNum,matrix);
    density.resize(bundle->typeNum);
}

void Worker::oneStep(){
    double t0=t;
    tStep=int(t0/speed)+1;

    int n=bundle->gillespieUpdate(pdum,t); // n is the chosen index to update in the rate vector
    if (flagKymograph){
        if(t>=double(tStep)*speed && t<(double(tStep)+1)*speed){
            bundle->kymographRecord(kymSize,true);// add a new row from bundleState
            tStep++;
        }
        else if(t>=(double(tStep)+1)*speed){
            while(t>=(double(tStep)+1)*speed){
                bundle->kymographRecord(kymSize,false);// add a row equal to the previous row
                tStep++;
            }
            bundle->kymographRecord(kymSize,true);
            tStep++;
        }
    }
    if (bundle->flagDen==true) {
        bundle->sumState(n);
        //bundle->sumStateParticles();
        densityCounter++;
    }
   /* time(&end);
    long tt=end-start;
    cout<<tt;*/
    if (bundle->flagCurrent){
        bundle->sumCurrent(n);
        bundle->currentCounter++;
        //cout<<bundle->currentCounter<<endl;
    }
    if (bundle->flagTipSize && step%200==0 ){//&& step>200 ){
            bundle->recordTipSize();
    }
    if (bundle->flagColocalizeTip) {
        bundle->sumColocalize();
    }

    step++;

    // 1 is the radius where rates are updated; also, rate update should be after updating the current vector
    bundle->updateRate(n,1);
    return;
}



void Worker::nsteps(int nn){
   /*for(int j=0;j<nn;j++){
        oneStep();
    }*/
}

void Worker::returnState(){
    bundle->copyState();
}


void Worker::projectEE(std::vector<int> & state){
    state.clear();
    state.resize(bundle->maxLength);
    for (vector<Particle>::size_type p=0;p<bundle->listOfParticles.size();p++){
        int site=bundle->listOfParticles[p].location[2];
        state[site]++;
        //cout<<state[site];
    }
}

void Worker::initializeColocalizeCalculation(){
    bundle->colocalization.clear();
    bundle->colocalization.resize(2);
    for (int i=0;i<2;i++) {
         bundle->colocalization[i].resize(bundle->maxLength);
    }

}

void Worker::initializeDensityCalculation(){
    bundle->denSum.clear();
    density.clear();
    vector<double> temp_L(bundle->maxLength,0);
    vector<vector<double> > den_type(bundle->maxWidth,temp_L);
    density.resize(bundle->typeNum);
    for (int i=0;i<bundle->typeNum;i++) {
        bundle->denSum.insert(bundle->denSum.end(),den_type);
        density[i].resize(bundle->maxWidth);
        for (int j=0;j<bundle->maxWidth;j++) density[i][j].resize(bundle->maxLength);
    }
    temp_L.clear();
    den_type.clear();

    bundle->copyState();
    densityCounter=0;
}

void Worker::densityCalculate(){
    if(densityCounter>0){
      for (vector<vector<vector<double> > >::size_type ty=0;ty<density.size();ty++){
        for (vector<vector<double> >::size_type la=0;la<density[ty].size();la++)
            for (vector<double>::size_type si=0;si<density[ty][la].size();si++){
                density[ty][la][si]=bundle->denSum[ty][la][si]/densityCounter;
            }
      }
    }


}

void Worker::initializeCurrentCalculation(){
    bundle->currentSum.clear();
    current.clear();
    vector<double> temp_L(bundle->maxLength,0);
    vector<vector<double> > den_type(bundle->maxWidth,temp_L);
    current.resize(bundle->typeNum);
    for (int i=0;i<bundle->typeNum;i++) {
        bundle->currentSum.insert(bundle->currentSum.end(),den_type);
        current[i].resize(bundle->maxWidth);
        for (int j=0;j<bundle->maxWidth;j++) current[i][j].resize(bundle->maxLength);
    }
    temp_L.clear();
    den_type.clear();
    bundle->currentCounter=0;
}

void Worker::initialzePauseCalculation(){
    extern double t;
    bundle->pauses.clear();
    for (int i=0;i<static_cast<int>(bundle->listOfParticles.size());i++){
        if (bundle->listOfParticles[i].labeled){
            double a[2]={double(bundle->listOfParticles[i].index),t};
            std::vector<double> temp(a,a+2);
            bundle->pauses.insert(bundle->pauses.end(),temp);
        }
    }

}

void Worker::currentCalculate(){
    if (bundle->currentCounter>0)
    for (vector<vector<vector<double> > >::size_type ty=0;ty<current.size();ty++){
        for (vector<vector<double> >::size_type la=0;la<current[ty].size();la++)
            for (vector<double>::size_type si=0;si<current[ty][la].size();si++){
                current[ty][la][si]=bundle->currentSum[ty][la][si]/(tmax-ttran);///double(bundle->currentCounter);
            }
    }
}


//ind=0 sum in type; ind=1 sum in direction
void Worker::sumDenCurr(vector<vector<double> > & den, const vector<vector<vector<double> > >& vectors, int ind){
    den.clear();
    den.resize(vectors.size());

    if (ind==0){
        for (vector<vector<double> >::size_type ty=0;ty<den.size();ty++){
            den[ty].resize(vectors[ty][0].size());
            for (vector<vector<vector<double> > >::size_type y=0;y<vectors[ty].size();y++){
                transform(vectors[ty][y].begin(),vectors[ty][y].end(),den[ty].begin(),den[ty].begin(),plus<double>());
            }
        }
    }
    else if(ind==1){
        for (int i=0;i<2;i++)            den[i].resize(vectors[i][0].size());
        for (int ty=0;ty<bundle->typeNum;ty++){
            int y=0;
            for (int track=0;track<bundle->trackNum;track++){
                if (bundle->getTrack(track)->plusEndLocations.size()==1){
                    int ori=(bundle->getTrack(track)->plusEndLocations[0]<bundle->getTrack(track)->minusEndLocations[0]);
                    int dir=(ori==ty); //dir=0---> left otherwise right going
                    //if ori==0 and type=1 int dir=0;
                    //else if ori==1 and type =1 int dir=1;
                    for (vector<vector<int> >::size_type la=0;la<bundle->bundleState[track].istate.size();la++){
                        transform(vectors[ty][y].begin(),vectors[ty][y].end(),den[dir].begin(),den[dir].begin(),plus<double>());
                        y++;
                    }
                }
                else if(bundle->getTrack(track)->plusEndLocations.size()==2){
                    int m=bundle->getTrack(track)->minusEndLocations[0];
                    for (vector<vector<int> >::size_type la=0;la<bundle->bundleState[track].istate.size();la++){
                        int dir=(1==ty); //dir=0---> left otherwise right going
                        transform(vectors[ty][y].begin(),vectors[ty][y].begin()+m,den[dir].begin(),den[dir].begin(),plus<double>());
                        dir=(0==ty); //dir=0---> left otherwise right going
                        transform(vectors[ty][y].begin()+m,vectors[ty][y].end(),den[dir].begin()+m,den[dir].begin()+m,plus<double>());
                        y++;
                    }
                }
            }
        }
    }
}

void Worker::generateRandomBundles(){
    char filename[256],dosline[256];
    sprintf(dosline,"mkdir random_ABC");
    system(dosline);
    for (int i=0;i<50;i++){
        sprintf(filename,"random_ABC/eg7_test%d.txt",i);
        std::vector<std::vector<int> > randBundle;
        int nb=setRandomBundle(randBundle,pdum,3);// return number of MTs
        exportGeometry(filename);
    }
}

// output a random bundle with mdouble double minus ends which have two plus end on two sides and meantime reinitialize the bundle state
int Worker::setRandomBundle(std::vector<std::vector<int> > & randBundle, long *idum,int mdouble){
    randBundle.clear();
    while(randBundle.size()==0)   bundle->setRandomBundleGeometry(randBundle,idum,mdouble);
    int nb=randBundle[0].size();
        initialzieWorkerVar();
        bundle->initializeBundleVar();
        clearCalculateVar();

        std::vector<int> minus=randBundle[1];
        std::sort(minus.begin(),minus.end());
        std::vector<int>::iterator it;
        it = std::unique(minus.begin(), minus.end());
        bundle->trackNum=int(std::distance(minus.begin(),it));
       /* if (*max_element(minus.begin(),minus.end())*bundle->stepSize>90 || *min_element(minus.begin(),minus.end())*bundle->stepSize<10){
            cout<<"minus not right!";
        }*/
        std::vector<std::vector<int> > MTs;
        for (int i=0;i<bundle->trackNum;i++){
            std::vector<int> MT;
            MT.insert(MT.end(),randBundle[1][0]);//minus
            MT.insert(MT.end(),randBundle[0][0]);//plus
            randBundle[0].erase(randBundle[0].begin());
            randBundle[1].erase(randBundle[1].begin());
            it = find (randBundle[1].begin(), randBundle[1].end(), MT[0]);
            if (it!=randBundle[1].end()){
                int ind=it-randBundle[1].begin();
                MT.insert(MT.end(),randBundle[0][ind]);
                std::sort(MT.begin()+1,MT.end());
                randBundle[0].erase(randBundle[0].begin()+ind);
                randBundle[1].erase(randBundle[1].begin()+ind);
            }
            MTs.insert(MTs.end(),MT);
            MT.clear();
        }


        bundle->bundleState.clear();
        for(int tr=0;tr<bundle->trackNum;tr++){

                    int m2;
                    m2=MTs[tr].size()-1;//number of plus ends;
                    int len=abs(MTs[tr][1]-MTs[tr][0]);
                    if (m2>1) len=len+abs(MTs[tr][2]-MTs[tr][0]);
                    Track temp(len,num_of_lane);//tempary track with length L and width list2[0].toInt()
                    bundle->bundleState.insert(bundle->bundleState.end(),temp);

                    bundle->getTrack(tr)->minusEndLocations.clear();
                    bundle->getTrack(tr)->plusEndLocations.clear();
                    bundle->getTrack(tr)->inhomoDistances.clear();
                    bundle->getTrack(tr)->dyneinNum.clear();

                    bundle->getTrack(tr)->minusEndLocations.insert(bundle->getTrack(tr)->minusEndLocations.end(),MTs[tr][0]);

                    for (int j=0;j<m2;j++){
                        bundle->getTrack(tr)->plusEndLocations.insert(bundle->getTrack(tr)->plusEndLocations.end(),MTs[tr][j+1]);
                        bundle->getTrack(tr)->inhomoDistances.insert(bundle->getTrack(tr)->inhomoDistances.end(),1);
                        bundle->getTrack(tr)->dyneinNum.insert(bundle->getTrack(tr)->dyneinNum.end(),num_of_dynein);
                    }
                    bundle->getTrack(tr)->getLength();
                    bundle->getTrack(tr)->getXCoordinate();
         }
        bundle->getDimension();
        bundle->tipRegion=int(0.85*double(bundle->maxLength));
        int nMAP=int(fMAP*double(bundle->getTotalSites()));
        bundle->addRandomBlockages(nMAP,idum);
        if(bundle->unaccessibility){
            bundle->setUnaccessibleRegion();
            bundle->setUnaccessibleBlockage();
        }
        bundle->randomInitializeState(idum);
        //bundle->defaultInitializeState();
        bundle->vkymograph.clear();
        bundle->kymographRecord(kymSize,true);

        bundle->initializeTransitionRates();
        bundle->initializeRateVector();
        return nb;
}



void Worker::setOutBundleBlockage(){
    std::vector<int> temp;
    bundle->getBundleRegion(temp);
    for (int tr=0;tr<bundle->trackNum;tr++){
        int x1=bundle->getXCoordinate(tr);
        if (x1<temp[0]){
            for (int la=0;la<bundle->getTrack(tr)->width;la++)
                for (int i=0;i<temp[0]-x1;i++)
                    bundle->getTrack(tr)->istate[la][i].unblocked=false;
        }
        int x2=x1+bundle->getTrack(tr)->length;
        if (x2>temp[1]){
            for (int la=0;la<bundle->getTrack(tr)->width;la++)
                for (int i=temp[1]+1-x1;i<x2-x1;i++)
                    bundle->getTrack(tr)->istate[la][i].unblocked=false;
        }

    }


}
