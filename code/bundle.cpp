/*
 *the main routine in the bundle object is to do gillespie update according to the algorithm in the paper
 *http://pubs.acs.org/doi/abs/10.1021/j100540a008
 *
 */
#include <time.h>
#include <vector>
#include <numeric>//for accumulate
#include <algorithm>// for sort, min_element and transform in vector // std::find // std::random_shuffle
#include <math.h> //for log
#include <stdlib.h>     /* srand, rand */

#include "bundle.h"
#include "parameter.h"

extern int kymSize;
extern int step;
double op_densum(int i, double j) {return double((i>0))+j;}
double op_multiply(double i) {return i;}
using namespace std;
extern double t;


int myrandom (int i) { return std::rand()%i;}

bool less (int i) { return ((i%2)==1); }

Bundle::Bundle(const int & Ini_tracknumber,const int & Ini_lanenumber, const int & Ini_sitenumber)
    : trackNum(Ini_tracknumber),laneNum(Ini_lanenumber),siteNum(Ini_sitenumber) {
    initializeBundleVar();
    particleNum=num_of_particle;
    typeNum=num_of_type;
    stepSize=hs;
    bundleRegion_1=B1/hs;//bundle region in [B1,B2]
    bundleRegion_2=B2/hs;
    initializeBundle();
    getDimension();
    tipRegion=int(0.95*double(maxLength));
    copyState();
    initializeTransitionRates();
    rMin=0.0001;
    //minRate();
}

Bundle::~Bundle(){}


void Bundle::addOneTrack(const Track & track){

}

void Bundle::addRandomBlockages(const int n, long *idum){
    int N_temp=0;
    while(N_temp<n){
        int track,lane,site;
        double r=ran2(idum);
        track=int(double(bundleState.size())*r);
        //int xco=getXCoordinate(track);
        r=ran2(idum);
        lane=int(double(bundleState[track].istate.size())*r);
        r=ran2(idum);
        site=int(double(bundleState[track].istate[lane].size())*r);
        if (bundleState[track].istate[lane][site].unblocked==true){// && site+xco<9*maxLength/10){
            bundleState[track].istate[lane][site].unblocked=false;
            N_temp=N_temp+1;
        }
    }
}

void Bundle::addRate(const int site, const int track){
    int x_coor=bundleState[track].xCoordinate;
    double r;
    vector<int> loc_1(3,0),loc_2;
    Unit* target;
    Unit* target1; // the initial Unit
    vector<vector<int> > neigh;
    getNeigh(site,1,neigh);// 1 is the radius where rates are updated
    vector<int> type_temp;
    for (int lane=0;lane<static_cast<int>(bundleState[track].istate.size());lane++){
        loc_1[1]=lane;
        loc_1[0]=track;
        target1=getUnit(track,lane,site);

        //hopping
        type_temp=target1->occupiedNum;
        if  (accumulate(type_temp.begin(),type_temp.end(),0)>0){
            loc_1[2]=site;
            for (int ty=0;ty<typeNum;ty++){
                if (type_temp[ty]>0){
                    for (vector<vector<int> >::size_type x=0; x!=neigh.size(); ++x){
                        for (int ty_tobe=0;ty_tobe<typeNum;ty_tobe++){
                            target=getUnit(neigh[x][0],neigh[x][1],neigh[x][2]);
                            if(target->unblocked && exclusion(target->occupiedNum,loc_1,neigh[x],maxAllowance(ty_tobe)) && !(ty==ty_tobe && equal(loc_1.begin(),loc_1.end(),neigh[x].begin()))){
                                r=getRateValue(ty,ty_tobe,loc_1,neigh[x]);
                                if (r>rMin/2){
                                    TransitionEvent temp={ty,loc_1,neigh[x],ty_tobe,r};
                                    rate.insert(rate.end(),temp);
                                }
                            }
                        }
                    }
                }
            }
        }
        if (bundleState[track].openBoundary){
            if(site==x_coor){
                for (int ty_tobe=0;ty_tobe<typeNum;ty_tobe++){
                    //inject
                    if( target1->unblocked){
                        loc_1[2]=site-1;
                        loc_2=loc_1;
                        loc_2[2]=site;
                        if (exclusion(target1->occupiedNum,loc_1,loc_2,maxAllowance(ty_tobe)) ){
                            r=getRateValue(ty_tobe,ty_tobe,loc_1,loc_2);
                            if (r>rMin/2){
                                TransitionEvent temp={ty_tobe,loc_1,loc_2,ty_tobe,r};
                                rate.insert(rate.end(),temp);
                            }
                        }
                    }
                    //exit
                    if(target1->occupiedNum[ty_tobe]>0){
                        loc_1[2]=site;
                        loc_2=loc_1;
                        loc_2[2]=site-1;
                        r=getRateValue(ty_tobe,ty_tobe,loc_1,loc_2);
                        if (r>rMin/2){
                            TransitionEvent temp={ty_tobe,loc_1,loc_2,ty_tobe,r};
                            rate.insert(rate.end(),temp);
                        }
                    }
                }
            }
            else if (site==static_cast<int>(bundleState[track].istate[0].size())-1+x_coor){
                for (int ty_tobe=0;ty_tobe<typeNum;ty_tobe++){
                    double r;
                    //inject
                    if( target1->unblocked){
                        loc_1[2]=site+1;
                        loc_2=loc_1;
                        loc_2[2]=site;
                        if (exclusion(target1->occupiedNum,loc_1,loc_2,maxAllowance(ty_tobe)) ){
                            r=getRateValue(ty_tobe,ty_tobe,loc_1,loc_2);
                            if(r>rMin/2){
                                TransitionEvent temp={ty_tobe,loc_1,loc_2,ty_tobe,r};
                                rate.insert(rate.end(),temp);
                            }
                        }
                    }
                    //exit
                    if(target1->occupiedNum[ty_tobe]>0){
                        loc_1[2]=site;
                        loc_2=loc_1;
                        loc_2[2]=site+1;
                        r=getRateValue(ty_tobe,ty_tobe,loc_1,loc_2);
                        if(r>rMin/2){
                            TransitionEvent temp={ty_tobe,loc_1,loc_2,ty_tobe,r};
                            rate.insert(rate.end(),temp);
                        }
                    }
                }
            }
        }
    }
}

void Bundle::addRunLength(Particle &particle){
    particle.labeled=true;
    double a[]={double(particle.index),double(particle.location[2]),t};
    vector<double>temp(a,a+3);
    runLength.insert(runLength.end(),temp);
}

void Bundle::addRunLengthVale(Particle & particle){
       particle.labeled=false;
       int ind=0;
       while(ind<static_cast<int>(runLength.size()) && runLength[ind][0]!=particle.index){
           ind++;
       }
       if (int(runLength[ind][1])>particle.location[2])
           runLenghValues.insert(runLenghValues.end(),int(runLength[ind][1])-particle.location[2]);
       else if (int(runLength[ind][1])==particle.location[2]){
           cout<<" no run length";
       }
       else{
           cout<<"wrong paris of particles!";
           system("pause");
       }
       runLength.erase(runLength.begin()+ind);
}

void Bundle::assignParticle(Particle &newParticle, int *particleParameters){
    newParticle.index=particleParameters[0];
    newParticle.labeled=particleParameters[1];
    newParticle.location.assign(particleParameters+2,particleParameters+5);
    newParticle.type=particleParameters[5];
}

void Bundle::bundledPosition(std::vector<int> & pp, const int site_1,const int site_2,const int track_1,const int track_2){
    // we only can have one minus end at each track and one to two plus ends

    vector<int> temp;
    temp.insert(temp.end(),getTrack(track_1)->minusEndLocations.begin(),getTrack(track_1)->minusEndLocations.end());

    if (getTrack(track_1)->plusEndLocations.size()==1)
    temp.insert(temp.end(),getTrack(track_1)->plusEndLocations.begin(),getTrack(track_1)->plusEndLocations.end());
    else {
        std::vector<int> vplus=getTrack(track_1)->plusEndLocations;
        sort(vplus.begin(),vplus.end());
        if (site_1<getTrack(track_1)->minusEndLocations[0]) temp.insert(temp.end(),vplus[0]);
        else temp.insert(temp.end(),vplus[1]);
    }
    sort(temp.begin(),temp.end());

    std::vector<int> temp2;
    temp2.insert(temp2.end(),getTrack(track_2)->minusEndLocations.begin(),getTrack(track_2)->minusEndLocations.end());
    if (getTrack(track_2)->plusEndLocations.size()==1)
    temp2.insert(temp2.end(),getTrack(track_2)->plusEndLocations.begin(),getTrack(track_2)->plusEndLocations.end());
    else {
        std::vector<int> vplus=getTrack(track_2)->plusEndLocations;
        sort(vplus.begin(),vplus.end());
        if (site_2<getTrack(track_2)->minusEndLocations[0]) temp2.insert(temp2.end(),vplus[0]);
        else temp2.insert(temp2.end(),vplus[1]);
    }
    sort(temp2.begin(),temp2.end());

    if (temp[1]<=temp2[0]){
        pp[0]=0;
        pp[1]=0;
    }
    else{
        temp.insert(temp.end(),temp2.begin(),temp2.end());
        sort(temp.begin(),temp.end());
        if (site_1>=temp[1]+1 && site_1<temp[temp.size()-2]-1) pp[0]=2;// interior of the bundle
        else if(site_1==temp[1] || site_1==temp[temp.size()-2]-1) pp[0]=1; // at the bundary of the bundle
        else pp[0]=0;//outside the bundle

        if (site_2>=temp[1]+1 && site_2<temp[temp.size()-2]-1) pp[1]=2;// interior of the bundle
        else if(site_2==temp[1] || site_2==temp[temp.size()-2]-1) pp[1]=1; // at the bundary of the bundle
        else pp[1]=0;//outside the bundle
    }
}

void Bundle::clearLable(){
    for (int i=0;i<trackNum;i++){
        for (int la=0;la<getTrack(i)->width;la++){
            for (int si=0;si<getTrack(i)->length;si++){
                for (int p=0;p<static_cast<int>(getTrack(i)->istate[la][si].occupiedParticle.size());p++)
                    getTrack(i)->istate[la][si].occupiedParticle[p].labeled=false;
            }
        }
    }
}

int Bundle::calculateTipSize(){
    int tiptemp=0;
    for (int tr=0;tr<trackNum;tr++){
        int xco=getXCoordinate(tr);
        for (int si=std::max(tipRegion,xco);si<getTrack(tr)->length+xco;si++){
            for (int la=0;la<getTrack(tr)->width;la++){
                tiptemp=tiptemp+(getTrack(tr)->istate[la][si-xco].getTotalOccupiedNum()>0);
            }
        }
    }
    return tiptemp;
}


void Bundle::copyState(){
    nstate.clear();
    getDimension();
    for (int track=0;track<trackNum;track++){
        vector<vector<vector<int> > > temp_3;
        for (int ty=0;ty<typeNum;ty++){
            vector<vector<int> > temp_2;
            for (vector<vector<int> >::size_type y=0;y<bundleState[track].istate.size();y++){
                vector<int> temp;
                for (int i=0;i<bundleState[track].xCoordinate;i++) temp.insert(temp.end(),0);
                for (vector<int>::size_type x=0;x<bundleState[track].istate[y].size();x++)
                    temp.insert(temp.end(),bundleState[track].istate[y][x].occupiedNum[ty]);
                int m=temp.size();
                for (int i=0;i<maxLength-m;i++) temp.insert(temp.end(),0);
                temp_2.insert(temp_2.end(),temp);
                temp.clear();
            }
            temp_3.insert(temp_3.end(),temp_2);
            temp_2.clear();
        }
        nstate.insert(nstate.end(),temp_3);
        temp_3.clear();
    }
}


void Bundle::defaultInitializeState(){
    int N_temp=0;
    particleIndexes.clear();
    listOfParticles.clear();
    runLength.clear();
    colocalization.clear();
    pauses.clear();
    for (int track=0;track<trackNum;track++){
        for (int j=0;j<static_cast<int>(bundleState[track].istate.size());j++){
            for (int i=static_cast<int>(bundleState[track].istate[j].size())-1;i>=0;--i){
                if (N_temp<particleNum && bundleState[track].istate[j][i].unblocked){
                     bundleState[track].istate[j][i].occupiedNum[0]=1;
                     Particle particle;
                     bool labeled=false;
                     int a[6]={N_temp,labeled,track,j,i+getXCoordinate(track),0};
                     assignParticle(particle,a);
                     listOfParticles.insert(listOfParticles.end(),particle);
                     bundleState[track].istate[j][i].occupiedParticle.clear();
                     bundleState[track].istate[j][i].occupiedParticle.push_back(particle);
                     N_temp+=1;
                     if (typeNum>1){
                         for (int ty=1;ty<typeNum;ty++){
                             bundleState[track].istate[j][i].occupiedNum[ty]=0;
                         }
                     }
                }
                else {
                    for (int ty=0;ty<typeNum;ty++)
                        bundleState[track].istate[j][i].occupiedNum[ty]=0;
                    bundleState[track].istate[j][i].occupiedParticle.clear();
                }
            }
        }
    }
}

double Bundle::getAvailableLane(const int track, const int site){
     int xco=getXCoordinate(track);
     double m=0;
     for (int i=0;i<static_cast<int>(getTrack(track)->istate.size());i++){
         if (getTrack(track)->istate[i][site-xco].unblocked)  m++;
     }
     return m;
}


void Bundle::getBundleRegion(std::vector<int> &temp){
    temp.clear();
    for (int track=0;track<trackNum;track++){
    temp.insert(temp.end(),getTrack(track)->minusEndLocations.begin(),getTrack(track)->minusEndLocations.end());
    temp.insert(temp.end(),getTrack(track)->plusEndLocations.begin(),getTrack(track)->plusEndLocations.end());
    }
    sort(temp.begin(),temp.end());
    temp.erase(temp.begin());
    temp.erase(temp.end()-1);

    while(temp[0]<bundleRegion_1){
        temp.erase(temp.begin());
        temp.insert(temp.begin(),bundleRegion_1);
    }
    while(temp[temp.size()-1]>bundleRegion_2){
        temp.erase(temp.end()-1);
        temp.insert(temp.begin(),bundleRegion_2);
    }
    sort(temp.begin(),temp.end());
    temp.erase(temp.begin()+1,temp.end()-1);
}



void Bundle::getDimension(){
    int temp=0;
    for (int track=0;track<trackNum;track++){
        if (bundleState[track].length+bundleState[track].xCoordinate>temp)
            temp=bundleState[track].length+bundleState[track].xCoordinate;
    }
    maxLength=temp;
    int temp2=0;
    for (int track=0;track<trackNum;track++){
        temp2=temp2+bundleState[track].istate.size();
    }
    maxWidth=temp2;
}


int Bundle::getForwardSite(const int orientation, const int type){
    if (orientation==1) return 2*type-1;
    else return 1-2*type;
}

int Bundle::getHeight(const int track,const int lane){
    int w=0;
    for (int tr=0;tr<track;tr++){
        w=w+getTrack(tr)->width;
    }
    w=w+lane;
    return w;
}

int Bundle::getParticleIndex(){
    if (particleIndexes.empty())
        return particleNum;
    else{
       int a=particleIndexes[0];
       particleIndexes.erase(particleIndexes.begin());
       return a;
    }
}

void Bundle::getKPercentage(vector<double> & per){
    per.clear();
    per.resize(3);
    int k=0;
    int d=0;
    copyState();
    for (int tr=0;tr<static_cast<int>(nstate.size());tr++){
        for (int l=0;l<static_cast<int>(nstate[tr][0].size());l++){
            k=k+accumulate(nstate[tr][0][l].begin(),nstate[tr][0][l].end(),0);
            d=d+accumulate(nstate[tr][1][l].begin(),nstate[tr][1][l].end(),0);
        }
    }
    per[0]=double(k);
    per[1]=double(d);
    per[2]=double(k)/double(d+k)*100;
}

void Bundle::getDirectMove(std::vector<int> & num,const std::vector<int> & region){
    num.clear();
    num.resize(2);
    for (int p=0;p<listOfParticles.size();p++){
        std::vector<int> loc=listOfParticles[p].location;
        if (loc[2]>=region[0] && loc[2]<region[1]){
            int tr=loc[0];
            int ori=getTrack(tr)->getOrientation(loc[2]);
            int ty=listOfParticles[p].type;
            int dir=(ori==ty); //dir=0---> left otherwise right going
            num[dir]++;
        }
    }
}

void Bundle::getNeigh(const int site,const int radius, vector<vector<int> > & neigh){// site relative to the whole bundle
    neigh.clear();
    for (int tr=0; tr<trackNum;++tr){
        int length=bundleState[tr].length;
        int width=bundleState[tr].istate.size();
        for (int i=max(bundleState[tr].xCoordinate,site-radius);i<min(site+1+radius,length+bundleState[tr].xCoordinate);i++){
            for (int la=0;la<width;la++){
                int temp[3]={tr,la,i};
                vector<int> temp_v(temp,temp+3);
                neigh.insert(neigh.end(),temp_v);
            }
        }
    }
}

int Bundle::getParticleNum(){
    int n=0;
    for (vector<Track>::size_type tr=0;tr<bundleState.size();tr++){
        for (vector<vector<Unit> >::size_type la=0;la<bundleState[tr].istate.size();la++){
            for (vector<Unit>::size_type si=0;si<bundleState[tr].istate[la].size();si++){
                for (vector<int>::size_type ty=0;ty<bundleState[tr].istate[la][si].occupiedNum.size();ty++)
                    n=n+(bundleState[tr].istate[la][si].occupiedNum[ty]);
            }

        }
    }
    return n;
}

int Bundle::getTotalSites(){
    int n=0;
    for (vector<Track>::size_type tr=0;tr<bundleState.size();tr++){
        for (vector<vector<Unit> >::size_type la=0;la<bundleState[tr].istate.size();la++){
            n=n+bundleState[tr].istate[la].size();
        }
    }
    return n;
}

int Bundle::getTipTrack(){
    std::vector<int> plus,minus;
    int tipTrack;
    for (int tr=0;tr<trackNum;tr++){
        plus=getTrack(tr)->plusEndLocations;
        minus=getTrack(tr)->minusEndLocations;
        if (plus[plus.size()-1]==maxLength || minus[minus.size()-1]==maxLength){
            tipTrack=tr;
            break;
        }
    }
    return tipTrack;
}

Track* Bundle::getTrack(const int track){
    return &bundleState[track];
}


Unit* Bundle::getUnit(const int track, const int lane, const int site){
    int x_coor=bundleState[track].xCoordinate;
    Unit* p=&(bundleState[track].istate[lane][site-x_coor]);
    return p;
}

int Bundle::getXCoordinate(const int track){
    bundleState[track].getXCoordinate();
    return bundleState[track].xCoordinate;
}


int Bundle::gillespieUpdate(long * idum, double & t){
    if (rate.size()>0){
        vector<double> rate_v;
        for (vector<TransitionEvent>::iterator x=rate.begin();x!=rate.end();++x)
            rate_v.insert(rate_v.end(),x->rateValue);
        int n;
        double r,s;
        s=accumulate(rate_v.begin(),rate_v.end(),0.0);
        r=ran2(idum);
        t=t-log(r)/s;
        r=ran2(idum);
        n=gillespieIndex(rate_v,r);

        // do the transtion accordingly
        vector<int> temp_1=rate[n].loc1;
        vector<int> temp_2=rate[n].loc2;
        int x_coor1=bundleState[temp_1[0]].xCoordinate;
        int x_coor2=bundleState[temp_2[0]].xCoordinate;
        int ty_to_move=rate[n].typeToMove;
        int ty_to_be=rate[n].typeToBe;
        if (temp_1[2]==x_coor1-1 || temp_1[2]==getTrack(temp_1[0])->length+x_coor1){
            Unit* C2=getUnit(temp_2[0],temp_2[1],temp_2[2]); //injec
            C2->occupiedNum[ty_to_be]=C2->occupiedNum[ty_to_be]+1;
            Particle particle;
            bool labeled=false;
            int id=getParticleIndex();
            int a[6]={id,labeled,temp_2[0],temp_2[1],temp_2[2],ty_to_be};
            assignParticle(particle,a);
            if (C2->occupiedParticle.empty())  C2->occupiedParticle.push_back(particle);
            else{
                cout<<"error!  step "<<step<<" "<<C2->occupiedParticle[0].index<<endl;
                cout<<"already have particles! Wrong";
                system("pause");
            }
            if (id<particleNum)
                listOfParticles[id]=particle;
            else listOfParticles.insert(listOfParticles.end(),particle);
            particleNum++;

        }
        else if(temp_2[2]==x_coor2-1 || temp_2[2]==getTrack(temp_2[0])->length+x_coor2){//exit
            Unit* C1=getUnit(temp_1[0],temp_1[1],temp_1[2]);//exit
            C1->occupiedNum[ty_to_move]=C1->occupiedNum[ty_to_move]-1;
            int i=0;
            while(i<static_cast<int>(C1->occupiedParticle.size()) && C1->occupiedParticle.at(i).type!=ty_to_move){
                i++;
            }
            particleNum--;

            Particle & particle=C1->occupiedParticle[i];
            if (flagRunLengthEnd || flagRunLengthBundle){  //REMOVE exit particles
               if (particle.labeled==true){
                   removeRunLength(particle);
               }
            }
            if(flagPause && particle.labeled) removePause(particle);
            C1->occupiedParticle.erase(C1->occupiedParticle.begin()+i);
            particleIndexes.insert(particleIndexes.begin(),particle.index);
            listOfParticles[particle.index].location.clear();
        }
        else{
            Unit* C1=getUnit(temp_1[0],temp_1[1],temp_1[2]);
            Unit* C2=getUnit(temp_2[0],temp_2[1],temp_2[2]);

            if (C1->occupiedNum[ty_to_move]<=0)     cout<<"Wrong occupied number is negative !"<<endl;

            //update units
            C1->occupiedNum[ty_to_move]=C1->occupiedNum[ty_to_move]-1;
            C2->occupiedNum[ty_to_be]=C2->occupiedNum[ty_to_be]+1;
            //which particle in the list to be updated
            int i=0;
            while(i<static_cast<int>(C1->occupiedParticle.size()) && C1->occupiedParticle.at(i).type!=ty_to_move){
                i++;
            }
            Particle & particle=C1->occupiedParticle[i];
            particle.moveTo(temp_2,ty_to_be);// care of reference here
            listOfParticles[particle.index].location=temp_2;
            listOfParticles[particle.index].type=ty_to_be;

            if (flagRunLengthEnd|| flagRunLengthBundle)
                updateRunLength(temp_1,temp_2,particle, ty_to_move);
            if(flagPause && particle.labeled) recordPause(particle);

            C2->occupiedParticle.insert(C2->occupiedParticle.end(),particle);
            C1->occupiedParticle.erase(C1->occupiedParticle.begin()+i);

            if (C1->occupiedNum[ty_to_move]<0 || C2->occupiedNum[ty_to_be]<0 || C1->occupiedNum[ty_to_move]>maxAllowance(ty_to_be) || C2->occupiedNum[ty_to_move]>maxAllowance(ty_to_be)){
                cout<<"error: update type "<<C1->occupiedNum[ty_to_move]<<" to "<<C2->occupiedNum[ty_to_be];
                showSingleRate(n);
                system("pause");
            }
        }
        return n;
    }
    else{
        cout<<"system is in stationary state as no events could occur; remind that only rate values greater than "<<rMin<<" is recorded!";
        t=tmax+1;
        showAllTypes();
    }
}



void Bundle::initializeBundle(){
        bundleState.clear();
        vector<int> trackwidth(trackNum,laneNum);// give lane number of each track;
        vector<int> tracklength(trackNum,siteNum);// give site number of each track;
        for (int i=0;i<trackNum;i++){
            Track temp(tracklength[i],trackwidth[i]);
            bundleState.insert(bundleState.end(),temp);
        }
        trackwidth.clear();
        tracklength.clear();
}

void Bundle::initializeBundleVar(){
    flagRunLengthEnd=false;//flagrunlengthend;
    flagRunLengthBundle=false;//flagrunlengthbundle;
    flagDen=false;//flagden;
    flagCurrent=false;//flagcurrent;
    currentCounter=0;
    flagInhomo=flaginhomo;
    flagTipSize=false;//flagtipsize;
    unaccessibility=flagunaccessible;//set unaccessible region as blockages
    flagColocalizeTip=false;//flagcolocalizetip;
    flagPause=false;//flagpause;
    flagDirPro=false;
    flagEEDis=false;
    flagEEDistriScore=false;
}

void Bundle::initializeRateVector(){
    rate.clear();
    for (int track=0;track<trackNum;track++){
        int x_coor=bundleState[track].xCoordinate;
        for (int site=0;site<static_cast<int>(bundleState[track].istate[0].size());site++){
            addRate(site+x_coor,track);
        }
    }
}


void Bundle::initializeTransitionRates(){
        pLoad=p_load;
        velocity.clear();
        runLengthBundled.clear();
        runLengthUnbundled.clear();
        transitionRates.forward.clear();
        transitionRates.forwardEnd.clear();
        transitionRates.laneChange.clear();
        transitionRates.trackSwitchAtEnd.clear();
        transitionRates.trackSwitchAtInner.clear();
        transitionRates.turning.clear();
        transitionRates.boundaryIn.clear();
        transitionRates.boundaryOut.clear();

        velocity.insert(velocity.end(),v_k);
        velocity.insert(velocity.end(),v_d);
        runLengthUnbundled.insert(runLengthUnbundled.end(),M_k);
        runLengthUnbundled.insert(runLengthUnbundled.end(),M_d);
        runLengthBundled.insert(runLengthBundled.end(),M_k);
        runLengthBundled.insert(runLengthBundled.end(),M_d);
        double p_d=v_d/stepSize-2*lane_d;
        double p_k=v_k/stepSize;

        transitionRates.forward.insert(transitionRates.forward.end(),p_k);
        if (p_d>0)
            transitionRates.forward.insert(transitionRates.forward.end(),p_d);
        else
            transitionRates.forward.insert(transitionRates.forward.end(),0);
        transitionRates.forward.insert(transitionRates.forward.end(),0);// crosing minus end on the same track

        transitionRates.turning.insert(transitionRates.turning.end(),w_k);//notbundled
        transitionRates.turning.insert(transitionRates.turning.end(),w_k);//isbundled
        transitionRates.turning.insert(transitionRates.turning.end(),w_d);
        transitionRates.turning.insert(transitionRates.turning.end(),w_d);//isbundled
        transitionRates.turning.insert(transitionRates.turning.end(),p_mturn1);//inhomogeneous at minus ends of its own track
        transitionRates.turning.insert(transitionRates.turning.end(),p_mturn2);//inhomogeneous at minus ends of another track

        //27_11_14 add inhomogeneous turning when meeting plus ends of another MT for dynein
        transitionRates.turning.insert(transitionRates.turning.end(),p_pturn2);


        transitionRates.laneChange.insert(transitionRates.laneChange.end(),lane_k);
        transitionRates.laneChange.insert(transitionRates.laneChange.end(),p_k/2);//blocked;
        if (p_d>0)
            transitionRates.laneChange.insert(transitionRates.laneChange.end(),lane_d);
        else
            transitionRates.laneChange.insert(transitionRates.laneChange.end(),v_d/stepSize/2);
        transitionRates.laneChange.insert(transitionRates.laneChange.end(),v_d/stepSize/2);//blocked;

        transitionRates.trackSwitchAtEnd.insert(transitionRates.trackSwitchAtEnd.end(),0);//kinesin unibundle into
        transitionRates.trackSwitchAtEnd.insert(transitionRates.trackSwitchAtEnd.end(),0); //E2  kinesin unibundle bundle leave & type change
        transitionRates.trackSwitchAtEnd.insert(transitionRates.trackSwitchAtEnd.end(),sek);//E1 kinesin unibundle bundle leave & no type chane
        transitionRates.trackSwitchAtEnd.insert(transitionRates.trackSwitchAtEnd.end(),0);// kinesin bipolar bundle into
        transitionRates.trackSwitchAtEnd.insert(transitionRates.trackSwitchAtEnd.end(),0);//E4 kinesin bipolar bundle leave & type change
        transitionRates.trackSwitchAtEnd.insert(transitionRates.trackSwitchAtEnd.end(),sek);//E3 kinesin bipolar bundle leave & no type chane
        transitionRates.trackSwitchAtEnd.insert(transitionRates.trackSwitchAtEnd.end(),0);//dynein to another track with the same directon
        transitionRates.trackSwitchAtEnd.insert(transitionRates.trackSwitchAtEnd.end(),sed);//E6 D-K
        transitionRates.trackSwitchAtEnd.insert(transitionRates.trackSwitchAtEnd.end(),0);//E7 D-D
        transitionRates.trackSwitchAtEnd.insert(transitionRates.trackSwitchAtEnd.end(),se);//E5 se bipolar D-K
        transitionRates.trackSwitchAtEnd.insert(transitionRates.trackSwitchAtEnd.end(),sed);//E8 q_2*v_d/hs);// D-K
        transitionRates.trackSwitchAtEnd.insert(transitionRates.trackSwitchAtEnd.end(),0);//E9 D-D

        transitionRates.trackSwitchAtInner.insert(transitionRates.trackSwitchAtInner.end(),0);//K-D//
        transitionRates.trackSwitchAtInner.insert(transitionRates.trackSwitchAtInner.end(),si);//K-K;//
        transitionRates.trackSwitchAtInner.insert(transitionRates.trackSwitchAtInner.end(),si);//D-K;/
        transitionRates.trackSwitchAtInner.insert(transitionRates.trackSwitchAtInner.end(),0);//D-D;/
       // transitionRates.forwardEnd.insert(transitionRates.forwardEnd.end(),0);//dynein cross minus end on the same track
        transitionRates.forwardEnd.insert(transitionRates.forwardEnd.end(),1);//kinesin cross a minus end of another track  unipolar bundle
        transitionRates.forwardEnd.insert(transitionRates.forwardEnd.end(),1);//kinesin cross a plus end of another track  bipolar
        transitionRates.forwardEnd.insert(transitionRates.forwardEnd.end(),1);//dynein cross a plus end of another track unipolar bundle
        transitionRates.forwardEnd.insert(transitionRates.forwardEnd.end(),q);//dynein cross a minus end of another track bipolar
        // forwardEnd rate will be multiplied by v/hs; switch rate will be multiplied by 1/lane numbers
        vector<double> alpha_k(trackNum,0);
        vector<double> alpha_d(trackNum,0);
        vector<double> beta_k(trackNum,0);
        vector<double> beta_d(trackNum,0);

        transitionRates.boundaryIn.insert(transitionRates.boundaryIn.end(),alpha_k);
        transitionRates.boundaryIn.insert(transitionRates.boundaryIn.end(),alpha_d);
        transitionRates.boundaryOut.insert(transitionRates.boundaryOut.end(),beta_k);
        transitionRates.boundaryOut.insert(transitionRates.boundaryOut.end(),beta_d);
}

bool Bundle::isBundled(const int site){
    vector<int> temp;
    for (int i=0;i<trackNum;i++){
         temp.insert(temp.end(),getTrack(i)->minusEndLocations.begin(),getTrack(i)->minusEndLocations.end());
         temp.insert(temp.end(),getTrack(i)->plusEndLocations.begin(),getTrack(i)->plusEndLocations.end());
    }
    sort(temp.begin(),temp.end());
    if (site>=temp[1] && site<temp[temp.size()-2]){
        return true;
    }
    else return false;
}



void Bundle::kymographRecord(const int time_size, bool current){
    vector<int> temp;
    
   if(current){
      for (int i=0;i<maxLength;i++){
        int k=0;
        int d=0;
        for (int tr=0;tr<trackNum;tr++){
            int xco=getXCoordinate(tr);
            int la=0;
            if (i>=xco && i<xco+bundleState[tr].length){
                while(k==0 && la<static_cast<int>(bundleState[tr].istate.size())){
                       //   cout<<bundleState[tr].length-bundleState[tr].istate[0].size()<<" "<<tr<<" "<<" "<<la<<" "<<i-xco<<"    ";
                    k=bundleState[tr].istate[la][i-xco].occupiedNum[0];
                    if (k>0){
                        for (int ty=0;ty<static_cast<int>(getTrack(tr)->istate[la][i-xco].occupiedParticle.size());ty++)
                            k=k+10*getTrack(tr)->istate[la][i-xco].occupiedParticle[ty].labeled;
                    }
                    la++;
                }              
               la=0;
                while(d==0 && la<static_cast<int>(bundleState[tr].istate.size())){
                   d=2*bundleState[tr].istate[la][i-xco].occupiedNum[1];
                   if (d>0){
                       for (int ty=0;ty<static_cast<int>(getTrack(tr)->istate[la][i-xco].occupiedParticle.size());ty++)
                           d=d+2*10*getTrack(tr)->istate[la][i-xco].occupiedParticle[ty].labeled;
                   }
                   la++;
               }

            }
        }
        temp.insert(temp.end(),k+d);
      }
      vkymograph.insert(vkymograph.begin(),temp);
      temp.clear();
    }
    else{
        if (vkymograph.size()>0){
            vkymograph.insert(vkymograph.begin(),vkymograph.front());
        }
    }

    if (static_cast<int>(vkymograph.size())>time_size) vkymograph.erase(vkymograph.begin()+time_size,vkymograph.end());
    else if (static_cast<int>(vkymograph.size())<time_size){
        vector<int> temp2(maxLength,0);
        while(static_cast<int>(vkymograph.size())<time_size){
            vkymograph.insert(vkymograph.end(),temp2);
        }
        temp2.clear();
    }
}

void Bundle::labelTip(const int type, const double per){
    int xstart=std::min(int(double(maxLength)*per),maxLength-1);
    int xend=maxLength;
    labelParticles(type,xstart,xend);
}


void Bundle::labelBundle(){
    /*vector<int> temp;
    getBundleRegion(temp);
    int xstart=int(double(temp[1]-temp[0])*labelPer);
    int xend=temp[1];
    labelParticles(1,xstart,xend);
*/
}


void Bundle::labelParticles(const int type, const int xstart, const int xend){
    bool checked=true;
    for (int track=0;track<trackNum;track++){
        int length=getXCoordinate(track)+getTrack(track)->length;
        if (xstart>=getXCoordinate(track) && xstart<length){
            for (int lane=0;lane<getTrack(track)->width;lane++){
                for (int site=xstart;site<min(xend,length);site++){
                    for (int i=0;i<static_cast<int>(getUnit(track,lane,site)->occupiedParticle.size());i++){
                        Particle & particle=getUnit(track,lane,site)->occupiedParticle[i];
                        if (particle.type==type){
                            setParticleLabel(particle,checked);
                            listOfParticles[particle.index].labeled=checked;
                            if (flagRunLengthEnd || flagRunLengthBundle){
                                double a[]={double(particle.index),double(particle.location[2]),t};
                                vector<double> temp(a,a+3);
                                runLength.insert(runLength.end(),temp);
                            }
                        }
                    }
                }

            }
        }
    }
}



/*void Bundle::minRate(){
    vector<double>rtemp;
    rtemp.insert(rtemp.begin(),transitionRates.forward.begin(),transitionRates.forward.end());
    rtemp.insert(rtemp.begin(),transitionRates.forwardEnd.begin(),transitionRates.forwardEnd.end());
    rtemp.insert(rtemp.begin(),transitionRates.laneChange.begin(),transitionRates.laneChange.end());
    rtemp.insert(rtemp.begin(),transitionRates.trackSwitchAtEnd.begin(),transitionRates.trackSwitchAtEnd.end());
    rtemp.insert(rtemp.begin(),transitionRates.trackSwitchAtInner.begin(),transitionRates.trackSwitchAtInner.end());
    rtemp.insert(rtemp.begin(),transitionRates.turning.begin(),transitionRates.turning.end());

    for (int ty=0;ty<typeNum;ty++){
       rtemp.insert(rtemp.begin(),transitionRates.boundaryIn[ty].begin(),transitionRates.boundaryIn[ty].end());
       rtemp.insert(rtemp.begin(),transitionRates.boundaryOut[ty].begin(),transitionRates.boundaryOut[ty].end());
    }
     for (std::vector<double>::size_type i=0;i<rtemp.size();){
        if (rtemp[i]==0){
            rtemp.erase(rtemp.begin()+i);
        }
        else i++;
    }
    rMin=*(min_element(rtemp.begin(),rtemp.end()));
    if (rMin<0) system("pause");
}
*/


void Bundle::oneendInitializeState(){
    int N_temp=0;
    particleIndexes.clear();
    listOfParticles.clear();
    runLength.clear();
    colocalization.clear();
    pauses.clear();
    for (int track=0;track<static_cast<int>(bundleState.size());track++){
        for (int i=static_cast<int>(bundleState[track].istate[0].size())-1;i>=0;--i){
            for (int j=0;j<static_cast<int>(bundleState[track].istate.size());j++){
                if (N_temp<particleNum && bundleState[track].istate[j][i].unblocked){
                    bundleState[track].istate[j][i].occupiedNum[1]=1;
                    bundleState[track].istate[j][i].occupiedNum[0]=0;
                    Particle particle;
                    bool labeled=false;
                    int a[6]={N_temp,labeled,track,j,i+getXCoordinate(track),1};
                    assignParticle(particle,a);
                    listOfParticles.insert(listOfParticles.end(),particle);
                    bundleState[track].istate[j][i].occupiedParticle.clear();
                    bundleState[track].istate[j][i].occupiedParticle.push_back(particle);
                    N_temp++;
                }
                else {
                    bundleState[track].istate[j][i].occupiedNum[1]=0;
                    bundleState[track].istate[j][i].occupiedNum[0]=0;
                    getTrack(track)->istate[j][i].occupiedParticle.clear();
                }
            }
        }
    }
}


void Bundle::tipendInitializeState(){
    int N_temp=0;
    particleIndexes.clear();
    listOfParticles.clear();
    runLength.clear();
    colocalization.clear();
    pauses.clear();

    int tipTrack=getTipTrack();

    for (int track=0;track<static_cast<int>(bundleState.size());track++){
        //track=bundleState.size()-track;
        //for (int i=50;i>=0;--i){ // for starting in the middle for specific movies
        for (int i=static_cast<int>(bundleState[track].istate[0].size())-1;i>=0;--i){
          for (int j=0;j<static_cast<int>(bundleState[track].istate.size());j++){
            int xco=getTrack(track)->xCoordinate;
           // if ((track==3 || track==2 || track==4) && (i+xco<=50 && ((i+xco>48) || i+xco>47 & track==4)) && N_temp<particleNum && bundleState[track].istate[j][i].unblocked){// for starting in the middle for specific movies
            if (track==tipTrack && N_temp<particleNum && bundleState[track].istate[j][i].unblocked){
                int xi=maxLength;//static_cast<int>(bundleState[track].istate[0].size())-1;
                int ty;
                if (getTrack(track)->getOrientation(xi)==0){
                bundleState[track].istate[j][i].occupiedNum[1]=0;//1;
                bundleState[track].istate[j][i].occupiedNum[0]=1;//0;
                ty=0;//1;
                }
                else{
                    bundleState[track].istate[j][i].occupiedNum[1]=0;
                    bundleState[track].istate[j][i].occupiedNum[0]=1;
                    ty=0;
                }
                Particle particle;
                bool labeled=false;
                int a[6]={N_temp,labeled,track,j,i+getXCoordinate(track),ty};
                assignParticle(particle,a);
                listOfParticles.insert(listOfParticles.end(),particle);
                bundleState[track].istate[j][i].occupiedParticle.clear();
                bundleState[track].istate[j][i].occupiedParticle.push_back(particle);
                N_temp++;
            }
            else {
                bundleState[track].istate[j][i].occupiedNum[1]=0;
                bundleState[track].istate[j][i].occupiedNum[0]=0;
                getTrack(track)->istate[j][i].occupiedParticle.clear();
            }
          }
        }
    }

}

void Bundle::randomInitializeState(long *idum){
    int N_temp=0;
    runLength.clear();
    particleIndexes.clear();
    listOfParticles.clear();
    colocalization.clear();
    pauses.clear();
    for (int track=0;track<static_cast<int>(bundleState.size());track++){
            for (int lane=0;lane<static_cast<int>(bundleState[track].istate.size());lane++){
                for (int site=0;site<static_cast<int>(bundleState[track].istate[lane].size());site++){
                    for (int ty=0;ty<static_cast<int>(bundleState[track].istate[lane][site].occupiedNum.size());ty++){
                         bundleState[track].istate[lane][site].occupiedNum[ty]=0;
                    }
                    getTrack(track)->istate[lane][site].occupiedParticle.clear();
                }
            }
     }
    while(N_temp<particleNum){
        int track,lane,site,ty;
        double r=ran2(idum);
        track=int(double(bundleState.size())*r);
        r=ran2(idum);
        lane=int(double(bundleState[track].istate.size())*r);
        r=ran2(idum);
        site=int(double(bundleState[track].istate[lane].size())*r);
        if (bundleState[track].istate[lane][site].getTotalOccupiedNum()==0){
            r=ran2(idum);
            ty=int(double(bundleState[track].istate[lane][site].occupiedNum.size())*r);
            if (bundleState[track].istate[lane][site].occupiedNum[ty]==0 && bundleState[track].istate[lane][site].unblocked){
                bundleState[track].istate[lane][site].occupiedNum[ty]=1;
                Particle particle;
                bool labeled=false;
                int a[6]={N_temp,labeled,track,lane,site+getXCoordinate(track),ty};
                assignParticle(particle,a);
                listOfParticles.insert(listOfParticles.end(),particle);
                bundleState[track].istate[lane][site].occupiedParticle.clear();
                bundleState[track].istate[lane][site].occupiedParticle.push_back(particle);
                N_temp=N_temp+1;
            }
        }
    }

}

void Bundle::removePause(const Particle & particle){
    if (pauses.size()>0){
        int i=0;
        while(i<static_cast<int>(pauses.size()) && !(int(pauses[i][0])==particle.index && pauses[i].size()==2)){
            i++;
        }
        pauses[i].clear();
    }
}

void Bundle::recordPause(const Particle & particle){
    if (pauses.size()>0){
        int i=0;
        while(i<static_cast<int>(pauses.size()) && !(int(pauses[i][0])==particle.index && pauses[i].size()==2)){
            i++;
        }
        if (i<static_cast<int>(pauses.size()) && (t-pauses[i][1]>=pausemin)){
            pauses[i].insert(pauses[i].end(),t);
        }
        else if(i<static_cast<int>(pauses.size()))
                pauses[i][1]=t;
        else{
            double a[2]={double(particle.index),t};
            std::vector<double>temp(a,a+2);
            pauses.insert(pauses.end(),temp);
        }
    }
    else{
        double a[2]={double(particle.index),t};
        std::vector<double>temp(a,a+2);
        pauses.insert(pauses.end(),temp);
    }
}

void Bundle::recordTipSize(){
    extern double t;
    int tip=calculateTipSize();
    double a[2]={t,double(tip)};
    vector<double> vtemp(a,a+2);
    tipSize.insert(tipSize.end(),vtemp);
}

void Bundle::removeUnaccessibleRegion(){
    for (int i=0;i<trackNum;i++){
        int xco=getTrack(i)->xCoordinate;
        int len=xco+getTrack(i)->length;
        for (int k=0;k<static_cast<int>(getTrack(i)->unaccessibleRegion.size());k++){
            std::vector<int>temp=getTrack(i)->unaccessibleRegion[k];
            for (int y=temp[0];y<=temp[1];y++){
                for (int x=xco;x<len;x++){
                //for (int x=max(temp[2]-2,xco);x<=min(temp[3]+2,len);x++){
                    getTrack(i)->istate[y][x-xco].unblocked=true;
                }
            }
        }
        getTrack(i)->unaccessibleRegion.clear();
    }
}

void Bundle::setParticleLabel(Particle & particle, bool & labeled){
    if (labeled) {
        particle.labeled=true;
    }
    else {
        particle.labeled=false;
    }
}

void Bundle::setRandomBlockages(const int n, long *idum){
    int N_temp=0;
    for (int track=0;track<static_cast<int>(bundleState.size());track++){
            for (int lane=0;lane<static_cast<int>(bundleState[track].istate.size());lane++){
                for (int site=0;site<static_cast<int>(bundleState[track].istate[lane].size());site++){
                      bundleState[track].istate[lane][site].unblocked=true;
                }
            }
    }
    if (unaccessibility){
        setUnaccessibleRegion();
        setUnaccessibleBlockage();
    }
    while(N_temp<n){
        int track,lane,site;
        double r=ran2(idum);
        track=int(double(bundleState.size())*r);
        //int xco=getXCoordinate(track);
        r=ran2(idum);
        lane=int(double(bundleState[track].istate.size())*r);
        r=ran2(idum);
        site=int(double(bundleState[track].istate[lane].size())*r);
        if (bundleState[track].istate[lane][site].unblocked==true){// && site+xco<9*maxLength/10){
            bundleState[track].istate[lane][site].unblocked=false;
            N_temp=N_temp+1;
        }
    }
}

void Bundle::setUnaccessibleRegion(){
    for (int i=0;i<trackNum;i++) getTrack(i)->unaccessibleRegion.clear();

    for (int i=0;i<trackNum-1;i++){
        int xco1=getTrack(i)->xCoordinate;
        int len1=getTrack(i)->length;
        int ml1=int(getTrack(i)->width/2);//middle lane
        int bw1=int(ceil(double(getTrack(i)->width*3)/double(13)));
        int l0=ml1-bw1;
        int l1=l0+bw1-1;
        int l2=l1+bw1-1;
        int l3=l2+bw1-1;
        for (int j=i+1;j<trackNum;j++){
            int xco2=getTrack(j)->xCoordinate;
            int len2=getTrack(j)->length;
            int ml2=int(getTrack(j)->width/2);//middle lane
            int bw2=int(ceil(double(getTrack(j)->width*3)/double(13)));
            int l4=ml2-bw2;
            int l5=l4+bw2-1;
            int l6=l5+bw2-1;
            int l7=l6+bw2-1;
            int xbegin=max(xco1,xco2);
            int xend=min(xco1+len1,xco2+len2)-1;
            // check intersection with previous set regions
            if (xend>=xbegin){
          //  std::vector<std::vector<int> > intersect1,intersect2;
            std::vector<int> vlane1,vlane2;
            for (int k=0;k<static_cast<int>(getTrack(i)->unaccessibleRegion.size());k++){
                std::vector<std::vector<int> > temp=getTrack(i)->unaccessibleRegion;
                int x1=max(xbegin,temp[k][2]);
                int x2=min(xend,temp[k][3]);
                if (x2>=x1){
                   // intersect1.insert(intersect1.end(),temp[k]);
                    vlane1.insert(vlane1.end(),temp[k][0]);
                }
            }
            if(find(vlane1.begin(),vlane1.end(),l1)==vlane1.end()){
                int reg[4]={l1,l2,xbegin,xend};
                std::vector<int> vreg(reg,reg+4);
                getTrack(i)->unaccessibleRegion.insert(getTrack(i)->unaccessibleRegion.end(),vreg);
            }
            else if(find(vlane1.begin(),vlane1.end(),l2)==vlane1.end()){
                int reg[4]={l2,l3,xbegin,xend};
                std::vector<int> vreg(reg,reg+4);
                getTrack(i)->unaccessibleRegion.insert(getTrack(i)->unaccessibleRegion.end(),vreg);
            }
            else if(find(vlane1.begin(),vlane1.end(),l0)==vlane1.end()){
                int reg[4]={l0,l1,xbegin,xend};
                std::vector<int> vreg(reg,reg+4);
                getTrack(i)->unaccessibleRegion.insert(getTrack(i)->unaccessibleRegion.end(),vreg);
            }

            for (int k=0;k<static_cast<int>(getTrack(j)->unaccessibleRegion.size());k++){
                std::vector<std::vector<int> > temp=getTrack(j)->unaccessibleRegion;
                int x1=max(xbegin,temp[k][2]);
                int x2=min(xend,temp[k][3]);
                if (x2>=x1){
                  //  intersect2.insert(intersect2.end(),temp[k]);
                    vlane2.insert(vlane2.end(),temp[k][0]);
                }
            }
            if(find(vlane2.begin(),vlane2.end(),l5)==vlane2.end()){
                int reg[4]={l5,l6,xbegin,xend};
                std::vector<int> vreg(reg,reg+4);
                getTrack(j)->unaccessibleRegion.insert(getTrack(j)->unaccessibleRegion.end(),vreg);
            }
            else if(find(vlane2.begin(),vlane2.end(),l6)==vlane2.end()){
                int reg[4]={l6,l7,xbegin,xend};
                std::vector<int> vreg(reg,reg+4);
                getTrack(j)->unaccessibleRegion.insert(getTrack(j)->unaccessibleRegion.end(),vreg);
            }
            else if(find(vlane2.begin(),vlane2.end(),l7)==vlane2.end()){
                int reg[4]={l4,l5,xbegin,xend};
                std::vector<int> vreg(reg,reg+4);
                getTrack(j)->unaccessibleRegion.insert(getTrack(j)->unaccessibleRegion.end(),vreg);
            }
            }
        }
    }
  /*  if (trackNum>6){
    // manual corrections for eg_13 copy.txt
    getTrack(3)->unaccessibleRegion[3][0]=7;
    getTrack(3)->unaccessibleRegion[3][1]=9;
    getTrack(3)->unaccessibleRegion[4][0]=5;
    getTrack(3)->unaccessibleRegion[4][1]=7;
    getTrack(3)->unaccessibleRegion[5][0]=9;
    getTrack(3)->unaccessibleRegion[5][1]=11;
    }*/
}



void Bundle::setRandomBundleGeometry(std::vector<std::vector<int> > & randBundle, long *idum,int mdouble){
    randBundle.clear();
    std::vector<int> plus,minus;
    bool flag=false;
    std::vector<int> temp,addm;
    std::vector<int>::iterator it;
    while (!flag){
        generatePMends(minus,plus,idum,mdouble);
        int np=plus.size();
        int nm=minus.size();
        int xmin=plus[0];
        int xmax=plus[np-1];
        for (int i=0;i<static_cast<int>(plus.size());i++)
           plus[i]=int(double(plus[i]-xmin)/double(xmax-xmin)*100.0/stepSize);
        plus[np-1]=100.0/stepSize+1;
        for (int i=0;i<static_cast<int>(minus.size());i++)
           minus[i]=int(double(minus[i]-xmin)*100.0/stepSize/double(xmax-xmin));

        temp.clear();
        temp.insert(temp.end(),plus.begin(),plus.end());
        temp.insert(temp.end(),minus.begin(),minus.end());
        std::sort(temp.begin(),temp.end());
        it = std::unique (temp.begin(), temp.end());

        //int mdouble=n;//check with plus number in generatePMends()

        if (it==temp.end() && minus[0]>plus[0] && plus[np-1]>minus[nm-1] && plus[np-3]<minus[nm-1] && plus[2]>minus[0]){
            std::vector<int> add=minus;
            std::random_shuffle(add.begin(), add.end(), myrandom);
            addm.clear();
            addm.insert(addm.end(),add.begin(),add.begin()+mdouble);
            std::sort(addm.begin(),addm.end());

            if (plus[1]<addm[0] && plus[mdouble-2]<addm[mdouble-1] && plus[np-(mdouble-1)]>addm[1] && plus[np-mdouble]>addm[0]){
                minus.insert(minus.end(),addm.begin(),addm.begin()+mdouble);
                std::sort(minus.begin(),minus.end());
                flag=true;
            }
        }

    }
    std::vector<std::vector<int> > MINUS;

    std::vector<int> perm(minus.begin(),minus.end());//(myints,myints+10);
    std::sort(perm.begin(),perm.end());
    int counts=0;

    std::vector<double> center,distri;
    for (int i=0;i<static_cast<int>(temp.size())-1;i++){
        center.insert(center.end(),double(temp[i+1]+temp[i])/2.0);
    }

    //std::sort(minus.begin(),minus.end());
   // it = std::unique(minus.begin(), minus.end());//minus has been changed
    //int nm=perm.size()-addm.size();//int(std::distance(minus.begin(),it));
    do {
        distri.clear();
        distri.resize(center.size());
        std::vector<int> MTlen;
        for (int i=0;i<static_cast<int>(perm.size());i++){
            MTlen.insert(MTlen.end(),abs(perm[i]-plus[i]));
        }
        for (int k=0;k<static_cast<int>(perm.size());k++){
            for (int m=0;m<static_cast<int>(center.size());m++){
                if ((center[m]-plus[k])*(center[m]-perm[k])<=0)
                    distri[m]++;
            }
        }
        bool flagmdouble=true;
        for (int j=0;j<static_cast<int>(addm.size());j++){
            int c=0;
            for (int k=0;k<static_cast<int>(perm.size());k++){
                if (((addm[j]-plus[k])*(addm[j]-perm[k]))<0){
                    c++;
                    break;
                }
            }
            if (c<1) {
                flagmdouble=false;
                break;
            }
        }


        bool flagbipolar=true;// to check bipolar for double minus end
        for (int i=0;i<static_cast<int>(addm.size());i++){
            it = find (perm.begin(),perm.end(), addm[i]);
            int n1=std::distance(perm.begin(),it);
            it = find (perm.begin()+n1+1,perm.end(), addm[i]);
            int n2=std::distance(perm.begin(),it);
            if ((plus[n1]-addm[i])*(plus[n2]-addm[i])>=0){
                flagbipolar=false;
                break;
            }

        }

      /* std::vector<std::vector<int> > MTs;
        MTs.insert(MTs.end(),plus);
        MTs.insert(MTs.end(),perm);
        int MT[2];

        for (int i=0;i<nm;i++){
            MT[1]=MTs[1][0];//minus
            MT[0]=MTs[0][0];//plus
            MTs[0].erase(MTs[0].begin());
            MTs[1].erase(MTs[1].begin());
            it = find (MTs[1].begin(), MTs[1].end(), MT[1]);
            if (it!=MTs[1].end()){
                int ind=it-MTs[1].begin();
                if ((MTs[0][ind]-MT[1])*(MT[0]-MT[1])>=0){
                    flagbipolar=false;
                    break;
                }
                MTs[0].erase(MTs[0].begin()+ind);
                MTs[1].erase(MTs[1].begin()+ind);
            }
        }*/
        counts++;
        //condition on MT length, distribution of #MTs, & double minus ends
        if (flagmdouble && flagbipolar && *min_element(distri.begin(),distri.end())>=1 && *max_element(distri.begin(),distri.end())<=3 && *min_element(MTlen.begin(),MTlen.end())>=int(5.0/stepSize) && *max_element(MTlen.begin(),MTlen.end())<=int(80.0/stepSize)){
            MINUS.insert(MINUS.end(),perm);
        }
    } while ( std::next_permutation(perm.begin(),perm.end()) );


    if (MINUS.size()>0){
       randBundle.insert(randBundle.end(),plus);
       int kk=int(ran2(idum)*MINUS.size());
       randBundle.insert(randBundle.end(),MINUS[kk]);
    }
   // counts=counts+1;


}

void Bundle::generatePMends(std::vector<int> & minus,std::vector<int> & plus,long * idum,int n){
    plus.clear();
    minus.clear();

    double x,r2,p;

    int i=1;
    int x0;
    double cg=0.2616;
    int np=7+n;
    while (i<=np){
        x=ran2(idum)*100;
        r2=ran2(idum);
        if ((x<10) & (x>=0))    p=-0.0176*(x-10)+(5-50*(0.0112+0.0176))/100;
        else if ((x>=10) & (x<=90))    p=0.0356;
        else if ((x>90) & (x<=100))  p=0.0112*(x-90)+(5-50*(0.0112+0.0176))/100;
        x0=100-int(floor(x));
        if (r2<=p/cg){
            // plus.insert(plus.end(),x0);
             std::vector<int>::iterator it=find(plus.begin(),plus.end(),x0);
             if (it==plus.end())
                 plus.insert(plus.end(),x0);
             i=plus.size();
             i=i+1;
        }
    }
    std::sort(plus.begin(),plus.end());

    double b[]={22.28,-1.583,0.5275,-0.02635,0.0005352,-4.847e-6,1.618e-8};//best fit for six order polynomial
    int nm=7;
    double cgm=0.0715;
    i=1;
    while (i<=nm){
        //x=10+ran2(idum)*80;//ran2(idum)*100;
        x=ran2(idum)*100;
        r2=ran2(idum);
        p=b[0]+b[1]*x+b[2]*pow(x,2)+b[3]*pow(x,3.0)+b[4]*pow(x,4.0)+b[5]*pow(x,5.0)+b[6]*pow(x,6.0);
        p=p/1000;
        if (r2<=p/cgm){
            x0=100-int(floor(x));
            if (find(minus.begin(), minus.end(),x0)==minus.end())  minus.insert(minus.end(),x0);
            i=minus.size();
            i=i+1;
        }
    }
    std::sort(minus.begin(),minus.end());
}

void Bundle::setUnaccessibleBlockage(){
    for (int i=0;i<trackNum;i++){
        std::vector<std::vector<int> > temp=getTrack(i)->unaccessibleRegion;
        int xco=getTrack(i)->xCoordinate;
        int len=getTrack(i)->length;
        for (int j=0;j<static_cast<int>(temp.size());j++){
            for (int y=temp[j][0];y<=temp[j][1];y++){
                for (int x=temp[j][2];x<=temp[j][3];x++)
                        getTrack(i)->istate[y][x-xco].unblocked=false;
            }
        }
        std::vector<std::vector<int> >intsec;
        std::vector<int> open;
        //left end
        for (int j=0;j<static_cast<int>(temp.size());j++){
            int k=0;
            bool flag=true;
            while(flag && k<static_cast<int>(temp.size())){
                if (k!=j){
                    if (temp[j][2]<=temp[k][3]+1 && temp[j][2]>=temp[k][2]) {
                        flag=false;
                        int a[2]={j,k};
                        std::vector<int> temp(a,a+2);
                        intsec.insert(intsec.end(),temp);
                    }
                }
                k++;
            }
            if (flag)    open.insert(open.end(),j);
        }
        for (std::vector<std::vector<int> >::size_type n=0;n<intsec.size();n++){
            int k=intsec[n][0];
            int j=intsec[n][1];
            int d=temp[k][1]-temp[k][0];
            if (temp[k][0]>temp[j][0]){
                for (int m=0;m<=d;m++){
                    int y1=temp[k][0];
                    int y2=temp[k][1]-1-m;
                    int x=max(temp[k][2]-m-1,xco);
                    for (int y=y1;y<=y2;y++){
                        getTrack(i)->istate[y][x-xco].unblocked=false;
                    }
                }
            }
            else if (temp[k][0]<temp[j][0]){
                for (int m=0;m<=d;m++){
                    int y1=temp[k][0]+1+m;
                    int y2=temp[k][1];
                    int x=max(temp[k][2]-m-1,xco);
                    for (int y=y1;y<=y2;y++){
                        getTrack(i)->istate[y][x-xco].unblocked=false;
                    }
                }
            }

        }
        for (std::vector<int>::size_type n=0;n<open.size();n++){
            int j=open[n];

            int l1=temp[j][0];
            int l2=temp[j][1];
            int x=max(temp[j][2],xco);
            int xx1=x;
            int xx2=x;
            while (xx1>=xco && !(getTrack(i)->istate[l1-1][xx1-xco].unblocked)){
                    xx1--;
            }
            while (xx2>=xco && !(getTrack(i)->istate[l2+1][xx2-xco].unblocked)){
                    xx2--;
            }
            if (min(xx1,xx2)>=xco){
                for (int m=x-1;m>=min(xx1,xx2);m--){
                    for (int y=l1;y<=l2;y++)
                        getTrack(i)->istate[y][m-xco].unblocked=false;
                }
            }
            int d=(l2-l1);
                for (int m=0;m<d/2;m++){
                    int x1=min(xx1,xx2)-m-1;
                    if (x1>=xco){
                        int y1=temp[j][0]+m+1;
                        int y2=temp[j][1]-m-1;
                        for (int y=y1;y<=y2;y++)
                           getTrack(i)->istate[y][x1-xco].unblocked=false;
                    }
                }
            }

        // right end
        open.clear();
        intsec.clear();
        for (int j=0;j<static_cast<int>(temp.size());j++){
            int k=0;
            bool flag=true;
            while(flag && k<static_cast<int>(temp.size())){
                if (k!=j){
                    if (temp[j][3]<=temp[k][3]+1 && temp[j][3]>=temp[k][2]) {
                        flag=false;
                        int a[2]={j,k};
                        std::vector<int> temp(a,a+2);
                        intsec.insert(intsec.end(),temp);
                    }
                }
                k++;
            }
            if (flag)    open.insert(open.end(),j);
        }
        for(std::vector<std::vector<int> >::size_type n=0;n<intsec.size();n++){
            int k=intsec[n][0];
            int j=intsec[n][1];
            int d=temp[k][1]-temp[k][0];
            if (temp[k][0]>temp[j][0]){
                for (int m=0;m<=d;m++){
                    int y1=temp[k][0];
                    int y2=temp[k][1]-1-m;
                    int x=min(temp[k][3]+m+1,len+xco-1);
                    for (int y=y1;y<=y2;y++){
                        getTrack(i)->istate[y][x-xco].unblocked=false;
                    }
                }
            }
            else if (temp[k][0]<temp[j][0]){
                for (int m=0;m<=d;m++){
                    int y1=temp[k][0]+1+m;
                    int y2=temp[k][1];
                    int x=min(temp[k][3]+m+1,xco+len-1);
                    for (int y=y1;y<=y2;y++){
                        getTrack(i)->istate[y][x-xco].unblocked=false;
                    }
                }
            }

        }
        for (std::vector<int>::size_type n=0;n<open.size();n++){
            int j=open[n];
            int l1=temp[j][0];
            int l2=temp[j][1];
            int x=min(temp[j][3],len+xco-1);
            int xx1=x;
            int xx2=x;
            while (xx1<xco+len && !(getTrack(i)->istate[l1-1][xx1-xco].unblocked)){
                    xx1++;
            }
            while (xx2<xco+len && !(getTrack(i)->istate[l2+1][xx2-xco].unblocked)){
                    xx2++;
            }
            if (max(xx1,xx2)<xco+len){
                for (int m=x+1;m<=max(xx1,xx2);m++){
                    for (int y=l1;y<=l2;y++)
                        getTrack(i)->istate[y][m-xco].unblocked=false;
                }
            }
            int d=(l2-l1);
            for (int m=0;m<d/2;m++){
                int x1=max(xx1,xx2)+1;
                if (x1<xco+len){
                    int y1=temp[j][0]+m+1;
                    int y2=temp[j][1]-m-1;
                    for (int y=y1;y<=y2;y++)
                       getTrack(i)->istate[y][x1-xco].unblocked=false;
                }
            }
            }

    }
}
void Bundle::sumCurrent(const int n){
    vector<int> temp_1=rate[n].loc1;
    vector<int> temp_2=rate[n].loc2;
    int x_coor1=bundleState[temp_1[0]].xCoordinate;
    int mlen=getTrack(temp_1[0])->length+x_coor1;
    int ty_to_move=rate[n].typeToMove;
   // int ty_to_be=rate[n].typeToBe;
    if(temp_1[2]>=x_coor1 &&  temp_1[2]<mlen){
            if (temp_2[2]!=temp_1[2]){
                int tr=0;
                int wid=0;
                while(tr<temp_1[0]){
                     wid=wid+getTrack(tr)->width;
                     tr++;
                }
                currentSum[ty_to_move][wid+temp_1[1]][temp_1[2]]++;
            }
    }
}

void Bundle::sumColocalize(){
    if (colocalization.size()==0){
              colocalization.resize(2);
              for (int i=0;i<2;i++)
              colocalization[i].resize(maxLength);
    }
    for (int i=0;i<static_cast<int>(listOfParticles.size());i++){
        std::vector<int> loc=listOfParticles[i].location;
        if (!loc.empty() && listOfParticles[i].labeled){
           // int w=getHeight(loc[0],loc[1]);
            colocalization[0][loc[2]]++;
            if (listOfParticles[i].type==1)
                colocalization[1][loc[2]]++;
        }
    }
}




void Bundle::sumState(const int n){
    if (denSum.size()==0){
          denSum.resize(typeNum);
          for (int ty=0;ty<typeNum;ty++){
              denSum[ty].resize(maxWidth);
              for (int j=0;j<maxWidth;j++)
                  denSum[ty][j].resize(maxLength);
          }
    }
    vector<int> temp_1=rate[n].loc1;
    vector<int> temp_2=rate[n].loc2;
    int x_coor1=bundleState[temp_1[0]].xCoordinate;
    int x_coor2=bundleState[temp_2[0]].xCoordinate;
    int ty_to_move=rate[n].typeToMove;
    int ty_to_be=rate[n].typeToBe;
    Unit* C1=getUnit(temp_1[0],temp_1[1],temp_1[2]);
    Unit* C2=getUnit(temp_2[0],temp_2[1],temp_2[2]);
    if (temp_1[2]==x_coor1-1 || temp_1[2]==getTrack(temp_1[0])->length+x_coor1){
            nstate[temp_2[0]][ty_to_be][temp_2[1]][temp_2[2]]+=1;
            if (nstate[temp_2[0]][ty_to_be][temp_2[1]][temp_2[2]]!=C2->occupiedNum[ty_to_be])
                cout<<"wrong injec"<<endl;
    }
    else if(temp_2[2]==x_coor2-1 || temp_2[2]==getTrack(temp_2[0])->length+x_coor2){//exit
            nstate[temp_1[0]][ty_to_move][temp_1[1]][temp_1[2]]-=1;
            if (nstate[temp_1[0]][ty_to_move][temp_1[1]][temp_1[2]]!=C1->occupiedNum[ty_to_move])
                cout<<"wrong exit"<<endl;
    }
    else{
            nstate[temp_1[0]][ty_to_move][temp_1[1]][temp_1[2]]-=1;
            nstate[temp_2[0]][ty_to_be][temp_2[1]][temp_2[2]]+=1;
            //error message
            if (nstate[temp_2[0]][ty_to_be][temp_2[1]][temp_2[2]]!=C2->occupiedNum[ty_to_be])
                cout<<"wrong in type to be"<<endl;
            if (nstate[temp_1[0]][ty_to_move][temp_1[1]][temp_1[2]]!=C1->occupiedNum[ty_to_move])
                cout<<"wrong in type to move"<<endl;
    }
        for (int ty=0;ty<typeNum;ty++){
            int y_co=0;
            for (int track=0;track<trackNum;track++){
                for (vector<vector<int> >::size_type y=0;y<bundleState[track].istate.size();y++){
                    transform(nstate[track][ty][y].begin(),nstate[track][ty][y].end(),denSum[ty][y_co].begin(),denSum[ty][y_co].begin(),op_densum);
                    y_co++;
                }
            }
        }
}


void Bundle::sumStateParticles(){
    for (std::vector<Particle>::size_type p=0;p<listOfParticles.size();p++){
        std::vector<int> loc=listOfParticles[p].location;
        int tr=loc[0];
        int wid=loc[1];
        int i=0;
        while(i<tr){
            wid=wid+getTrack(i)->width;
            i++;
        }
        int si=loc[2];
        int ty=listOfParticles[p].type;
        denSum[ty][wid][si]++;
    }
}

void Bundle::updateRate(const int ind, const int radius){
    bool flag_add=true;
    int site_1=rate[ind].loc1[2];
    int site_2=rate[ind].loc2[2];

    //delete rows in previous rate and add new ones
    for (vector<TransitionEvent>::size_type x=0;x<rate.size();){
        int track=rate[x].loc1[0];
        int x_coor=bundleState[track].xCoordinate;
        if (abs(rate[x].loc1[2]-site_1)<=radius || abs(rate[x].loc1[2]-site_2)<=radius){// care of boundary inject exit
            rate.erase(rate.begin()+x);
            flag_add=false;
        }
        else if ((site_1-radius ==x_coor|| site_2-radius==x_coor) && min(rate[x].loc1[2], rate[x].loc2[2])==x_coor-1 ){//left boundary
            rate.erase(rate.begin()+x);
            flag_add=false;
        }
        else {
            if (((site_1==static_cast<int>(bundleState[track].istate[0].size())+x_coor-1-radius) || (site_2==static_cast<int>(bundleState[track].istate[0].size())+x_coor-1-radius) ) && max(rate[x].loc1[2],rate[x].loc2[2])==static_cast<int>(bundleState[track].istate[0].size())+x_coor && rate[x].loc1[0]==track && rate[x].loc1[0]==rate[x].loc2[0]){
                rate.erase(rate.begin()+x);
                flag_add=false;
            }
            else flag_add=true;

        }
        if (flag_add)     x=x+1;
    }

    // add new rates
    for (int tr=0;tr<trackNum;tr++){
        int x_coor=bundleState[tr].xCoordinate;
        for (int i=max(min(site_1,site_2)-radius,x_coor);i<=min(max(site_2,site_1)+radius,static_cast<int>(bundleState[tr].istate[0].size())-1+x_coor);i++){
            addRate(i,tr);
        }
    }
}


void Bundle::updateRunLength(const std::vector<int> & temp_1, const std::vector<int> & temp_2, Particle & particle, int ty_to_move){
    int x_coor1=getTrack(temp_1[0])->xCoordinate;
    int mlen=getTrack(temp_1[0])->length+x_coor1;
    if(temp_1[2]>=x_coor1 &&  temp_1[2]<mlen){
        if (flagRunLengthEnd){
            //if (temp_1[2]==temp_2[2] && particle.labeled==true)// original
            //if (ty_to_be!=ty_to_move && temp_1[2]==temp_2[2] && particle.labeled==true){ modify 04_02_14 CP
            if (particle.type!=ty_to_move && particle.labeled==true) // modified to look at run length by dynein only  07/02/2018
             addRunLengthVale(particle);
            //only relabel dyneins at the tip
            else if (temp_2[2]>=labelRegion[0] && temp_2[2]<labelRegion[1] && particle.labeled==false){
                int track=getTipTrack();
                int xi=maxLength;//bundleState[track].istate[0].size()-1;
                if ( (particle.type==1 && getTrack(track)->getOrientation(xi)==0) ||  (particle.type==0 && getTrack(track)->getOrientation(xi)==1) )
                addRunLength(particle);
            }
        }
        else if (flagRunLengthBundle){
            std::vector<int>region;
            getBundleRegion(region);
            if (temp_1[2]==temp_2[2] && particle.labeled==true)
            //if (ty_to_be!=ty_to_move && temp_1[2]==temp_2[2] && particle.labeled==true){ modify 04_02_14 CP
              addRunLengthVale(particle);
            //only relabel dyneins at the tip
            else if (temp_2[2]>=labelRegion[0] && temp_2[2]<labelRegion[1] &&  particle.labeled==false){
                int track=getTipTrack();
                int xi=maxLength;//bundleState[track].istate[0].size()-1;
                if ( (particle.type==1 && getTrack(track)->getOrientation(xi)==0) ||  (particle.type==0 && getTrack(track)->getOrientation(xi)==1) )
                addRunLength(particle);
            }
            // remove particle label when out of the bundle
            else if ((temp_2[2]<region[0] || temp_2[2]>=region[1]) && particle.labeled){
                removeRunLength(particle);
                particle.labeled=false;
                listOfParticles[particle.index].labeled=false;
            }
        }
    }
}




void Bundle::removeRunLength(Particle &particle){
    int j=0;
    while(j<static_cast<int>(runLength.size()) && runLength[j][0]!=particle.index){
       j++;
    }
    if (j<static_cast<int>(runLength.size()))   runLength.erase(runLength.begin()+j);
    else{
        cout<<" Wrong, no element to delete!";
        system("pause");
    }
}

