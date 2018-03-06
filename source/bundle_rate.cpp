 /*
 *this file constructs the transition rate for each event including hopping, track switching, turning, lane changes ect
 */
#include <vector>
#include <numeric>
#include <algorithm>
#include "bundle.h"
#include "parameter.h"

//this file needs to be amended when the transition rules are changed.
//13/12/12 including lane independent boundary rates; track switching for type-1 when meet minus only; steping, neighood lane changes;

// Meet_minus check for whether the stepping of particles with a type at a site in a track will meet a minus of another track
// in Inhomo_wk(), using gridlength rather than gridwidth 16_01_13 updated

using namespace std;
bool Bundle::isMinus(const int site, const int track){
    vector<int>::iterator it;
    int siteM=site+bundleState[track].getOrientation(site);// corrected 20-5-2016 to accouont MT case
    it=find(getTrack(track)->minusEndLocations.begin(),getTrack(track)->minusEndLocations.end(),siteM);
    if(it<bundleState[track].minusEndLocations.end())
        return true;
    else return false;
}

bool Bundle::isPlus(const int site, const int track){
    vector<int>::iterator it;
    it=find(getTrack(track)->plusEndLocations.begin(),getTrack(track)->plusEndLocations.end(),site);
    if(it<bundleState[track].plusEndLocations.end())
        return true;
    else return false;
}

/*
 *the minus/plus end are located at its right hand site.
 meetMinus(site,track,type) return ture if type=1(only for dynei), the site forward to (if right onging)/ the orinial site (if left going)
 is a minsu end of some track.
 This "some track" can be the track as the particle in on.
 */

//27/11/14 add condition for meetMinus to be ture: for dynein, the minus end needs to be a single minus end and opposite orientation to the MT where it was
// also require they are in the bundle region characterized by bundleRegion_1 & bundleRegion_2
bool Bundle::meetMinus(const int site, const int track, const int type){// site---initial site
    bool flag=false;
        int orientation_1=bundleState[track].getOrientation(site);
        int siteM;
        if (type==1){// for dynein
            int track_2=0;
            siteM=site+orientation_1;
            if ((max(siteM,site)<bundleRegion_1 || min(site,siteM)>bundleRegion_2)) flag=false;
            else {
              while (track_2<trackNum && !flag ){
                  int ori_2=0;
                  if (bundleState[track_2].plusEndLocations[0]<bundleState[track_2].minusEndLocations[0])
                      ori_2=1;
                //int ori_2=bundleState[track_2].getOrientation(siteM);// newly added 27/11
                int np=bundleState[track_2].plusEndLocations.size();//
                if (track_2!=track && isMinus(siteM,track_2) && ori_2!=orientation_1 && np==1 && (site==getTrack(track_2)->xCoordinate-1 ||  site==getTrack(track_2)->xCoordinate+getTrack(track_2)->length)){
                //if (isMinus(siteM,track_2) && ori_2!=orientation_1&& np==1 && (site>=getTrack(track_2)->xCoordinate-1 &&  site<=getTrack(track_2)->xCoordinate+getTrack(track_2)->length)){
                //if (isMinus(siteM,track_2) && (site>=getTrack(track_2)->xCoordinate-1 &&  site<=getTrack(track_2)->xCoordinate+getTrack(track_2)->length)){
                //if (isMinus(siteM,track_2) && (site<getTrack(track_2)->xCoordinate || site>=getTrack(track_2)->xCoordinate+getTrack(track_2)->length)){// amdend for one minus with two plus 13_04
                    flag=true;
                    break;
                }
                else track_2++;
              }
            }
        }
        else if(type==0){
            int track_2=0;
            siteM=site+1-orientation_1;
            if ((max(siteM,site)<bundleRegion_1 || min(site,siteM)>bundleRegion_2)) flag=false;
            else {
              while (track_2<trackNum && !flag){
                //if (isMinus(siteM,track_2) && (site>=getTrack(track_2)->xCoordinate-1 &&  site<=getTrack(track_2)->xCoordinate+getTrack(track_2)->length)){
                if (isMinus(siteM,track_2) && (site<getTrack(track_2)->xCoordinate || site>=getTrack(track_2)->xCoordinate+getTrack(track_2)->length)){
                    flag=true;
                    break;
                }
                else track_2++;
              }
            }
        }
    return flag;

}


//27_11_14 add inhomogeneous turning for dynein when meeting plus end of other MT of the same orientation

bool Bundle::meetPlus(const int site, const int track, const int type){// site---initial site
    bool flag=false;

        int orientation_1=bundleState[track].getOrientation(site);
        int siteM;

        if (type==1){// for dynein
            siteM=site+orientation_1;
            int track_2=0;
            if ((max(siteM,site)<bundleRegion_1 || min(site,siteM)>bundleRegion_2)) flag=false;
            else {
              while (track_2<trackNum && !flag){
                int ori_2=bundleState[track_2].getOrientation(siteM);// newly added
                //if (isPlus(siteM,track_2) && (site<getTrack(track_2)->xCoordinate || site>=getTrack(track_2)->xCoordinate+getTrack(track_2)->length)){
                if (isPlus(siteM,track_2) && ori_2==orientation_1 && (site<getTrack(track_2)->xCoordinate || site>=getTrack(track_2)->xCoordinate+getTrack(track_2)->length)){
                    flag=true;
                    break;
                }
                else track_2++;
              }
            }
        }
        else if(type==0){
            int track_2=0;
            if ((max(siteM,site)<bundleRegion_1 || min(site,siteM)>bundleRegion_2)) return false;
            else {
              while (track_2<trackNum && !flag){
                siteM=site+1-orientation_1;
                if (isPlus(siteM,track_2) && (site<getTrack(track_2)->xCoordinate || site>=getTrack(track_2)->xCoordinate+getTrack(track_2)->length)){
                    flag=true;
                    break;
                }
                else track_2++;
              }
            }
       }
       return flag;
}

double Bundle::inhomoWk(const int site, const int track){// site=x/hs to the plus end // for type 0
    int orien=getTrack(track)->getOrientation(site);
    sort(getTrack(track)->plusEndLocations.begin(),getTrack(track)->plusEndLocations.end());
    double pk=transitionRates.forward[0];// inhomogeneity for first type
    if (flagInhomo){
        if(orien==0){
            int i=0;
            while (i<static_cast<int>(getTrack(track)->plusEndLocations.size()) && getTrack(track)->plusEndLocations[i]<site){
                i++;
            }
            int plussite=getTrack(track)->plusEndLocations[i];
            int dyneinnum=getTrack(track)->dyneinNum.at(i);
            int dis;
            if (site>=plussite) dis=site-plussite;
            else dis=plussite-1-site;
            int distance_temp=getTrack(track)->inhomoDistances.at(i);
            if((dis==0 && distance_temp>0)){// && (plussite==0 || plussite==1000)) {//if(dis<=distance_temp && distance_temp>0){
               double tempr=pLoad*double(dyneinnum);//(2*double(dyneinnum)/double(distance_temp)-fluxin/velocity[0]*stepSize-double(dis)/double(distance_temp)*(2*double(dyneinnum)/double(distance_temp)-2*fluxin/velocity[0]*stepSize));
               return tempr;
               //return hs*pLoad*((2*dyneinnum-fluxin*distance_temp/(pk*hs))/distance_temp+fluxin/(pk*hs)-double(dis)*hs/(distance_temp*distance_temp)*(2*dyneinnum-fluxin*distance_temp/(pk*hs)));
            }
            else if(distance_temp==0 && dis<=distance_temp){
                double tempr=transitionRates.turning[0]*pLoad;
                return tempr;
            }
            else if (site<50 || site >= 950)//(!isBundled(site))
                return transitionRates.turning[0];//w_k;
            else return transitionRates.turning[1];
        }
        else if(orien==1){
            int i=getTrack(track)->plusEndLocations.size()-1;
            while (i<static_cast<int>(getTrack(track)->plusEndLocations.size()) && getTrack(track)->plusEndLocations[i]>site){
                i--;
            }
            int plussite=getTrack(track)->plusEndLocations[i];
            int dyneinnum=bundleState[track].dyneinNum[i];
            int dis=abs(site-plussite);
            int distance_temp=getTrack(track)->inhomoDistances.at(i);
            if((dis==0 && distance_temp>0) && (plussite==0 || plussite==1000)){//if (dis<=distance_temp && distance_temp>0){
                double tempr=pLoad*double(dyneinnum);//(2*double(dyneinnum/distance_temp)-fluxin/velocity[0]*stepSize-double(dis)/double(distance_temp)*(2*double(dyneinnum/distance_temp)-2*fluxin/velocity[0]*stepSize));
                return tempr;
                // return hs*pLoad*((2*dyneinnum-fluxin*distance_temp/(pk*hs))/distance_temp+fluxin/(pk*hs)-double(dis)*hs/(distance_temp*distance_temp)*(2*dyneinnum-fluxin*distance_temp/(pk*hs)));
            }
            else if(distance_temp==0 && dis<=distance_temp){
                double tempr=transitionRates.turning[0]*pLoad;
                return tempr;
            }
            else if(site<50 || site >= 950)//(!isBundled(site))
                return transitionRates.turning[0];//w_k;
            else return transitionRates.turning[1];
        }
        else return transitionRates.turning[0];
    }
    else{
        if(isBundled(site))
        return transitionRates.turning[1];//w_k;
        else
            return transitionRates.turning[0];
    }
}

double Bundle::getRateValue(const int type_to_move, const int type_to_be, const vector<int> & loc_1, const vector<int> & loc_2){
    int track_1=loc_1[0];
    int track_2=loc_2[0];
    int site_1=loc_1[2];
    int site_2=loc_2[2];
    int lane_1=loc_1[1];
    int lane_2=loc_2[1];
    int orientation_1=bundleState[track_1].getOrientation(site_1);
    int orientation_2=bundleState[track_2].getOrientation(site_2);

    if  (track_1==track_2){
        int track=track_1;
        int orientation=orientation_1;
        //boundary
        if ((bundleState[track].xCoordinate-1==min(site_1,site_2) || bundleState[track].xCoordinate+bundleState[track].length ==max(site_1,site_2))){
            // injection
            if (type_to_move==0 && ((site_1==bundleState[track].xCoordinate-1 && orientation_2==0)|| (orientation_2==1 && site_1==bundleState[track].xCoordinate+bundleState[track].length))){
                return transitionRates.boundaryIn[type_to_move][track];//alpha_k[track];
            }
            else if (type_to_move==1 && ((site_1==bundleState[track].xCoordinate-1 && orientation_2==1)|| (orientation_2==0 && site_1==bundleState[track].xCoordinate+bundleState[track].length))){
                return transitionRates.boundaryIn[type_to_move][track];//alpha_d[track];
            }
            // exit
            else if (type_to_move==0 && ((site_2==bundleState[track].xCoordinate-1 && orientation_1==1)|| (orientation_1==0 && site_2==bundleState[track].xCoordinate+bundleState[track].length))){
                return transitionRates.boundaryOut[type_to_move][track];//beta_k[track];
            }
            else if (type_to_move==1 && ((site_2==bundleState[track].xCoordinate-1 && orientation_1==0)|| (orientation_1==1 && site_2==bundleState[track].xCoordinate+bundleState[track].length))){
                return transitionRates.boundaryOut[type_to_move][track];//beta_d[track];
            }
            else return 0;
        }
        else{//inner part
            if(orientation==orientation_2 && orientation <2){
            //step forward
               if (type_to_move==type_to_be && abs(site_1-site_2)==1 && lane_1==lane_2){
                  if (type_to_move==0 && ((orientation==0 && site_2>site_1)||(orientation==1 && site_2<site_1))){//kinesin               
                      if (meetMinus(site_1,track,type_to_move)){
                          return  transitionRates.forward[0]*transitionRates.forwardEnd[0];//p_rd;
                      }
                      else if(meetPlus(site_1,track,type_to_move)){
                          return  transitionRates.forward[0]*transitionRates.forwardEnd[1];//p_rd;
                      }
                      else return transitionRates.forward[0];//p_k;
                  }
                  else if (type_to_move==1 && ((orientation==0 && site_2<site_1)||(orientation==1 && site_2>site_1))){//dynein
                      if (meetPlus(site_1,track,type_to_move)){
                          return  transitionRates.forward[1]*transitionRates.forwardEnd[2];
                      }
                      else if (meetMinus(site_1,track,type_to_move)){
                          return  transitionRates.forward[1]*transitionRates.forwardEnd[3];
                      }
                      else  return transitionRates.forward[1];//p_d;
                  }
                  else return 0;
               }
            // type change between 0 and 1
               else if (type_to_move!=type_to_be && site_1==site_2 && lane_1==lane_2){
                  if (type_to_move==0 && type_to_be==1){//kineins ->dynein
                    return inhomoWk(site_1,track);// may site dependent associate with dynein density
                }
                  else if (type_to_move==1 && type_to_be==0){//dynein-->kinesin                      
                      if(isMinus(site_1,track))
                          return transitionRates.turning[4];// inhomogeneity for dynein truning at minus end of its own track
                      else if(meetMinus(site_1,track,type_to_move))
                          return transitionRates.turning[5];// inhomogeneity for dynein truning at minus end of another track
                      else if (meetPlus(site_1,track,type_to_move))//27_11_14 add inhomogeneous turning for dynein when meeting plus end of other MT of the same orientation
                          return transitionRates.turning[6];
                      else if(isBundled(site_1))
                            return transitionRates.turning[3];//w_d;
                      else return transitionRates.turning[2];
                }
                  else return 0;
               }
            // lane changes; keep the type, site forward
               else if (type_to_move==type_to_be && abs(site_1-site_2)==1 && (abs(lane_1-lane_2)%bundleState[track].istate.size()==1 ||abs(lane_1-lane_2)%bundleState[track].istate.size()==bundleState[track].istate.size()-1)){
                    if ((type_to_move==0) && ( ( (orientation==0) && (site_2>site_1) )||( (orientation==1) && (site_2<site_1) ) ) ){//kinesin
                       Unit * target=getUnit(track,lane_1,site_2);
                       if (accumulate(target->occupiedNum.begin(),target->occupiedNum.end(),0)==0 && target->unblocked){//unblocked
                           if (meetMinus(site_1,track,type_to_move)){
                               return  transitionRates.laneChange[0]*transitionRates.forwardEnd[0];
                           }
                           else if(meetPlus(site_1,track,type_to_move)){
                               return  transitionRates.laneChange[0]*transitionRates.forwardEnd[1];
                           }
                           return transitionRates.laneChange[0];//lane_k;
                       }
                       else{
                               if (meetMinus(site_1,track,type_to_move)){
                                   return  transitionRates.laneChange[1]*transitionRates.forwardEnd[0];//p_rd;
                               }
                               else if (meetPlus(site_1,track,type_to_move)){
                                 return  transitionRates.laneChange[1]*transitionRates.forwardEnd[1];//p_rd;
                               }
                               return transitionRates.laneChange[1];
                           }
                   }
                    else if (type_to_move==1 && ((orientation==0 && site_2<site_1)||(orientation==1 && site_2>site_1))){//dynein
                       Unit * target=getUnit(track,lane_1,site_2);
                       if (accumulate(target->occupiedNum.begin(),target->occupiedNum.end(),0)==0 && target->unblocked){//unblocked
                           if (meetPlus(site_1,track,type_to_move)){
                               return  transitionRates.laneChange[2]*transitionRates.forwardEnd[2];
                           }
                           else if(meetMinus(site_1,track,type_to_move)){
                               return  transitionRates.laneChange[2]*transitionRates.forwardEnd[3];
                           }
                           return transitionRates.laneChange[2];//lane_d;
                       }
                       else{
                           if (meetPlus(site_1,track,type_to_move)){
                               return  transitionRates.laneChange[3]*transitionRates.forwardEnd[2];
                           }
                           else if(meetMinus(site_1,track,type_to_move)){
                               return  transitionRates.laneChange[3]*transitionRates.forwardEnd[3];
                           }
                           return transitionRates.laneChange[3];//p_d/2;
                       }
                  }
                    else return 0;
               }
               else return 0;
            }
            else if(orientation!=orientation_2 && orientation_1<2 && orientation_2<2){
                if (type_to_move==1 && type_to_be==0  && (site_1-site_2)*(orientation_1-orientation_2)<0){
                    return transitionRates.forward[2];//cross minus end
                }
                else return 0;
            }

            else return 0;
        }
    }
    // track switch
    else if (max(site_1,site_2)>=bundleRegion_1 && min(site_1,site_2)<=bundleRegion_2) {//track_1!=track_2// assuming maximum one site difference in one step ---bundle->updateRate(n,1);
        //int p1=bundledPosition(site_1,track_1,track_2);
        //int p2=bundledPosition(site_2,track_1,track_2);
        std::vector<int> pp(2,0);
        bundledPosition(pp,site_1,site_2,track_1,track_2);
        int p1=pp[0];
        int p2=pp[1];
        double m=getAvailableLane(track_2,site_2);
        int forwardsite=site_1+getForwardSite(orientation_1,type_to_move);
        bool possibleforward=true;
        Track * track_temp=getTrack(track_1);
        if (forwardsite>=track_temp->xCoordinate && forwardsite<track_temp->length+track_temp->xCoordinate){
            if (track_temp->getOrientation(forwardsite)!=track_temp->getOrientation(site_1))
            possibleforward=false;
        }
        else possibleforward=false;


        if(std::max(p1,p2)==2){// at least one in the interior of bundle
            if (type_to_move==0 && type_to_be==1){
                if ((site_1-site_2)*(orientation_1-orientation_2)>0 || (site_1==site_2 && orientation_1==orientation_2))
                return transitionRates.trackSwitchAtInner[0]/m;
                else return 0;
            }
            else if (type_to_move==0 && type_to_be==0){
                if ((orientation_1!=orientation_2 && site_1==site_2) || (orientation_1==orientation_2 && site_1-site_2==-1+2*orientation_1) )
                return transitionRates.trackSwitchAtInner[1]/m;
                else return 0;
            }
            else if (type_to_move==1 && type_to_be==0){
                if (((orientation_1-orientation_2)*(site_1-site_2)<0) || (site_1==site_2 && orientation_1==orientation_2))
                return transitionRates.trackSwitchAtInner[2]/m;
                else return 0;
            }
            else if (type_to_move==1 && type_to_be==1) {
                if((orientation_1!=orientation_2 && site_1==site_2)||(orientation_1==orientation_2 && site_1-site_2==1-2*orientation_1))
                return transitionRates.trackSwitchAtInner[3]/m;
                else return 0;
            }
            else return 0;
        }
        else if(p1==0 && p2==1 && site_2==forwardsite){// into ends
            if(orientation_1==orientation_2){
                if (type_to_move==0 && type_to_be==0)   return transitionRates.trackSwitchAtEnd[0]/m;
                else if(type_to_move==1 && type_to_be==1) return transitionRates.trackSwitchAtEnd[6]/m;
                else return 0;
            }
            else{
                if (type_to_move==0 && type_to_be==1)   return transitionRates.trackSwitchAtEnd[3]/m;
                else if(type_to_move==1 && type_to_be==0) return transitionRates.trackSwitchAtEnd[9]/m;
                else return 0;
            }
        }
        else if (p1<2 && p2==0 && site_2==forwardsite){// out of ends
        //else if (p1==1 && p2==0 && site_2==site_1+getForwardSite(orientation_1,type_to_move)){// out of ends
            if (orientation_1==orientation_2){
                if(type_to_move==0 && type_to_be==0)   return transitionRates.trackSwitchAtEnd[2]/m;
                else if(type_to_move==1 && type_to_be==1) return transitionRates.trackSwitchAtEnd[8]/m;
                else return 0;
            }
            else{
                if (type_to_move==0 && type_to_be==1)   return transitionRates.trackSwitchAtEnd[4]/m;
                else if(type_to_move==1 && type_to_be==0) return transitionRates.trackSwitchAtEnd[10]/m;
                else return 0;
            }
        }
        else if(p1==1 && p2==1 && !possibleforward){//type switch at ends with no movement forwards//possibleforwd added 15_04
            if(orientation_1==orientation_2){
                if (type_to_move==0 && type_to_be==1)   return transitionRates.trackSwitchAtEnd[1]/m;
                else if(type_to_move==1 && type_to_be==0) return transitionRates.trackSwitchAtEnd[7]/m;
                else return 0;
            }
            else{
                if (type_to_move==0 && type_to_be==0)   return transitionRates.trackSwitchAtEnd[5]/m;
                else if(type_to_move==1 && type_to_be==1) return transitionRates.trackSwitchAtEnd[11]/m;
                else return 0;
            }
        }
        else //no track switching is possible
           return 0;
     }       
}

