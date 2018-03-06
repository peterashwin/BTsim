/*
 *this file consits of member functions doing the display of bundle state or possible transition in update with corresponding transition rates
 */


#include <vector>
#include "bundle.h"

using namespace std;
void Bundle::showSingleType(const int type){
    if (type<typeNum){
        for (int n=0;n<trackNum;n++){
            cout<<"On the "<<n<<"th track";
            bundleState[n].showSingleType(type);
        }
    }
    else {
        cout<<"wrong input of type index, which should be less than the number of types"<<endl;
    }
}

void Bundle::showAllTypes(){
    for (int n=0;n<trackNum;n++){
        cout<<"On the "<<n<<"th track"<<endl;
        bundleState[n].showAllTypes();
    }
}

void Bundle::showAllRates(){
    int n=0;
    for (vector<TransitionEvent>::size_type x=0;x<rate.size();x++){
        vector<int>loc=rate[x].loc1;
        cout<<n<<" type change from: "<<rate[x].typeToMove<<" into: "<<rate[x].typeToBe<<" from location ( ";
        for (vector<int>::size_type i=0;i<rate[x].loc1.size();i++) cout<<rate[x].loc1[i]<<" ";
        cout<<") to ( ";
        for (vector<int>::size_type i=0;i<rate[x].loc2.size();i++) cout<<rate[x].loc2[i]<<" ";
        cout<<") with rate "<<rate[x].rateValue<<endl;//" and type occy "<<Get_cell(loc[0],loc[1],loc[2])->type[0]<<" "<<Get_cell(loc[0],loc[1],loc[2])->type[1]<<endl;
        n=n+1;
    }
}

void Bundle::showSingleRate(const int x){
    vector<int>loc=rate[x].loc1;
    cout<<"type change from: "<<rate[x].typeToMove<<" into: "<<rate[x].typeToBe<<" from location ( ";
    for (vector<int>::size_type i=0;i<rate[x].loc1.size();i++) cout<<rate[x].loc1[i]<<" ";
    cout<<") to ( ";
    for (vector<int>::size_type i=0;i<rate[x].loc2.size();i++) cout<<rate[x].loc2[i]<<" ";
    cout<<") with rate "<<rate[x].rateValue<<" and type occy "<<getUnit(loc[0],loc[1],loc[2])->occupiedNum[0]<<" "<<getUnit(loc[0],loc[1],loc[2])->occupiedNum[1]<<endl;
}

