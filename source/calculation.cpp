/*
 *this file consists of functions (not belong to any class in this project)
 * which are used in gillespie calculation and density calculation in the bundle class
 */
#include <vector>
#include <iostream>
#include <algorithm>// for find_first_of in vector
#include <numeric>// for partial_sum and accumulate in vector
#include <stdio.h>

using namespace std;
bool comp (double a,double b) {return (a>=b);}

bool checkBundleConnect(const std::vector<std::vector<int> > & tracks){
    int m=tracks.size();
    std::vector<int> vtr(1,0);
    std::vector<int> vremain;
    for (int i=1;i<m;i++)vremain.insert(vremain.end(),i);
    int tr=0;
    while (vtr.size()<m && tr<vtr.size()){
        // search for tracks connected with track[k] and add to the list tr; when finish add k;
        int x0=tracks[tr][0];
        int x1=tracks[tr][1];
        for (std::vector<int>::size_type k=0;k<vremain.size();){
            int j=vremain[k];
            int y0=tracks[j][0];
            int y1=tracks[j][1];
            if ((y1<=x1 && y1>=x0) || (y0>=x0 && y0<=x1)) {
                vtr.push_back(j);
                vremain.erase(vremain.begin()+k);
            }
            else{
                k++;
            }
        }
        tr++;

    }
    if (vremain.empty())return true;
    else return false;
}

int gillespieIndex(const vector<double> & rate, const double & r){
    vector<double> part_sum;
    part_sum.resize(rate.size());
    vector <double>::iterator LIterend;
    LIterend=partial_sum(rate.begin(),rate.end(),part_sum.begin());
    vector<double> r_vec;
    if (rate.size()>0) r_vec.push_back(r*part_sum[rate.size()-1]);
    else{
        cout<<"system is in stationary state as no events could occur";
        system("pause");
    }
    vector<double>::iterator it;
    it=find_first_of(part_sum.begin(), part_sum.end(),r_vec.begin(),r_vec.end(), comp);
    return int(it-part_sum.begin());
}

int maxAllowance(const int & type){
    return 1;
}

bool exclusion(const vector<int> & v_type, const vector<int> & loc_1, const vector<int> & loc_2, const int max){
    if (equal(loc_1.begin(),loc_1.end(),loc_2.begin())){
        return true;
    }
    else{
        if (accumulate(v_type.begin(),v_type.end(),0)>=max)
            return false;
        else return true;
    }
}

int sum(vector<vector<vector<vector<int> > > > & state){
    int s=0;
    for (vector<vector<vector<vector<int> > > >::size_type t=0;t<state.size();++t)
        for (vector<vector<vector<int> > >::size_type i=0;i<state[0].size();++i)
            for (vector<vector<int> > ::size_type j=0;j<state[0][0].size();++j)
                for (vector<int>::size_type k=0;k<state[0][0][0].size();++k)
                    s+=state[t][i][j][k];
    return s;
}

void writeToFile1(const vector<int> & vectors, char * pathname, char *filename){
    FILE * f;
    char filename_4[256];
    sprintf(filename_4,"%s/%s",pathname,filename);
    f=fopen(filename_4,"w");
    for (vector<int>::size_type s=0;s<vectors.size();s++){
        fprintf(f,"%d ",vectors[s]);
    }
    fclose(f);
}

void writeToFile2(const vector<vector<double> > & vectors, char * pathname, char *filename){
    FILE * f_dencur;
    char filename_4[256];
    sprintf(filename_4,"%s/%s",pathname,filename);
    f_dencur=fopen(filename_4,"w");
    for (vector<vector<double> >::size_type x=0;x<vectors.size();x++){
        for (vector<double>::size_type y=0;y<vectors[0].size();y++){
            fprintf(f_dencur,"%f ",vectors[x][y]);
        }
        fprintf(f_dencur,"\n");
    }
    fclose(f_dencur);
}

void writeToFile3(const vector<vector<vector<double> > >& vectors, char * pathname, char *filename){
    FILE * f;
    char filename_4[256];
    sprintf(filename_4,"%s/%s",pathname,filename);
    f=fopen(filename_4,"w");
    for (vector<vector<vector<double> > >::size_type x=0;x<vectors.size();x++){
        for (vector<vector<double> >::size_type y=0;y<vectors[0].size();y++){
            for (vector<double>::size_type s=0;s<vectors[0][0].size();s++){
                fprintf(f,"%f ",vectors[x][y][s]);
            }
            fprintf(f,"\n");
        }
    }
    fclose(f);
}

void evendistribution(vector<double> & distribution, const int n, const int N){
    distribution.clear();
    distribution.resize(n);
    for (int i=0;i<n;i++)
        distribution[i]=double(i*N)/double(n);
}
