#include "mainwindow.h"
#include "worker.h"
#include "parameter.h"
#include <time.h>
#include <iostream>
#include <stdlib.h>// for system command
#include <stdio.h>
#include <fstream>

int step;
double t;
long * pdum;
int tStep;
int kymSize;
time_t start,end;

double t_movie;
using namespace std;

bool flagmovie=true;


int main(int argc, char *argv[])
{
    long idum;
    pdum=&idum;
    idum=-(long)time(NULL);
    start = 0,end = 0;
    time(&start);
    tStep=1;
    kymSize=200;
    t=0;
    step=0;    

    srand (time(NULL));

    if (false){
        char dosline[256];
        char path[256];
        char filename[256];
        FILE *f_pro;
        FILE *f_eedis;
        FILE *f_eestate;
        FILE * f_eedistriscore;
       if(flagdefault)  sprintf(path, "eg7_test39_q%2.2f_pmturn%2.2f_pload%2.2f_default",q,p_mturn1,p_load);
       else  sprintf(path, "eg7_test39_q%2.2f_pmturn%2.2f_pload%2.2f_random",q,p_mturn1,p_load);
        sprintf(dosline,"mkdir %s",path);
        system(dosline);
        char f[256];
        sprintf(f,"hs%2.2f_se%2.2f_sek%2.2f_sed%2.2f.txt",hs,se,sek,sed);
       // f_pro=fopen("test.txt","w");


       // for (int run=0;run<num_of_run;run++){
            tStep=1;
            kymSize=1;
            t=0;
            step=0;
            Worker *worker_server = new Worker();
            bool flag=false;
            char fbundle[256];
            cout<<"bundle region "<<worker_server->bundle->bundleRegion_1<<" "<<worker_server->bundle->bundleRegion_2<<" ";

            sprintf(fbundle,"eg7_test20.txt");
            //worker_server->inputCompleteSetUp(fbundle);
            worker_server->inputGeometrySetUp(fbundle);
            worker_server->bundle->initializeTransitionRates();
           // worker_server->bundle->minRate();
            worker_server->bundle->initializeRateVector();

            /*char fpara[256];
            sprintf(fpara,"para.txt");
            worker_server->inputParameterSetup(fpara);*/
           vector<int> num;
           vector<int> region;
          worker_server->bundle->getBundleRegion(region);
          // worker_server->setOutBundleBlockage();
          cout<<"maxlength "<<worker_server->bundle->maxLength<<endl;
          cout<<"bundle region "<<worker_server->bundle->bundleRegion_1<<" "<<worker_server->bundle->bundleRegion_2<<" ";
           if (flagdirpro){
               char filename_4[256];
               sprintf(filename_4,"%s/pro_%s",path,f);
               f_pro=fopen(filename_4,"w");
           }
           if (flageedis){
               char filename_4[256];
               sprintf(filename_4,"%s/eedis_%s",path,f);
               f_eedis=fopen(filename_4,"w");
           }
           if(flageeprojection){
               char filename_4[256];
               sprintf(filename_4,"%s/EEstate_%s",path,f);
               f_eestate=fopen(filename_4,"w");
           }
           if(flageedistriscore){
               char filename_4[256];
               sprintf(filename_4,"%s/EEDistriScore_%s",path,f);
               f_eedistriscore=fopen(filename_4,"w");
           }
           int STEP=100000;
           int den_t=1;
           double tt=5;
            while(t<tmax && step<STEP){
             //   cout<<"step & time "<<step<<" "<<t<<endl;
                if ((t>=ttran) && !flag){
                    worker_server->initializeCalculate();
                    flag=true;
                    //STEP=step+1;
                }
                if ((t>ttran) && (step%3000==0)) {
                    if (flagdirpro){
                        worker_server->bundle->getDirectMove(num,region);
                        fprintf(f_pro," %2.2f %d %d\n",t,num[0],num[1]);
                    }
                    if (flageedis){
                        worker_server->calculateEEDistance(num,region);
                        for (vector<int>::size_type i=0;i<num.size();i++)
                            fprintf(f_eedis,"%d ", num[i]);
                        fprintf(f_eedis,"\n");
                    }
                    if (flageeprojection){
                        worker_server->projectEE(num);
                        for (vector<int>::size_type i=0;i<num.size();i++){
                            fprintf(f_eestate,"%d ", num[i]);
                            //cout<<num[i]<<" ";
                        }
                        fprintf(f_eestate,"\n");
                    }
                }

                if((step%50==0) & (t>ttran && flageedistriscore)){
                      vector<int>state;
                      worker_server->projectEE(state);
                      double score=worker_server->calculateEvendistributionSore(state);
                      fprintf(f_eedistriscore, "%2.2f %2.2f \n", t,score);
                }
                //cout<<"step & time "<<STEP-step<<" "<<tmax-t<<endl;
                worker_server->oneStep();

            //}
         //   std::cout<<"tdensity counter"<<worker_server->densityCounter<<" density_t "<<den_t<<std::endl;

            if ((t-ttran)>double(den_t)*tt) {
            if (worker_server->bundle->flagDen){
                worker_server->densityCalculate();
                /*sprintf(filename,"den%d%s",den_t,f);
                writeToFile3(worker_server->density,path,filename);*/
                sprintf(filename,"sumden%d%s",den_t,f);
                std::vector<std::vector<double> > sumden;
                worker_server->sumDenCurr(sumden,worker_server->density,1);
                writeToFile2(sumden,path,filename);
                std::cout<<"tdensity counter"<<worker_server->densityCounter<<" density_t "<<den_t<<std::endl;

            }
            if(worker_server->bundle->flagCurrent){
                worker_server->currentCalculate();
                sprintf(filename,"curr%s",f);
                writeToFile3(worker_server->current,path,filename);
                sprintf(filename,"sumcurr%s",f);
                vector<vector<double> > sumcurr;
                worker_server->sumDenCurr(sumcurr,worker_server->current,1);
                writeToFile2(sumcurr,path,filename);

            }

            if (worker_server->bundle->flagRunLengthEnd){
               // std::cout<<worker_server->bundle->runLength.size()<<" "<<worker_server->bundle->runLenghValues.size();
                    sprintf(filename,"run%s",f);
                    writeToFile1(worker_server->bundle->runLenghValues,path,filename);
            }
            if (worker_server->bundle->flagTipSize){
                    sprintf(filename,"tipsize%s",f);
                    writeToFile2(worker_server->bundle->tipSize,path,filename);
            }
           /* if (worker_server->bundle->flagEEDistriScore){
                    sprintf(filename,"eedistriScore%s",f);
                    writeToFile2(worker_server->score_t,path,filename);
            }*/
            //delete worker_server;
            worker_server->initializeCalculate();
            den_t=den_t+1;
            }// end of den_t


            }// end of t<ttran
            if (t>tmax){
                if (flagdirpro) fclose(f_pro);
                if (flageedis) fclose(f_eedis);
                if (flageeprojection) fclose(f_eestate);
                if(flageedistriscore) fclose(f_eedistriscore);
            }

      }

    /* if running in service machine, please comment #include "mainwindow.h" and lines below*/

   QApplication a(argc, argv);
    tStep=1;
    kymSize=600;
    t=0;
    step=0;

    Worker *worker = new Worker();
   // cout<<kymSize<<endl;

/*    char fbundle[256];
    sprintf(fbundle,"eg7_test20.txt");
    worker->inputGeometrySetUp(fbundle);
    worker->bundle->initializeTransitionRates();
    worker->bundle->initializeRateVector();
   vector<int> num;
   vector<int> region;
   worker->bundle->getBundleRegion(region);*/

    MainWindow w;
    w.resize(1000,700);
    w.setWindowTitle("BTSim (Bundled Transport Simulation)");
    w.bindWorker(worker);
    w.setFlags();
    w.setPlusLabels();
    w.setBundle3DBackground();
    //worker->setOutBundleBlockage();

    worker->bundle->tipendInitializeState();
    worker->bundle->initializeTransitionRates();

    worker->bundle->initializeRateVector();

    w.show();
    w.showImage();
    worker->flagKymograph=true;
    worker->initializeCalculate();



    cout<<"kymsize"<<kymSize;
    return a.exec();


    time(&end);
    cout<<"time spent "<<(end-start)<<" second"<<endl;
    system("pause");
    return 0;
}
