/*
 * this is the mainwindow souce file. the main window will display the bundle in three ways:
 *the detailed state in the bundle -- bundle view
 * the compressed MT view- the state in each track is represented by one line;
 *kymograph view- the state in the whole bundle is represented by a single line meanwhile prevous states are also displayed
 *
 *the main window also has a menubar which allow to change the geometry structure, the parameters in the system
 *as well as do the calculations of runlength, density, current and save data and files
 */
#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "worker.h"
#include "parameter.h"
#include "mydialog.h"
#include "setting.h"
#include "plotwindow.h"
#include <numeric>// for partial_sum and accumulate in vector

#include "Qpainter.h"

#include <QtWidgets>
#include <QPixmap>
#include <QPainter>
#include <QLabel>
#include <QImage>
#include <QPushButton>
#include <QApplication>
#include <QObject>
#include <QFile>
#include <QTextStream>
#include <QFont>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QGraphicsTextItem>
#include <QPointF>
#include <QVector>



extern long* pdum;
extern double t;
extern int step;
extern int tStep;
extern int kymSize;

using namespace std;



MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{

    ui->setupUi(this);
     this->setStyleSheet("QMainWindow{ background: rgb(200,200,250)}");
    this->createStatusBar();

    assistant = new Assistant;

    setWindow=new Setting();
    setTransWindow=new TransWindow();

    plotWindow=new PlotWindow();
    plusLabel.resize(12);
    bundleImage=QImage(200,26,QImage::Format_RGB16);
    previousKymographImage=QImage(kymSize,200,QImage::Format_RGB16);
    for (int i=0;i<plusLabel.size();i++){
        plusLabel[i]=new QLabel(" ",this);
        plusLabel[i]->setFont(QFont("Helvetica", 14, QFont::Bold));
        plusLabel[i]->setGeometry(0,0,1,1);
        // uncomment the following to include minus labels
        /*minusLabel[i]=new QLabel("1",this);
        minusLabel[i]->setFont(QFont("Helvetica", 14, QFont::Bold));
        minusLabel[i]->setGeometry(0,0,1,1);*/
    }
    imageLabel = new QLabel("",this);
    shrinkLabel =new QLabel("",this);
    kymographLabel=new QLabel("",this);
    timeLabel=new QLabel("",this);
    timeLabel1=new QLabel(" ",this);
    kPercentageLabel=new QLabel("",this);
    imageLabel->setGeometry(0,0,1,1);
    shrinkLabel->setGeometry(0,0,1,1);
    kymographLabel->setGeometry(0,0,1,1);
    timeLabel->setGeometry(0,0,1,1);
    timeLabel1->setGeometry(0,0,1,1);
    kPercentageLabel->setGeometry(0,0,1,1);

    QLabel* aLabel = new QLabel("MT view", this);
    QLabel* bLabel = new QLabel("bundle view", this);
    aLabel->setGeometry(50,80,100,20);
    bLabel->setGeometry(50,280,100,20);

    QPushButton* bstart=new QPushButton("start",this);
    bstart->setGeometry(50,50,100,20);
    QPushButton* bstop=new QPushButton("stop",this);
    bstop->setGeometry(150,50,100,20);
    QPushButton* blotssteps=new QPushButton("do one update",this);
    blotssteps->setGeometry(250,50,100,20);
    QLabel* stepsLabel = new QLabel(" steps per update", this);
    stepsLabel->setGeometry(350,50,100,20);
    nStepEdit=new QLineEdit("1",this);
    nStepEdit->setGeometry(450,50,50,20);

    connect(bstart,SIGNAL(clicked()),this,SLOT(startRun()));
    connect(bstop,SIGNAL(clicked()),this,SLOT(stopRun()));
    connect(blotssteps,SIGNAL(clicked()),this,SLOT(lotsUpdate()));

    createActions();
    createMenus();
    densityView=NULL;
    //worker=NULL;



}

MainWindow::~MainWindow()
{
    delete ui;
    //delete all pointers defined by new
}

void MainWindow::about()
{
   QMessageBox::about(this, "About Application", "The <b>Application</b> shows a simulation of bidirectional transport of two species.");
}

void MainWindow::addBlock(){
    int Nmax=worker->bundle->getTotalSites();
    bool ok;
    int value = QInputDialog::getInt(this, tr("QInputDialog::getInteger()"),
                                 tr("set the number of blockages to be random located: "), 0, 0, Nmax, 1, &ok);//innitial, min, max ,step
    worker->bundle->addRandomBlockages(value,pdum);
    chooseInitialize();
    worker->bundle->initializeRateVector();
    worker->bundle->blockedNum=value;
    if (worker->blocking==true){
        blockChanged();
    }
}

void MainWindow::activation(){// active a single site
    if (activationAct->isChecked()){
        activateTipAct->setChecked(false);
        pauseAct->setChecked(false);
        //runLengthFromBundleAct->setChecked(false);
        runLengthFromEndAct->setChecked(false);
        bool ok;
        Bundle* temp_bundle=worker->bundle;
        int Nmax=temp_bundle->particleNum;
        int value = QInputDialog::getInt(this, tr("QInputDialog::getInteger()"), tr("which particle to label:"), 0, 0, Nmax, 1, &ok);
        for (int tr=0;tr<temp_bundle->trackNum;tr++){
            for (int la=0;la<temp_bundle->getTrack(tr)->width;la++){
                for (int site=0;site<temp_bundle->getTrack(tr)->length;site++){
                    for (int par=0;par<static_cast<int>(temp_bundle->getTrack(tr)->istate[la][site].occupiedParticle.size());par++)
                        if (temp_bundle->getTrack(tr)->istate[la][site].occupiedParticle[par].index==value)
                            temp_bundle->getTrack(tr)->istate[la][site].occupiedParticle[par].labeled=true;
                        else temp_bundle->getTrack(tr)->istate[la][site].occupiedParticle[par].labeled=false;
                }
            }
        }
        showLatticeBackup();
    }
    else{
        Bundle* temp_bundle=worker->bundle;
        for (int tr=0;tr<temp_bundle->trackNum;tr++){
            for (int la=0;la<temp_bundle->getTrack(tr)->width;la++){
                for (int site=0;site<temp_bundle->getTrack(tr)->length;site++){
                    for (int par=0;par<static_cast<int>(temp_bundle->getTrack(tr)->istate[la][site].occupiedParticle.size());par++)
                        temp_bundle->getTrack(tr)->istate[la][site].occupiedParticle[par].labeled=false;
                }
            }
        }
        showLatticeBackup();
    }
}

void MainWindow::activateTip(){
    bool checked=activateTipAct->isChecked();
    if (checked){
        activationAct->setChecked(false);
        pauseAct->setChecked(false);
        //runLengthFromBundleAct->setChecked(false);
        runLengthFromEndAct->setChecked(false);
        bool ok=setLabelRegion();
        if (ok)
        worker->bundle->labelParticles(1,worker->bundle->labelRegion[0],worker->bundle->labelRegion[1]);
        else activateTipAct->setChecked(false);
        showLatticeBackup();
    }
    else{
        worker->bundle->clearLable();
        showLatticeBackup();
    }
}

void MainWindow::addTrack(){
    int Nmax=worker->bundle->trackNum;
    bool ok;
    MyDialog *dialog = new MyDialog;
    dialog->show();
    if (dialog->exec()==QDialog::Accepted){
    }
    DialogPara para;
    dialog->setPara(para,ok);

    if (ok){
        std::vector<int> plustemp=para.plus;
        std::vector<int> vminustemp=para.minus;
        std::vector<std::vector<int> > endlist;
        for (int i=0;i<Nmax;i++){
            int xco=worker->bundle->getTrack(i)->xCoordinate;
            int a[2]={xco,xco+worker->bundle->getTrack(i)->length};
            std::vector<int> temp(a,a+2);
            endlist.insert(endlist.end(),temp);
        }
        std::vector<int> endtemp;
        endtemp.insert(endtemp.end(),plustemp.begin(),plustemp.end());
        endtemp.insert(endtemp.end(),vminustemp.begin(),vminustemp.end());
        sort(endtemp.begin(),endtemp.end());
        if (endtemp.size()>2)    endtemp.erase(endtemp.begin()+1,endtemp.end()-1);
        endlist.insert(endlist.end(),endtemp);
        bool flag=checkBundleConnect(endlist);

        if (flag){
        int tr=para.trackindex;
        int width=para.width;
        //std::vector<int> plustemp=para.plus;
        std::vector<int> distemp=para.distance;
        std::vector<int> dyneintemp=para.dynein;
        int length=max(para.minus.back(),para.plus.back())-min(para.minus[0],para.plus[0]);
        if(width*length>0 && tr<Nmax+1 && distemp.size()==dyneintemp.size() && dyneintemp.size()==plustemp.size()){
            worker->bundle->removeUnaccessibleRegion();
            worker->bundle->trackNum=Nmax+1;
            Track temp(length,width);
        worker->bundle->bundleState.insert(worker->bundle->bundleState.begin()+tr,temp);
       worker->bundle->bundleState[tr].minusEndLocations=para.minus;
        worker->bundle->bundleState[tr].plusEndLocations=para.plus;
        worker->bundle->bundleState[tr].inhomoDistances=para.distance;
        worker->bundle->bundleState[tr].dyneinNum=para.dynein;
       // worker->bundle->bundleState[tr].initializeInhomoPara();
        worker->bundle->bundleState[tr].getLength();
        worker->bundle->bundleState[tr].getXCoordinate();
        worker->bundle->getDimension();

        evendistribution(worker->accum_dis,worker->bundle->maxLength,num_of_particle);
        qDebug()<<"new track created";
        worker->bundle->tipRegion=int(9.5*double(worker->bundle->maxLength)/10);
        worker->bundle->initializeTransitionRates();
        if (worker->bundle->unaccessibility){
            checkBundleAccess();
            worker->bundle->setUnaccessibleRegion();
            worker->bundle->setUnaccessibleBlockage();
        }
        chooseInitialize();

        }

        else if(tr>Nmax){
            QMessageBox msgBox;
            QString s("the track index should be non-zero integer and no more than ");
            s.append(QString::number(Nmax));
            msgBox.setText(s);
            msgBox.exec();
        }
        else if(width*length==0){
            QMessageBox msgBox;
            QString s("Both width and length of a track should be a positive integer!");
           // s.append(QString::number(Nmax));
            msgBox.setText(s);
            msgBox.exec();
        }
        else if ((plustemp.size()!=distemp.size() || (plustemp.size()!=dyneintemp.size()))){
            QMessageBox msgBox;
            QString s("vector size of plus location, distance, and dynein not consistant!");
           // s.append(QString::number(Nmax));
            msgBox.setText(s);
            msgBox.exec();
        }
    }
        else if (!flag){
            QMessageBox msgBox;
            QString s("Bundle Not Connected!");
           // s.append(QString::number(Nmax));
            msgBox.setText(s);
            msgBox.exec();
    }
    }

}

void MainWindow::bindWorker(Worker *wk){
    worker=wk;
    qDebug() << "Mainwindow binding worker "<<worker;
}

void MainWindow::blockChanged(){
    if(showBlockageAct->isChecked()){
        worker->blocking=true;
        show3DEffect();
    }
    else{
        worker->blocking=false;
        setBundle3DBackground();
        if (worker->ending) showPMEnds();
        showLatticeBackup();
    }
}


void MainWindow::blockOutBundle(){
    if (blockOutBundleAct->isChecked())
        worker->setOutBundleBlockage();
    else{
        if (flagdefault)
            worker->bundle->defaultInitializeState();
        else if (flagrandom)
            worker->bundle->randomInitializeState(pdum);
        else if (flagtipend)
            worker->bundle->tipendInitializeState();
        else
            worker->bundle->oneendInitializeState();
    }
}

void MainWindow::breakTrack(){
     vector<vector<int> > breakages;
     //worker->getBreakages(breakages);
    // bool flag=worker->breakhealMT(breakages,pdum,1);
   //  while(!flag)
    //     flag=worker->breakhealMT(breakages,pdum,1);
    // showInitialDisplay();
}

void MainWindow::changeTrack(){
    int Nmax=worker->bundle->trackNum-1;
    bool ok;
    MyDialog *dialog = new MyDialog;
    dialog->show();
    dialog->bindWorkerDialog(worker);
    if (dialog->exec()==QDialog::Accepted){
    }
    DialogPara para;
    dialog->setPara(para,ok);


    if (ok){
        std::vector<int> plustemp=para.plus;
        std::vector<int> vminustemp=para.minus;
        std::vector<std::vector<int> > endlist;
        for (int i=0;i<=Nmax;i++){
            if (i!=para.trackindex){
                int xco=worker->bundle->getTrack(i)->xCoordinate;
                int a[2]={xco,xco+worker->bundle->getTrack(i)->length};
                std::vector<int> temp(a,a+2);
                endlist.insert(endlist.end(),temp);
            }
        }
        std::vector<int> endtemp;
        endtemp.insert(endtemp.end(),plustemp.begin(),plustemp.end());
        endtemp.insert(endtemp.end(),vminustemp.begin(),vminustemp.end());
        sort(endtemp.begin(),endtemp.end());
        if (endtemp.size()>2)    endtemp.erase(endtemp.begin()+1,endtemp.end()-1);
        endlist.insert(endlist.end(),endtemp);
        bool flag=checkBundleConnect(endlist);
        if (flag){
        int tr=para.trackindex;
        int width=para.width;
           if (tr<worker->bundle->trackNum){
            worker->bundle->removeUnaccessibleRegion();
            worker->bundle->bundleState[tr].minusEndLocations=para.minus;
            worker->bundle->bundleState[tr].plusEndLocations=para.plus;
            worker->bundle->bundleState[tr].getLength();
            int length=worker->bundle->bundleState[tr].length;
            if (length>0 && width>0 && para.distance.size()==para.dynein.size() && para.distance.size()==para.plus.size()){
                worker->bundle->bundleState[tr].getXCoordinate();
                worker->bundle->bundleState[tr].istate.resize(width);
                worker->bundle->bundleState[tr].inhomoDistances=para.distance;
                worker->bundle->bundleState[tr].dyneinNum=para.dynein;
                for (int i=0;i<width;i++){
                    if (length>static_cast<int>(worker->bundle->bundleState[tr].istate[i].size())) {
                        vector<int> temp(worker->bundle->typeNum,0);
                        Unit C_temp(temp);
                        worker->bundle->bundleState[tr].istate[i].resize(length,C_temp);//may not need C_temp here as default
                    }
                    else {
                        worker->bundle->bundleState[tr].istate[i].erase(worker->bundle->bundleState[tr].istate[i].begin()+length,worker->bundle->bundleState[tr].istate[i].end());
                    }
                 }
                 worker->bundle->getDimension();
                 evendistribution(worker->accum_dis,worker->bundle->maxLength,num_of_particle);
                 worker->bundle->tipRegion=int(9.5*double(worker->bundle->maxLength)/10);
                 worker->bundle->initializeTransitionRates();
                 if (worker->bundle->particleNum>worker->bundle->getTotalSites())
                     worker->bundle->particleNum=worker->bundle->getTotalSites();
                 for (int i=0;i<min((1+Nmax)*2,10);i++){
                     plusLabel[i]->setGeometry(0,0,1,1);
                 }
                 if (worker->bundle->unaccessibility){
                     checkBundleAccess();
                     worker->bundle->setUnaccessibleRegion();
                     worker->bundle->setUnaccessibleBlockage();
                 }
                 chooseInitialize();
            }
            else if (length<=0 || width<=0){
                QMessageBox msgBox;
                QString s("lenth and width of a track should be a positive integer!");
                msgBox.setText(s);
                msgBox.exec();
            }
            else if (para.distance.size()!=para.dynein.size() || para.distance.size()!=para.plus.size()){
                QMessageBox msgBox;
                QString s("vector size of a plus, distance, dynein should be the same!");
                msgBox.setText(s);
                msgBox.exec();
            }
        }
           else if(tr>Nmax){
            QMessageBox msgBox;
            QString s("the track index should be non-zero integer and no more than ");
            s.append(QString::number(Nmax));
            msgBox.setText(s);
            msgBox.exec();
        }
        }
        else if(!flag){
        QMessageBox msgBox;
        QString s("Bundle Not Connected! Please try again");
       // s.append(QString::number(Nmax));
        msgBox.setText(s);
        msgBox.exec();
    }

    }
}


void MainWindow::chooseInitialize(){
    worker->clearCalculateVar();
    if (defaultInitializerAct->isChecked()){
        oneendInitializerAct->setChecked(false);
        randomInitializerAct->setChecked(false);
        defaultInitialize();
    }
    else if(oneendInitializerAct->isChecked()){
        defaultInitializerAct->setChecked(false);
        randomInitializerAct->setChecked(false);
        oneendInitialize();
    }
    else if(randomInitializerAct->isChecked()){
        defaultInitializerAct->setChecked(false);
        oneendInitializerAct->setChecked(false);
        randomInitialize();
    }
    else if(tipendInitializerAct->isChecked()){
        defaultInitializerAct->setChecked(false);
        randomInitializerAct->setChecked(false);
        tipendInitialize();
    }
    else defaultInitialize();
}

void MainWindow::closeEvent(QCloseEvent *)
{
    delete assistant;
    delete setWindow;
    delete setTransWindow;
    delete plotWindow;
    qDebug()<<"delete pointers defined by new";
}

void MainWindow::check(){
    bool ok=true;
    Bundle* temp_bundle=worker->bundle;
    for (int i=0;i<temp_bundle->trackNum;i++){
        for (int j=0;j<temp_bundle->getTrack(i)->width;j++){
            for (int k=0;k<temp_bundle->getTrack(i)->length;k++){
                int a[3]={i,j,k+temp_bundle->getXCoordinate(i)};
                vector<int> temp_loc(a,a+3);
                vector<Particle> & vparticle=temp_bundle->getTrack(i)->istate[j][k].occupiedParticle;
                if (vparticle.size()>0){
                  //  int id=vparticle[0].index;
                    if ( !equal(temp_loc.begin(),temp_loc.end(),vparticle[0].location.begin()) ){
                    //if (vparticle[0].labeled!=temp_bundle->particles[id].labeled || (!equal(temp_loc.begin(),temp_loc.end(),vparticle[0].location.begin())) || (!equal(vparticle[0].location.begin(),vparticle[0].location.end(),temp_bundle->particles[id].location.begin())) ){
                       ok=false;
                       QMessageBox msgBox;
                       QString s("WRONG! not coinsist!");
                       msgBox.setText(s);
                       msgBox.exec();
                    }
                }
            }
        }
    }
    if (ok){
        QMessageBox msgBox;
        QString s("Coinsist! GOOD!");
        msgBox.setText(s);
        msgBox.exec();
    }
}

void MainWindow::checkBundleAccess(){
    int maxMT=3;
    int tr=0;
    int k=0;
    while(tr<worker->bundle->trackNum && worker->bundle->trackNum>maxMT){
        int xco=worker->bundle->getTrack(tr)->xCoordinate;
       // int len=worker->bundle->getTrack(tr)->length;
        k=0;
        int j=0;
        while (k<maxMT && j<worker->bundle->trackNum){
            int xco1=worker->bundle->getTrack(j)->xCoordinate;
            int len1=worker->bundle->getTrack(j)->length;
            if (j!=tr && xco>=xco1 && xco<xco1+len1) k++;
            j++;
        }
        if (k>=maxMT){
            QMessageBox msgBox;
            QString s("Bundle organization are not valid as there are regions with more than ");
            s.append(QString::number(maxMT));
            s.append(" MTs!");
            msgBox.setText(s);
            msgBox.exec();
            break;
        }
        tr++;
    }
    if (k<maxMT){
        QMessageBox msgBox;
        QString s("Bundle organization are valid !");
        msgBox.setText(s);
        msgBox.exec();
    }
}

void MainWindow::colocalize(){
    if (colocalizeTipAct->isChecked()){
        activateTipAct->setChecked(false);
        activationAct->setChecked(false);
        bool ok=setLabelRegion();
        if (ok){
        worker->bundle->labelParticles(1,worker->bundle->labelRegion[0],worker->bundle->labelRegion[1]);        
        worker->bundle->labelParticles(0,worker->bundle->labelRegion[0],worker->bundle->labelRegion[1]);
        worker->bundle->flagColocalizeTip=true;
        worker->initializeColocalizeCalculation();
        }
        else colocalizeTipAct->setChecked(false);
        showLatticeBackup();
    }
    else{
        worker->bundle->clearLable();
        showLatticeBackup();
        worker->bundle->flagColocalizeTip=false;
        worker->bundle->colocalization.clear();
    }
}

void MainWindow::createActions(){
     saveAct = new QAction(tr("&Start Data"), this);
     saveAct->setShortcuts(QKeySequence::Save);
     saveAct->setStatusTip(tr("Save the kymograph data to disk"));
     connect(saveAct, SIGNAL(triggered()), this, SLOT(saveKymograph()));


     stopAct = new QAction(tr("&Stop Save"), this);
     stopAct->setShortcuts(QKeySequence::Close);
     stopAct->setStatusTip(tr("Stop saving the kymograph data to disk"));
     connect(stopAct, SIGNAL(triggered()), this, SLOT(stopSave()));

     recordBundleAct=new QAction(tr("&Bundle"),this);
     recordBundleAct->setStatusTip(tr("Record bundles into a folder"));
     recordBundleAct->setCheckable(true);
     connect(recordBundleAct,SIGNAL(changed()),this,SLOT(recordBundle()));

     recordMTAct=new QAction(tr("&MTs"),this);
     recordMTAct->setStatusTip(tr("Record MTs into a folder"));
     recordMTAct->setCheckable(true);
     connect(recordMTAct,SIGNAL(changed()),this,SLOT(recordMT()));

     recordKymographAct=new QAction(tr("&Kymograph"),this);
     recordKymographAct->setStatusTip(tr("Record kymographs into a folder"));
     recordKymographAct->setCheckable(true);
     connect(recordKymographAct,SIGNAL(changed()),this,SLOT(recordKymograph()));

     exportAct=new QAction(tr("Export"),this);
     exportAct->setStatusTip(tr("Export bundle structure, bundle state and current transition rate"));
     connect(exportAct,SIGNAL(triggered()),this,SLOT(exportBundleGeometry()));

     importGeomAct=new QAction(tr("Import Geometry"),this);
     importGeomAct->setStatusTip(tr("Import bundle structure with default transition rates and initial state"));
     connect(importGeomAct,SIGNAL(triggered()),this,SLOT(importGeometry()));

     importCompleAct=new QAction(tr("Import Whole Set"),this);
     importCompleAct->setStatusTip(tr("Import bundle structure, bundle state and current transition rate"));
     connect(importCompleAct,SIGNAL(triggered()),this,SLOT(importBundle()));

     exitAct=new QAction(tr("&Exit"),this);
     exitAct->setStatusTip("Exit the program!");
     connect(exitAct,SIGNAL(triggered()),qApp,SLOT(quit()));


     startDensityAct=new QAction(tr("&Start Calculate"),this);
     startDensityAct->setCheckable(true);
     startDensityAct->setStatusTip(tr("Start to calculate density (time-average occupancy) if ticked"));
     connect(startDensityAct,SIGNAL(triggered()),this,SLOT(startDensityCalculate()));

     showDensityAct=new QAction(tr("Show &Density"),this);
     showDensityAct->setStatusTip(tr("Show Density &Plots"));
     connect(showDensityAct,SIGNAL(triggered()),this,SLOT(showDensity()));

     savePlotsAct=new QAction(tr("&Save Plots"),this);
     savePlotsAct->setStatusTip(tr("Save plots of density to jpg file"));
     connect(savePlotsAct,SIGNAL(triggered()),this, SLOT(savePlots()));

     saveDensityAct=new QAction(tr("&Save Data"),this);
     saveDensityAct->setStatusTip(tr("Save the density data into file"));
     connect(saveDensityAct,SIGNAL(triggered()),this,SLOT(saveDensity()));


     startCurrentAct=new QAction(tr("&Start Calculate"),this);
     startCurrentAct->setCheckable(true);
     startCurrentAct->setChecked(false);
     startCurrentAct->setStatusTip(tr("Current calculation current (time-average flux) will turn on if ticked"));
     connect(startCurrentAct,SIGNAL(triggered()),this,SLOT(startCurrentCalculate()));

     showCurrentAct=new QAction(tr("&Show Current"),this);
     showCurrentAct->setStatusTip(tr("Show current plots"));
     connect(showCurrentAct,SIGNAL(triggered()),this,SLOT(showCurrent()));


     saveCurrentAct=new QAction(tr("&Save Data"),this);
     saveCurrentAct->setStatusTip(tr("Save the current data into file"));
     connect(saveCurrentAct,SIGNAL(triggered()),this,SLOT(saveCurrent()));

     runLengthFromEndAct=new QAction(tr("&From End"),this);
//     runLengthFromEndAct=new QAction(tr("&From End"),this);
     runLengthFromEndAct->setStatusTip(tr("Calculate runlength of particles starting from the tip end"));
     runLengthFromEndAct->setCheckable(true);
     runLengthFromEndAct->setChecked(false);
     connect(runLengthFromEndAct,SIGNAL(changed()),this,SLOT(runLengthFromEnd()));
/*
     runLengthFromBundleAct=new QAction(tr("&From Bundle"),this);
     runLengthFromBundleAct->setStatusTip(tr("Calculate runlength of particles starting from the beginning"));
     runLengthFromBundleAct->setCheckable(true);
     runLengthFromBundleAct->setChecked(false);
     //connect(runLengthFromBundleAct,SIGNAL(changed()),this,SLOT(runLengthFromBundle()));
*/
     showDistributionAct=new QAction(tr("&Show Distribution"),this);
     showDistributionAct->setStatusTip(tr("Show distribution of run length"));
     connect(showDistributionAct,SIGNAL(triggered()),this,SLOT(showDistribution()));

     saveRunLengthAct=new QAction(tr("&Save Data"),this);
     saveRunLengthAct->setStatusTip(tr("Save run length data"));
     connect(saveRunLengthAct,SIGNAL(triggered()),this,SLOT(saveRunLength()));

     checkAct=new QAction(tr("Check"),this);
     checkAct->setStatusTip("Check if particles labels are consistent! ");
     connect(checkAct,SIGNAL(triggered()),this,SLOT(check()));

     colocalizeTipAct=new QAction(tr("From Tip"),this);
     colocalizeTipAct->setStatusTip(tr("Start calculating colocalization of dynein-carried EE from tip within a short period"));
     colocalizeTipAct->setCheckable(true);
     colocalizeTipAct->setChecked(false);
     connect(colocalizeTipAct,SIGNAL(changed()),this,SLOT(colocalize()));

     saveColocalizeAct=new QAction(tr("Save"),this);
     saveColocalizeAct->setStatusTip(tr("Save colocalization with dynein density"));
     connect(saveColocalizeAct,SIGNAL(triggered()),this,SLOT(saveColocalize()));

     pauseAct=new QAction(tr("Start Record"),this);
     pauseAct->setCheckable(true);
     pauseAct->setChecked(false);
     pauseAct->setStatusTip(tr("Labelled particles at the tip and start record pauses of labelled particles"));
     connect(pauseAct,SIGNAL(changed()),this,SLOT(pauseCalculate()));

     savePauseAct=new QAction(tr("Save"),this);
     savePauseAct->setStatusTip("Save list of pausing times");
     connect(savePauseAct,SIGNAL(triggered()),this,SLOT(savePause()));

     tipRegionAct=new QAction(tr("&Set Tip Region"),this);
     tipRegionAct->setStatusTip(tr("Set up tip region by an interger for tipsize calculation "));
     connect(tipRegionAct,SIGNAL(triggered()),this,SLOT(setTipRegion()));


     startTipSizeAct=new QAction(tr("&Start Calculation"),this);
     startTipSizeAct->setCheckable(true);
     startTipSizeAct->setStatusTip(tr("Start calculation of tip size"));
     connect(startTipSizeAct,SIGNAL(changed()),this,SLOT(startTipSize()));

     showTipsizeAct=new QAction(tr("&Show Tipsize"),this);
     showTipsizeAct->setStatusTip(tr("Show the tipsize evolution against time"));
     connect(showTipsizeAct,SIGNAL(triggered()),this,SLOT(showTipSize()));

     saveTipSizeAct=new QAction(tr("&Save Data"),this);
     saveTipSizeAct->setStatusTip(tr("Save the tip size data into a file"));
     connect(saveTipSizeAct,SIGNAL(triggered()),this,SLOT(saveTipSize()));

     eeDistanceAct=new QAction(tr("EE &Distance"),this);
     eeDistanceAct->setStatusTip(tr("Show instantaneous EE distance distribution"));
     connect(eeDistanceAct,SIGNAL(triggered()),this,SLOT(showEEdistance()));

     eeNumAct=new QAction(tr("EE &Number"),this);
     eeNumAct->setStatusTip(tr("Show EE number as a function of length to the septum."));
     connect(eeNumAct, SIGNAL(triggered()),this,SLOT(showEENum()));

     startEEDistriScoreAct=new QAction(tr("Start & Dsitribution"),this);
     startEEDistriScoreAct->setCheckable(true);
     //startEEDistriScoreAct->setChecked(false);
     startEEDistriScoreAct->setStatusTip(tr("Start record the time dependence of EE score on distribution difference to even distribution"));
     connect(startEEDistriScoreAct,SIGNAL(triggered()),this,SLOT(startDistributionScore()));

     showEEDistriScoreAct=new QAction(tr("Show &Scores"),this);
     showEEDistriScoreAct->setStatusTip(tr("Show time dependence of EE score on distribution difference to even distribution"));
     connect(showEEDistriScoreAct,SIGNAL(triggered()),this,SLOT(showDistributionScore()));

     saveEEStateAct=new QAction(tr("Save Data"),this);
     saveEEStateAct->setStatusTip(tr("Save projected EE state"));
     connect(saveEEStateAct,SIGNAL(triggered()),this,SLOT());

     particleNumAct=new QAction(tr("&Set #Particle"),this);
     particleNumAct->setStatusTip("Change particle numbers and reinitialize");
     connect(particleNumAct,SIGNAL(triggered()),this,SLOT(particleNumChange()));

     kymographSpeedAct=new QAction(tr("&Set Kymorgraph Speed"),this);
     kymographSpeedAct->setStatusTip("Change time scale in the kymograph view!");
     connect(kymographSpeedAct,SIGNAL(triggered()),this,SLOT(setSpeed()));

     setUpBundleAct=new QAction(tr("Set Up Bundle"),this);
     setUpBundleAct->setStatusTip("Set up a bundle structure randomly");
     connect(setUpBundleAct,SIGNAL(triggered()),this,SLOT(setUpBundle()));

     openBoundaryAct=new QAction(tr("Open Boundary"),this);
     openBoundaryAct->setCheckable(true);
     openBoundaryAct->setChecked(false);
     openBoundaryAct->setStatusTip("Set boundary to be open if check! Please check transition rates for open boundary!");
     connect(openBoundaryAct,SIGNAL(changed()),this,SLOT(setOpenBoundary()));

     addTrackAct=new QAction(tr("&Add Track"),this);
     addTrackAct->setStatusTip("Add one track to the bundle.");
     connect(addTrackAct,SIGNAL(triggered()),this, SLOT(addTrack()));

     removeTrackAct=new QAction(tr("&Remove Track"),this);
     removeTrackAct->setStatusTip("Remove one track from the bundle.");
     connect(removeTrackAct,SIGNAL(triggered()),this, SLOT(removeTrack()));

     changeTrackAct=new QAction(tr("&Change Track"),this);
     changeTrackAct->setStatusTip("Change one track in the bundle.");
     connect(changeTrackAct,SIGNAL(triggered()),this, SLOT(changeTrack()));

     breakTrackAct=new QAction(tr("&Break Track"),this);
     breakTrackAct->setStatusTip("Break or heal MTs randomly.");
     connect(breakTrackAct,SIGNAL(triggered()),this,SLOT(breakTrack()));

     setUnaccessibleAct=new QAction(tr("&Set Unaccesibility"),this);
     setUnaccessibleAct->setCheckable(true);
     setUnaccessibleAct->setStatusTip("Set unaccessible regions in the bundle as blockages.");
     connect(setUnaccessibleAct,SIGNAL(changed()),this,SLOT(setUnaccessible()));

     checkBundleAccessAct=new QAction(tr("&Check Bundle"),this);
     checkBundleAccessAct->setStatusTip("Check bundle organization: if the track numbers in the bundle over the limit");
     connect(checkBundleAccessAct,SIGNAL(triggered()),this,SLOT(checkBundleAccess()));


     geometryShowAct=new QAction(tr("&Show Setup"),this);
     geometryShowAct->setStatusTip("Show the bundle strucutre");
     connect(geometryShowAct,SIGNAL(triggered()),this,SLOT(showSetupWindow()));

     activationAct=new QAction(tr("&Activate Particles"),this);
     activationAct->setCheckable(true);
     activationAct->setChecked(false);
     activationAct->setStatusTip("Show labelled particles in different color");
     connect(activationAct,SIGNAL(changed()),this,SLOT(activation()));

     activateTipAct=new QAction(tr("&Activate Tip"),this);
     activateTipAct->setCheckable(true);
     activateTipAct->setChecked(false);
     activateTipAct->setStatusTip("Activate particles near the tip (right end)");
     connect(activateTipAct,SIGNAL(changed()),this,SLOT(activateTip()));

     blockOutBundleAct=new QAction(tr("&BlockOutBundle"),this);
     blockOutBundleAct->setStatusTip("Get regions out side the bundle neat the two ends blcoked!");
     blockOutBundleAct->setCheckable(true);
     blockOutBundleAct->setChecked(false);
     connect(blockOutBundleAct,SIGNAL(changed()),this,SLOT(blockOutBundle()));

     showBlockageAct=new QAction(tr("&Blockages"),this);
     showBlockageAct->setCheckable(true);
     showBlockageAct->setChecked(false);
     showBlockageAct->setStatusTip("Display of blocked units if checked");
     connect(showBlockageAct,SIGNAL(changed()),this,SLOT(blockChanged()));

     blockNumAct=new QAction(tr("&Set #Blockages"),this);
     blockNumAct->setStatusTip("Change the number of blockages and set blockages randomly. Reinitialze will follow after this action.");
     connect(blockNumAct,SIGNAL(triggered()),this,SLOT(setBlockNum()));

     addBlockAct=new QAction(tr("&Add Blockages"),this);
     addBlockAct->setStatusTip("Add additional blockages to the bundle!");
     connect(addBlockAct,SIGNAL(triggered()),this,SLOT(addBlock()));

     showPMendAct=new QAction(tr("&Plus/Minus Ends"),this);
     showPMendAct->setCheckable(true);
     showPMendAct->setChecked(false);
     showPMendAct->setStatusTip("Display of +/- ends if checked");
     connect(showPMendAct,SIGNAL(changed()),this,SLOT(showPMendChanged()));

     threeDimensionEffctAct=new QAction(tr("&3D Effct"),this);
     threeDimensionEffctAct->setCheckable(true);
     threeDimensionEffctAct->setChecked(true);
     threeDimensionEffctAct->setStatusTip("Display 3D effect for the bundle");
     connect(threeDimensionEffctAct,SIGNAL(triggered()),this,SLOT(show3DEffect()));

     showBundleAct=new QAction(tr("&MT View"),this);
     showBundleAct->setCheckable(true);
     showBundleAct->setChecked(true);
     showBundleAct->setStatusTip("Display the MT view");
     connect(showBundleAct,SIGNAL(triggered()),this,SLOT(showLatticeBackup()));

     showMTAct=new QAction(tr("&Bundle View"),this);
     showMTAct->setStatusTip("Display the bundle view");
     showMTAct->setCheckable(true);
     showMTAct->setChecked(true);
     connect(showMTAct,SIGNAL(triggered()),this,SLOT(showShrink()));

     showKymographAct=new QAction(tr("&Kymograph View"),this);
     showKymographAct->setStatusTip("Display kymograph view.");
     showKymographAct->setCheckable(true);
     showKymographAct->setChecked(true);
     connect(showKymographAct,SIGNAL(triggered()),this,SLOT(showKymographTimeDown()));

    //  inhomoAct=new QAction(tr("&inhomogeneous turning (+)"),this);
    //  inhomoAct->setCheckable(true);
    //  inhomoAct->setChecked(true);
     // inhomoAct->setStatusTip("set inhomogeneous turning rate for kinesin-carried EE near plus ends !");
      //connect(inhomoAct,SIGNAL(triggered()),this, SLOT(inhomoChanged()));

      defaultInitializerAct=new QAction(tr("&Default"),this);
      defaultInitializerAct->setStatusTip("Set particles from the first lane of the first track !");
      defaultInitializerAct->setCheckable(true);
      defaultInitializerAct->setChecked(false);
      connect(defaultInitializerAct,SIGNAL(triggered()),this, SLOT(defaultInitialize()));


      oneendInitializerAct=new QAction(tr("Oneend"),this);
      oneendInitializerAct->setStatusTip("Set particles at the right end of the first track");
      oneendInitializerAct->setCheckable(true);
      oneendInitializerAct->setChecked(false);
      connect(oneendInitializerAct,SIGNAL(triggered()),this,SLOT(oneendInitialize()));

      tipendInitializerAct=new QAction(tr("Tipend"),this);
      tipendInitializerAct->setStatusTip("Set particles at the tip end");
      tipendInitializerAct->setCheckable(true);
      tipendInitializerAct->setChecked(false);
      connect(tipendInitializerAct,SIGNAL(triggered()),this,SLOT(tipendInitialize()));


      randomInitializerAct=new QAction(tr("Random"),this);
      randomInitializerAct->setStatusTip(tr("Set particles randomly"));
      randomInitializerAct->setCheckable(true);
      randomInitializerAct->setChecked(true);
      connect(randomInitializerAct,SIGNAL(triggered()),this,SLOT(randomInitialize()));

      transitionRateAct=new QAction(tr("Transition Rates"),this);
      transitionRateAct->setStatusTip("Allow to modify the transition rates in the system.");
      connect(transitionRateAct,SIGNAL(triggered()),this,SLOT(setTransitionRates()));

      aboutAct=new QAction(tr("&About"),this);
      connect(aboutAct,SIGNAL(triggered()),this,SLOT(about()));

      helpAct=new QAction(tr("&Content"),this);
      connect(helpAct,SIGNAL(triggered()),this,SLOT(showdocumentation()));
}

void MainWindow::creatDir(){
    dirname = QInputDialog::getText(this, "Enter New Foler Name", "Enter name");
    QDir temp;
    temp.mkdir(dirname);
}

void MainWindow::createMenus()
{
    fileMenu=menuBar()->addMenu(tr("&File"));
    dataMenu =fileMenu->addMenu(tr("Save &Data"));
    dataMenu->addAction(saveAct);
    dataMenu->addAction(stopAct);
    imageMenu=fileMenu->addMenu(tr("Save &Image"));
    imageMenu->addAction(recordBundleAct);
    imageMenu->addAction(recordMTAct);
    imageMenu->addAction(recordKymographAct);
    fileMenu->addAction(exportAct);
    importMenu=fileMenu->addMenu(tr("&Import"));
    importMenu->addAction(importGeomAct);
    //importMenu->addAction(importCompleAct);
    fileMenu->addAction(exitAct);


    calculateMenu=menuBar()->addMenu(tr("&Calculations"));
    densityMenu=calculateMenu->addMenu(tr("&Density"));
    densityMenu->addAction(startDensityAct);
    densityMenu->addAction(showDensityAct);
    densityMenu->addAction(savePlotsAct);
    densityMenu->addAction(saveDensityAct);


    currentMenu=calculateMenu->addMenu(tr("&Current"));
    currentMenu->addAction(startCurrentAct);
    currentMenu->addAction(showCurrentAct);
    currentMenu->addAction(saveCurrentAct);

    runlengthMenu=calculateMenu->addMenu(tr("&Run Length"));
    runlengthMenu->addAction(runLengthFromEndAct);
    //runlengthMenu->addAction(runLengthFromBundleAct);// not finished yet
    runlengthMenu->addAction(showDistributionAct);
    runlengthMenu->addAction(saveRunLengthAct);
    runlengthMenu->addAction(checkAct);

    colocalizeMenu=calculateMenu->addMenu(tr("Colocalization"));
    colocalizeMenu->addAction(colocalizeTipAct);
    colocalizeMenu->addAction(saveColocalizeAct);

    pauseMenu=calculateMenu->addMenu(tr("Pause"));
    pauseMenu->addAction(pauseAct);
    pauseMenu->addAction(savePauseAct);

    tipSizeMenu=calculateMenu->addMenu(tr("&Tip Size"));
    tipSizeMenu->addAction(tipRegionAct);
    tipSizeMenu->addAction(startTipSizeAct);
    tipSizeMenu->addAction(showTipsizeAct);
    tipSizeMenu->addAction(saveTipSizeAct);

    projectMenu=calculateMenu->addMenu(tr("Projection"));
    projectMenu->addAction(eeDistanceAct);
    projectMenu->addAction(eeNumAct);
    projectMenu->addAction(startEEDistriScoreAct);
    projectMenu->addAction(showEEDistriScoreAct);
    projectMenu->addAction(saveEEStateAct);

    displayMenu=menuBar()->addMenu(tr("Geometry"));    

    displayMenu->addAction(addTrackAct);
    displayMenu->addAction(removeTrackAct);
    displayMenu->addAction(changeTrackAct);
    displayMenu->addAction(setUpBundleAct);
    //displayMenu->addAction(breakTrackAct);
    displayMenu->addAction(setUnaccessibleAct);
    displayMenu->addSeparator();

    displayMenu->addAction(checkBundleAccessAct);
    displayMenu->addSeparator();

    displayMenu->addAction(geometryShowAct);
    displayMenu->addSeparator();

    displayMenu->addAction(activationAct);
    displayMenu->addAction(activateTipAct);

    //displayMenu->addAction(blockOutBundleAct);
    displayMenu->addAction(showBlockageAct);
    displayMenu->addSeparator();
    displayMenu->addAction(showPMendAct);
    displayMenu->addSeparator();
    displayMenu->addAction(threeDimensionEffctAct);

    displayMenu->addAction(showBundleAct);
    displayMenu->addAction(showMTAct);
    displayMenu->addAction(showKymographAct);
    displayMenu->addAction(openBoundaryAct);

    paraMenu=menuBar()->addMenu(tr("Parameters"));
    //paraMenu->addAction(inhomoAct);
    paraMenu->addAction(particleNumAct);
    paraMenu->addAction(blockNumAct);
    paraMenu->addAction(addBlockAct);
    paraMenu->addAction(kymographSpeedAct);
    paraMenu->addAction(transitionRateAct);

    initializerMenu=menuBar()->addMenu(tr("Initializer"));
    initializerMenu->addAction(defaultInitializerAct);
    initializerMenu->addAction(oneendInitializerAct);
    initializerMenu->addAction(randomInitializerAct);
    initializerMenu->addAction(tipendInitializerAct);

    helpMenu=menuBar()->addMenu(tr("Help"));
    helpMenu->addAction(aboutAct);
    //helpMenu->addAction(helpAct);
}


void MainWindow::createStatusBar()
{
    statusBar()->showMessage("Ready");
}


void MainWindow::defaultInitialize(){
    oneendInitializerAct->setChecked(false);
    randomInitializerAct->setChecked(false);
    tipendInitializerAct->setChecked(false);
    //defaultInitializerAct->setChecked(true);
    statusBar()->showMessage("default Initializer");
    worker->bundle->defaultInitializeState();
    worker->bundle->vkymograph.clear();
    worker->bundle->kymographRecord(kymSize,true);
    worker->bundle->initializeRateVector();
    worker->bundle->flagDen=false;
    tStep=1;
    step=0;
    t=0;
    showInitialDisplay();
    showImage();
    statusBar()->showMessage("Ready");
}
void MainWindow::exportBundleGeometry(){
    QString filename = QFileDialog::getSaveFileName(this, tr("Export bundle"),QDir::currentPath(),tr("Text files (*.txt)"));
    QFile * file=new QFile(filename,this);
    file->open( QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(file);
    QString s;
    Bundle *temp_bundle=worker->bundle;

    //bundle structures
    s.append(QString::number(temp_bundle->trackNum));
    s.append("\n");

    for (int tr=0;tr<temp_bundle->trackNum;tr++){
            s.append(QString::number(temp_bundle->getTrack(tr)->minusEndLocations.size())+" ");
            s.append(QString::number(temp_bundle->getTrack(tr)->plusEndLocations.size())+"\n");

            s.append(QString::number(temp_bundle->getTrack(tr)->width)+" ");
            for (vector<double>::size_type m=0;m<temp_bundle->getTrack(tr)->minusEndLocations.size();m++)
                s.append(QString::number(temp_bundle->getTrack(tr)->minusEndLocations[m])+" ");
            for (vector<double>::size_type p=0;p<temp_bundle->getTrack(tr)->plusEndLocations.size();p++)
                s.append(QString::number(temp_bundle->getTrack(tr)->plusEndLocations[p])+" ");
            for (vector<double>::size_type p=0;p<temp_bundle->getTrack(tr)->inhomoDistances.size();p++)
                s.append(QString::number(temp_bundle->getTrack(tr)->inhomoDistances[p])+" ");
            for (vector<double>::size_type p=0;p<temp_bundle->getTrack(tr)->dyneinNum.size();p++)
                s.append(QString::number(temp_bundle->getTrack(tr)->dyneinNum[p])+" ");
            s.append("\n\n\n");
    }
    out <<s;
}

void MainWindow::exportBundle(){
    QString filename = QFileDialog::getSaveFileName(this, tr("Export bundle"),QDir::currentPath(),tr("Text files (*.txt)"));
    QFile * file=new QFile(filename,this);
    file->open( QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(file);
    QString s;
    Bundle *temp_bundle=worker->bundle;

    //bundle structure
    s.append(QString::number(temp_bundle->trackNum));
    s.append("\n");

    for (int tr=0;tr<temp_bundle->trackNum;tr++){
            s.append(QString::number(temp_bundle->getTrack(tr)->minusEndLocations.size())+" ");
            s.append(QString::number(temp_bundle->getTrack(tr)->plusEndLocations.size())+"\n");

            s.append(QString::number(temp_bundle->getTrack(tr)->width)+" ");
            for (vector<double>::size_type m=0;m<temp_bundle->getTrack(tr)->minusEndLocations.size();m++)
                s.append(QString::number(temp_bundle->getTrack(tr)->minusEndLocations[m])+" ");
            for (vector<double>::size_type p=0;p<temp_bundle->getTrack(tr)->plusEndLocations.size();p++)
                s.append(QString::number(temp_bundle->getTrack(tr)->plusEndLocations[p])+" ");
            for (vector<double>::size_type p=0;p<temp_bundle->getTrack(tr)->inhomoDistances.size();p++)
                s.append(QString::number(temp_bundle->getTrack(tr)->inhomoDistances[p])+" ");
            for (vector<double>::size_type p=0;p<temp_bundle->getTrack(tr)->dyneinNum.size();p++)
                s.append(QString::number(temp_bundle->getTrack(tr)->dyneinNum[p])+" ");
            s.append("\n\n\n");
    }

    // bundle state
    int Big=5;
    for (int tr=0;tr<temp_bundle->trackNum;tr++){
         for (vector<vector<int> >::size_type y=0;y<temp_bundle->bundleState[tr].istate.size();y++){
             for (int ty=0;ty<temp_bundle->typeNum;ty++){
                for (int i=0;i<temp_bundle->bundleState[tr].xCoordinate;i++) s.append("0 ");
                for (vector<int>::size_type x=0;x<temp_bundle->bundleState[tr].istate[y].size();x++){
                    if (temp_bundle->getTrack(tr)->istate[y][x].unblocked)
                        s.append(QString::number(temp_bundle->bundleState[tr].istate[y][x].occupiedNum[ty])+" ");
                    else
                        s.append(QString::number(Big)+" ");
                }
                int m=temp_bundle->bundleState[tr].xCoordinate+temp_bundle->bundleState[tr].length;
                for (int i=0;i<temp_bundle->maxLength-m;i++) s.append("0 ");
                s.append("\n");
            }
        }
         s.append("\n");
    }
 s.append("\n\n");
 // transition rates

    s.append(QString::number(temp_bundle->pLoad)+"\n");
    for (std::vector<double>::size_type k=0;k<temp_bundle->velocity.size();k++)
        s.append(QString::number(temp_bundle->velocity[k])+"\n");
    for (std::vector<double>::size_type k=0;k<temp_bundle->runLengthUnbundled.size();k++)
        s.append(QString::number(temp_bundle->runLengthUnbundled[k])+"\n");
    for (std::vector<double>::size_type k=0;k<temp_bundle->runLengthBundled.size();k++)
        s.append(QString::number(temp_bundle->runLengthBundled[k])+"\n");
    for (std::vector<double>::size_type k=0;k<temp_bundle->transitionRates.forward.size();k++)
        s.append(QString::number(temp_bundle->transitionRates.forward[k])+"\n");
    for (std::vector<double>::size_type k=0;k<temp_bundle->transitionRates.turning.size();k++)
        s.append(QString::number(temp_bundle->transitionRates.turning[k])+"\n");
    for (std::vector<double>::size_type k=0;k<temp_bundle->transitionRates.laneChange.size();k++)
        s.append(QString::number(temp_bundle->transitionRates.laneChange[k])+"\n");
    for (std::vector<double>::size_type k=0;k<temp_bundle->transitionRates.trackSwitchAtEnd.size();k++)
        s.append(QString::number(temp_bundle->transitionRates.trackSwitchAtEnd[k])+"\n");
    for (std::vector<double>::size_type k=0;k<temp_bundle->transitionRates.trackSwitchAtInner.size();k++)
        s.append(QString::number(temp_bundle->transitionRates.trackSwitchAtInner[k])+"\n");
    for (std::vector<double>::size_type k=0;k<temp_bundle->transitionRates.forwardEnd.size();k++)
        s.append(QString::number(temp_bundle->transitionRates.forwardEnd[k])+"\n");
    for (std::vector<double>::size_type k=0;k<temp_bundle->transitionRates.boundaryIn.size();k++)
        s.append(QString::number(temp_bundle->transitionRates.boundaryIn[k][0])+"\n");
    for (std::vector<double>::size_type k=0;k<temp_bundle->transitionRates.boundaryOut.size();k++)
        s.append(QString::number(temp_bundle->transitionRates.boundaryOut[k][0])+"\n");

    out <<s;
}

void MainWindow::importBundle(){
    tStep=1;
  //  kymSize=250;
    t=0;
    step=0;
    worker->initialzieWorkerVar();
    worker->bundle->initializeBundleVar();
    worker->clearCalculateVar();
    setActChecked();
    QString filename = QFileDialog::getOpenFileName(this, tr("Import bundle"),QDir::currentPath(),tr("Text files (*.txt)"));
    QFile * file=new QFile(filename,this);
    file->open( QIODevice::ReadOnly| QIODevice::Text);
    if(true){
        QByteArray ba = filename.toLocal8Bit();
        char *c_file = ba.data();
        worker->inputCompleteSetUp(c_file);
    }
    // uncomment below to use Qstring to import bundle
    else{
     QTextStream in(file);
     Bundle *temp_bundle=worker->bundle;
        QString trStr = in.readLine();
        temp_bundle->trackNum=trStr.toInt();
        temp_bundle->bundleState.clear();
        int tr=0;
        while(tr<temp_bundle->trackNum){
           QString mp = in.readLine();
           while (mp.isEmpty())   mp = in.readLine();
           QStringList list = mp.split(" ", QString::SkipEmptyParts);

           QString trackIndex = in.readLine();
           while(trackIndex.isEmpty())  trackIndex = in.readLine();
           QStringList list2 = trackIndex.split(" ", QString::SkipEmptyParts);

           int m=list[0].toInt()+1;
           std::vector<int> vmp;
           for (int j=0;j<list[1].toInt();j++)
           vmp.insert(vmp.end(),list2[j+m].toInt());
           for (int i=0;i<list[0].toInt();i++)
           vmp.insert(vmp.end(),list2[i+1].toInt());
           int len=*(max_element(vmp.begin(),vmp.end()))-*(min_element(vmp.begin(),vmp.end()));

           Track temp(len,list2[0].toInt());//tempary track with length L and width list2[0].toInt()
           temp_bundle->bundleState.insert(temp_bundle->bundleState.end(),temp);

           //temp_bundle->getTrack(tr)->width=list2[0].toInt();
           temp_bundle->getTrack(tr)->minusEndLocations.clear();
           temp_bundle->getTrack(tr)->plusEndLocations.clear();
           temp_bundle->getTrack(tr)->inhomoDistances.clear();
           temp_bundle->getTrack(tr)->dyneinNum.clear();
           for (int i=0;i<list[0].toInt();i++)
               temp_bundle->getTrack(tr)->minusEndLocations.insert(temp_bundle->getTrack(tr)->minusEndLocations.end(),list2[i+1].toInt());

           for (int j=0;j<list[1].toInt();j++){
               temp_bundle->getTrack(tr)->plusEndLocations.insert(temp_bundle->getTrack(tr)->plusEndLocations.end(),list2[j+m].toInt());
               temp_bundle->getTrack(tr)->inhomoDistances.insert(temp_bundle->getTrack(tr)->inhomoDistances.end(),list2[j+list[1].toInt()+m].toInt());
               temp_bundle->getTrack(tr)->dyneinNum.insert(temp_bundle->getTrack(tr)->dyneinNum.end(),list2[j+list[1].toInt()*2+m].toInt());
           }
           temp_bundle->getTrack(tr)->getLength();
           temp_bundle->getTrack(tr)->getXCoordinate();
           tr++;
        }

       // temp_bundle->bundleState.erase(temp_bundle->bundleState.begin()+temp_bundle->trackNum,temp_bundle->bundleState.end());
        temp_bundle->getDimension();

        // bundle state
        int parNum=0;
        for (int tr=0;tr<temp_bundle->trackNum;tr++){
           // qDebug()<<tr;
            int xco=temp_bundle->getTrack(tr)->xCoordinate;
            int len=temp_bundle->getTrack(tr)->length;
            for (int la=0;la<temp_bundle->getTrack(tr)->width;la++){
                //qDebug()<<la;
                    QString ty0Index = in.readLine();
                    while(ty0Index.isEmpty())  ty0Index = in.readLine();
                    QStringList list0 = ty0Index.split(" ", QString::SkipEmptyParts);

                    QString ty1Index = in.readLine();
                    while(ty1Index.isEmpty())  ty1Index = in.readLine();
                    QStringList list1 = ty1Index.split(" ", QString::SkipEmptyParts);

                    for (int x=0;x<len;x++){
                        int occu0=list0[x+xco].toInt();
                        int occu1=list1[x+xco].toInt();
                        if (occu0==1 && occu1==0){
                            temp_bundle->bundleState[tr].istate[la][x].occupiedNum[0]=1;
                            temp_bundle->bundleState[tr].istate[la][x].occupiedNum[1]=0;
                            Particle particle;
                            bool labeled=false;
                            int a[6]={parNum,labeled,tr,la,x+xco,0};
                            temp_bundle->assignParticle(particle,a);
                            temp_bundle->listOfParticles.insert(temp_bundle->listOfParticles.end(),particle);
                            temp_bundle->bundleState[tr].istate[la][x].occupiedParticle.clear();
                            temp_bundle->bundleState[tr].istate[la][x].occupiedParticle.push_back(particle);
                            parNum++;
                        }
                        else if(occu0==0 && occu1==1){
                            temp_bundle->bundleState[tr].istate[la][x].occupiedNum[0]=0;
                            temp_bundle->bundleState[tr].istate[la][x].occupiedNum[1]=1;
                            Particle particle;
                            bool labeled=false;
                            int a[6]={parNum,labeled,tr,la,x+xco,1};
                            temp_bundle->assignParticle(particle,a);
                            temp_bundle->listOfParticles.insert(temp_bundle->listOfParticles.end(),particle);
                            temp_bundle->bundleState[tr].istate[la][x].occupiedParticle.clear();
                            temp_bundle->bundleState[tr].istate[la][x].occupiedParticle.push_back(particle);
                            parNum++;
                        }
                        else if(occu0==0 && occu1==0){
                            temp_bundle->bundleState[tr].istate[la][x].occupiedNum[0]=0;
                            temp_bundle->bundleState[tr].istate[la][x].occupiedNum[1]=0;
                            temp_bundle->bundleState[tr].istate[la][x].occupiedParticle.clear();
                        }
                        else if (occu0==1 && occu1==1){
                            QMessageBox msgBox;
                            QString s("not satisfy the exclusion principle at site !");
                            s.append(QString::number(x));
                            msgBox.setText(s);
                            msgBox.exec();
                        }
                        else{//blocked
                            temp_bundle->getTrack(tr)->istate[la][x].unblocked=false;
                            temp_bundle->getTrack(tr)->istate[la][x].occupiedNum.clear();
                            temp_bundle->getTrack(tr)->istate[la][x].occupiedParticle.clear();
                        }
                  }
            }
        }
        temp_bundle->copyState();
        temp_bundle->particleNum=parNum;


// set cell property

        // transition rates

       QString rateIndex = in.readLine();
       while(rateIndex.isEmpty())  rateIndex = in.readLine();
       temp_bundle->pLoad=rateIndex.toDouble();

       for (int i=0;i<2;i++){
           rateIndex = in.readLine();
           temp_bundle->velocity.insert(temp_bundle->velocity.end(),rateIndex.toDouble());
       }
       for (int i=0;i<2;i++){
           rateIndex = in.readLine();
           temp_bundle->runLengthUnbundled.insert(temp_bundle->runLengthUnbundled.end(),rateIndex.toDouble());
       }
       for (int i=0;i<2;i++){
           rateIndex = in.readLine();
           temp_bundle->runLengthBundled.insert(temp_bundle->runLengthBundled.end(),rateIndex.toDouble());
       }
       for (int i=0;i<3;i++){// have included the cross minus end on the same track (which is set to be 0)
           rateIndex = in.readLine();
           temp_bundle->transitionRates.forward.insert(temp_bundle->transitionRates.forward.end(),rateIndex.toDouble());
       }
       for (int i=0;i<5;i++){
           rateIndex = in.readLine();
           temp_bundle->transitionRates.turning.insert(temp_bundle->transitionRates.turning.end(),rateIndex.toDouble());
       }
       for (int i=0;i<4;i++){
           rateIndex = in.readLine();
           temp_bundle->transitionRates.laneChange.insert(temp_bundle->transitionRates.laneChange.end(),rateIndex.toDouble());
       }
       for (int i=0;i<12;i++){
           rateIndex = in.readLine();
           temp_bundle->transitionRates.trackSwitchAtEnd.insert(temp_bundle->transitionRates.trackSwitchAtEnd.end(),rateIndex.toDouble());
       }
       for (int i=0;i<4;i++){
           rateIndex = in.readLine();
           temp_bundle->transitionRates.trackSwitchAtInner.insert(temp_bundle->transitionRates.trackSwitchAtInner.end(),rateIndex.toDouble());
       }
       for (int i=0;i<4;i++){
           rateIndex = in.readLine();
           temp_bundle->transitionRates.forwardEnd.insert(temp_bundle->transitionRates.forwardEnd.end(),rateIndex.toDouble());
       }

       rateIndex = in.readLine();
       vector<double> alpha_k(temp_bundle->trackNum,rateIndex.toDouble());
       rateIndex = in.readLine();
       vector<double> alpha_d(temp_bundle->trackNum,rateIndex.toDouble());
       rateIndex = in.readLine();
       vector<double> beta_k(temp_bundle->trackNum,rateIndex.toDouble());
       rateIndex = in.readLine();
       vector<double> beta_d(temp_bundle->trackNum,rateIndex.toDouble());

       temp_bundle->transitionRates.boundaryIn.insert(temp_bundle->transitionRates.boundaryIn.end(),alpha_k);
       temp_bundle->transitionRates.boundaryIn.insert(temp_bundle->transitionRates.boundaryIn.end(),alpha_d);
       temp_bundle->transitionRates.boundaryOut.insert(temp_bundle->transitionRates.boundaryOut.end(),beta_k);
       temp_bundle->transitionRates.boundaryOut.insert(temp_bundle->transitionRates.boundaryOut.end(),beta_d);
       //temp_bundle->minRate();
       temp_bundle->initializeRateVector();
    }
    worker->flagKymograph=true;
    showInitialDisplay();
}

/*void MainWindow::inhomoChanged(){
    if (inhomoAct->isChecked()){
        worker->bundle->flagInhomo=true;
        showImage();
    }
    else{
        worker->bundle->flagInhomo=false;
        showImage();
    }
}
*/

void MainWindow::importGeometry(){
    tStep=1;
    //kymSize=250;
    t=0;
    step=0;
    worker->initialzieWorkerVar();
    worker->bundle->initializeBundleVar();
    worker->clearCalculateVar();
    setActChecked();
    QString filename = QFileDialog::getOpenFileName(this, tr("Import bundle"),QDir::currentPath(),tr("Text files (*.txt)"));
    QFile * file=new QFile(filename,this);
    file->open( QIODevice::ReadOnly| QIODevice::Text);

        QByteArray ba = filename.toLocal8Bit();
        char *c_file = ba.data();
        worker->inputGeometrySetUp(c_file);


        worker->bundle->initializeTransitionRates();
        worker->bundle->initializeRateVector();

        //qDebug()<<worker->previousStep;
        worker->flagKymograph=true;
        for (int i=0;i<10;i++){
            plusLabel[i]->setGeometry(0,0,1,1);
        }
    showInitialDisplay();
}

void MainWindow::initialKymographImage(){//time_size=worker->bundle->vkymograph.size();
    previousKymographImage=QImage(worker->bundle->maxLength,kymSize,QImage::Format_RGB16);
    previousKymographImage.fill(Qt::white);
    for (int i=0;i<previousKymographImage.width();i++){
          // the pixel value can be 0(empty),1(kinesin only),2(dynein only),3(both)
          if (worker->bundle->vkymograph[0][i]>0)
          previousKymographImage.setPixel(i,0,qRgb(255.*(worker->bundle->vkymograph[0][i]%2==1),255.*(worker->bundle->vkymograph[0][i]>1),255.*0));
    }

}

void MainWindow::initialKymographImageTimeDown(){//time_size=worker->bundle->vkymograph.size();
    previousKymographImage=QImage(worker->bundle->maxLength,kymSize,QImage::Format_RGB16);
    previousKymographImage.fill(Qt::black);
    for (int i=0;i<previousKymographImage.width();i++){
          // the pixel value can be 0(empty),1(kinesin only),2(dynein only),3(both)
          if (worker->bundle->vkymograph[0][i]>0)
          previousKymographImage.setPixel(i,kymSize-1,qRgb(255.*(worker->bundle->vkymograph[0][i]%2==1),255.*(worker->bundle->vkymograph[0][i]>1),255.*0));
    }

}
void MainWindow::lotsUpdate()
{
    int n=nStepEdit->text().toInt();

    for (int i=0;i<n;i++){
        statusBar()->showMessage("Running for "+QString::number(n)+" steps!"+" Time elaps "+QString::number(t)+" s");
        worker->oneStep();

    }
    //worker->nsteps(n);
    showImage();
    statusBar()->showMessage("Ready");
}


void MainWindow::oneendInitialize(){
    defaultInitializerAct->setChecked(false);
    randomInitializerAct->setChecked(false);
    tipendInitializerAct->setChecked(false);
    //oneendInitializerAct->setChecked(true);
    statusBar()->showMessage("Reinitialize");
    worker->bundle->oneendInitializeState();
    worker->bundle->vkymograph.clear();
    worker->bundle->kymographRecord(kymSize,true);
    worker->bundle->initializeRateVector();
    worker->bundle->flagDen=false;
    tStep=1;
    step=0;
    t=0;
    showInitialDisplay();
    statusBar()->showMessage("Ready");

}


void MainWindow::particleNumChange(){
    int Nmax=worker->bundle->getTotalSites();//-worker->bundle->blockedNum;
    int initial=worker->bundle->particleNum;
    bool ok;
    QString str(" number of particles: (maximum ");
    str.append(QString::number(Nmax));
    str.append(")");
    int value = QInputDialog::getInt(this, tr("QInputDialog::getInteger()"),str, initial, 0, Nmax, 1, &ok);
    worker->bundle->particleNum=value;
    qDebug()<<"particleNum changed to be "<<worker->bundle->particleNum;
    chooseInitialize();
}


void MainWindow::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing, true);
    QPen myPen(Qt::black, 2, Qt::SolidLine);
    painter.setPen(myPen);

    // time bar =250ms=0.25s
    int len=int(115*2*1.0/kymSize/worker->speed);//kymographLabel(50,280,300*3,250);
    //int len=int(250*1.0/kymSize/worker->speed);//kymographLabel(50,280,300*3,250);
    painter.drawLine(40,410,40,410+len);
    painter.drawText(QRectF(18,410,30,30),"1 s");

    int spa=int(5.0*900.0/worker->bundle->stepSize/worker->bundle->maxLength);// 900 corresponds to the size of kymograph in window
    painter.drawLine(950,350,950-spa,350);
    painter.drawText(QRectF(920,330,35,25),"5 um");
}

void MainWindow::pauseCalculate(){
    if (pauseAct->isChecked()){
        bool ok=setLabelRegion();
        if (ok){
            worker->bundle->flagPause=true;
            worker->bundle->labelParticles(1,worker->bundle->labelRegion[0],worker->bundle->labelRegion[1]);
            worker->bundle->labelParticles(0,worker->bundle->labelRegion[0],worker->bundle->labelRegion[1]);
            worker->initialzePauseCalculation();
        }
        else pauseAct->setChecked(false);
        showLatticeBackup();
    }
    else{
        worker->bundle->flagPause=false;
        worker->bundle->clearLable();
        worker->bundle->pauses.clear();
        showLatticeBackup();
    }
}

void MainWindow::quitMw()
{
    statusBar()->showMessage("Quitting");
    worker->running=false;
    QCoreApplication::processEvents();
    qApp->quit();
}




void MainWindow::randomInitialize(){
    defaultInitializerAct->setChecked(false);
    oneendInitializerAct->setChecked(false);
    tipendInitializerAct->setChecked(false);
   // randomInitializerAct->setChecked(true);
    statusBar()->showMessage("Reinitialize");
    worker->bundle->randomInitializeState(pdum);
    worker->bundle->vkymograph.clear();
    worker->bundle->kymographRecord(kymSize,true);
    worker->bundle->initializeRateVector();
    worker->bundle->flagDen=false;
    tStep=1;
    step=0;
    t=0;

    showInitialDisplay();
    statusBar()->showMessage("Ready");
}


void MainWindow::recordBundle(){
    if (recordBundleAct->isChecked()){
        creatDir();
        worker->recordBundle=true;
    }
    else worker->recordBundle=false;
}

void MainWindow::recordKymograph(QFile* f){// record data
    if (worker->saving==true && worker->running==true){
        QTextStream out(f);
        QString s;
        for (vector<vector<int> >::size_type x=0;x<worker->bundle->vkymograph[0].size();x++){
            s.append(QString::number(worker->bundle->vkymograph[0][x]));
            s.append(" ");
        }
        s.append("\n");
        out <<s;
    }
}


void MainWindow::recordMT(){
    if (recordMTAct->isChecked()){
        creatDir();
        worker->recordMT=true;
    }
    else worker->recordMT=false;
}

void MainWindow::recordKymograph(){
    if (recordKymographAct->isChecked()){
        creatDir();
        worker->recordKymograph=true;
    }
    else worker->recordKymograph=false;
}

void MainWindow::removeTrack(){
    int Nmax=worker->bundle->trackNum-1;
    bool ok;
    QString str("which track to delete: (maximum index ");
    str.append(QString::number(Nmax));
    str.append(")");
    int value = QInputDialog::getInt(this, tr("QInputDialog::getInteger()"), str, 0, 0, Nmax, 1, &ok);
    std::vector<std::vector<int> > endlist;
    for (int i=0;i<=Nmax;i++){
        if (i!=value){
            int xco=worker->bundle->getTrack(i)->xCoordinate;
            int a[2]={xco,xco+worker->bundle->getTrack(i)->length};
            std::vector<int> temp(a,a+2);
            endlist.insert(endlist.end(),temp);
        }
    }
    bool flag=checkBundleConnect(endlist);
    if (ok && flag){
        worker->bundle->removeUnaccessibleRegion();
        worker->bundle->bundleState.erase(worker->bundle->bundleState.begin()+value);
        worker->bundle->trackNum=Nmax;
        worker->bundle->getDimension();        
        evendistribution(worker->accum_dis,worker->bundle->maxLength,num_of_particle);
        worker->bundle->tipRegion=int(9.5*double(worker->bundle->maxLength)/10);
        worker->bundle->initializeTransitionRates();

        for (int i=0;i<min((1+Nmax)*2,10);i++){
            plusLabel[i]->setGeometry(0,0,1,1);
        }
        if (worker->bundle->particleNum>worker->bundle->getTotalSites())
            worker->bundle->particleNum=worker->bundle->getTotalSites();

        if (worker->bundle->unaccessibility){
            checkBundleAccess();
            worker->bundle->setUnaccessibleRegion();
            worker->bundle->setUnaccessibleBlockage();
        }
        chooseInitialize();
    }
    else if (!flag){
        QMessageBox msgBox;
        QString s("Bundle Not Connected!");
        msgBox.setText(s);
        msgBox.exec();
    }
}


void MainWindow::runLengthFromEnd(){
    Bundle *temp_bundle=worker->bundle;

    bool checked=runLengthFromEndAct->isChecked();
    if (checked){
        bool ok=setLabelRegion();
        if (ok){
        temp_bundle->flagRunLengthEnd=true;
        //runLengthFromBundleAct->setChecked(false);
        std::vector<int> plus,minus;
        int tr=temp_bundle->getTipTrack();
        plus=temp_bundle->getTrack(tr)->plusEndLocations;
        minus=temp_bundle->getTrack(tr)->minusEndLocations;
        if (plus[plus.size()-1]==temp_bundle->maxLength)
        temp_bundle->labelParticles(1,temp_bundle->labelRegion[0],temp_bundle->labelRegion[1]);
        else
            temp_bundle->labelParticles(0,temp_bundle->labelRegion[0],temp_bundle->labelRegion[1]);
        //temp_bundle->labelTip(1,per);
        }
        else runLengthFromEndAct->setChecked(false);
        showLatticeBackup();
    }
    else{
        temp_bundle->flagRunLengthEnd=false;
        temp_bundle->runLenghValues.clear();
        temp_bundle->runLength.clear();
        temp_bundle->clearLable();
        showLatticeBackup();
    }
}

/*void MainWindow::runLengthFromBundle(){// not checked yet
    Bundle *temp_bundle=worker->bundle;
    bool checked=runLengthFromBundleAct->isChecked();
    if (checked){
        temp_bundle->flagRunLengthBundle=true;
        runLengthFromEndAct->setChecked(false);
        setLabelRegion();
        temp_bundle->labelParticles(1,temp_bundle->labelRegion[0],temp_bundle->labelRegion[1]);
        showLatticeBackup();
    }
    else{
        temp_bundle->flagRunLengthEnd=false;
        temp_bundle->runLenghValues.clear();
        temp_bundle->runLength.clear();
        temp_bundle->clearLable();
        showLatticeBackup();
    }
}*/

void MainWindow::saveKymograph()
{
    QString filename = QFileDialog::getSaveFileName(this, tr("Save File"),QDir::currentPath(),tr("Text files (*.txt)"));
    f=new QFile(filename,this);
    f->open( QIODevice::WriteOnly | QIODevice::Text);
     worker->saving=true;
     recordKymograph(f);
}


void MainWindow::saveCurrent(){
    worker->currentCalculate();
    for (vector<vector<vector<double> > >::size_type x=0;x<worker->current.size();x++){
        QString filename = QFileDialog::getSaveFileName(this, tr("Save file for current of different types"),QDir::currentPath(),tr("Text files (*.txt)"));
        QFile * file=new QFile(filename,this);
        file->open( QIODevice::WriteOnly | QIODevice::Text);
        QTextStream out(file);
        QString s;
        for (vector<vector<double> >::size_type l=0;l<worker->current[0].size();l++){
            for (vector<double>::size_type si=0;si<worker->current[x][l].size();si++){
                s.append(QString::number(worker->current[x][l][si]));
                 s.append(" ");
            }
            s.append("\n");
        }
        out <<s;
    }
}

void MainWindow::saveColocalize(){
    QString filename = QFileDialog::getSaveFileName(this, tr("Save file for colocalization with dynein"),QDir::currentPath(),tr("Text files (*.txt)"));
    QFile * file=new QFile(filename,this);
    file->open( QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(file);
    QString s;
    for (int i=0;i<2;i++){
        for (int j=0;j<worker->bundle->maxLength;j++){
            s.append(QString::number(worker->bundle->colocalization[i][j]));
            s.append(" ");
        }
        s.append("\n");
    }
    out <<s;
}

void MainWindow::saveDensity(){
    worker->densityCalculate();
    for (vector<vector<vector<double> > >::size_type x=0;x<worker->density.size();x++){
        QString filename = QFileDialog::getSaveFileName(this, tr("Save file for density of different types"),QDir::currentPath(),tr("Text files (*.txt)"));
        QFile * file=new QFile(filename,this);
        file->open( QIODevice::WriteOnly | QIODevice::Text);
        QTextStream out(file);
        QString s;
        for (vector<vector<double> >::size_type l=0;l<worker->density[0].size();l++){
            for (vector<double>::size_type si=0;si<worker->density[x][l].size();si++){
                s.append(QString::number(worker->density[x][l][si]));
                 s.append(" ");
            }
            s.append("\n");
        }
        out <<s;
    }
}


void MainWindow::saveEEState(){
       std::vector<int>state;
       worker->projectEE(state);
       QString filename = QFileDialog::getSaveFileName(this, tr("Save file for EE state"),QDir::currentPath(),tr("Text files (*.txt)"));
       QFile *file=new QFile(filename,this);
       file->open( QIODevice::WriteOnly | QIODevice::Text);
       QTextStream out(file);
       QString s;
       for (int i=0;i<static_cast<int>(state.size());i++){
               s.append(QString::number(state[i]));
       }
       out<<s;
}

void MainWindow::savePause(){
    QString filename = QFileDialog::getSaveFileName(this, tr("Save file for pauses of labeled particles"),QDir::currentPath(),tr("Text files (*.txt)"));
    QFile * file=new QFile(filename,this);
    file->open( QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(file);
    QString s;
    for (int i=0;i<static_cast<int>(worker->bundle->pauses.size());i++){
        for (int j=0;j<static_cast<int>(worker->bundle->pauses[i].size());j++)
            s.append(QString::number(worker->bundle->pauses[i][j])+" ");
        s.append("\n");
    }

    out<<s;
}

void MainWindow::savePlots(){
    if (densityView!=NULL){
      QPixmap pixMap = QPixmap::grabWidget(densityView);
      pixMap.save("density.jpg");
    }
    else{
        QMessageBox msgBox;
        QString s("there is no plots to save! ");
        msgBox.setText(s);
        msgBox.exec();
   }
   /* QImage image(densityScene->sceneRect().size().toSize(), QImage::Format_ARGB32);
    QPainter painter(&image);
    painter.setRenderHint(QPainter::Antialiasing);
    densityScene->render(&painter);
    image.save("file_name.png");*/
}

void MainWindow::saveRunLength(){
    QString filename = QFileDialog::getSaveFileName(this, tr("Save file for recorded tip size"),QDir::currentPath(),tr("Text files (*.txt)"));
    QFile * file=new QFile(filename,this);
    file->open( QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(file);
    QString s;
    for (vector<int>::size_type l=0;l<worker->bundle->runLenghValues.size();l++){
            s.append(QString::number(worker->bundle->runLenghValues[l]));
        s.append("\n");
    }
    out <<s;
}

void MainWindow::saveTipSize(){
        QString filename = QFileDialog::getSaveFileName(this, tr("Save file for recorded tip size"),QDir::currentPath(),tr("Text files (*.txt)"));
        QFile * file=new QFile(filename,this);
        file->open( QIODevice::WriteOnly | QIODevice::Text);
        QTextStream out(file);
        QString s;
        for (vector<vector<double> >::size_type l=0;l<worker->bundle->tipSize.size();l++){
            for (vector<double>::size_type i=0;i<worker->bundle->tipSize[0].size();i++){
                s.append(QString::number(worker->bundle->tipSize[l][i]));
                s.append(" ");
            }
            s.append("\n");
        }
        out <<s;
}

void MainWindow::setActChecked(){
    if(worker->threeDEffect)threeDimensionEffctAct->setChecked(true);//show 3D effect on the bundle
    else threeDimensionEffctAct->setChecked(false);
    if(worker->blocking)showBlockageAct->setChecked(true);
    else showBlockageAct->setChecked(false);//show blocked units in different colors
    if(worker->ending) showPMendAct->setChecked(true);
    else showPMendAct->setChecked(false);//show plus/minus ends in view


    if(worker->recordBundle) recordBundleAct->setChecked(true);
    else recordBundleAct->setChecked(false);// save imagers of bundles
    if(worker->recordMT)recordMTAct->setChecked(true);
    else recordMTAct->setChecked(false);//save images of MT
    if(worker->recordKymograph)recordKymographAct->setChecked(true);
    else recordKymographAct->setChecked(false);//save images of kymograph
    if(worker->bundle->flagRunLengthEnd) runLengthFromEndAct->setChecked(true);
    else runLengthFromEndAct->setChecked(false);

    if (worker->bundle->flagDen) startDensityAct->setChecked(true);
    else startDensityAct->setChecked(false);//flagden;
    if (worker->bundle->flagCurrent) startCurrentAct->setChecked(true);
    else startCurrentAct->setChecked(false);//flagcurrent;
    if (worker->bundle->flagTipSize) startTipSizeAct->setChecked(true);
    else startTipSizeAct->setChecked(false);//flagtipsize;
    if (worker->bundle->unaccessibility) setUnaccessibleAct->setChecked(true);
    else setUnaccessibleAct->setChecked(false);//flagunaccessible;//set unaccessible region as blockages
    if (worker->bundle->flagColocalizeTip) colocalizeTipAct->setChecked(true);
    else colocalizeTipAct->setChecked(false);//flagcolocalizetip;
    if (worker->bundle->flagPause) pauseAct->setChecked(true);
    else pauseAct->setChecked(false);//flagpause;
}

void MainWindow::setBlockNum(){
    int Nmax=worker->bundle->getTotalSites();
    bool ok;
    int value = QInputDialog::getInt(this, tr("QInputDialog::getInteger()"),
                                 tr("set the number of blockages to be random located: "), 0, 0, Nmax, 1, &ok);//innitial, min, max ,step
    worker->bundle->setRandomBlockages(value,pdum);
    chooseInitialize();
    worker->bundle->initializeRateVector();
    worker->bundle->blockedNum=value;
    if (worker->blocking==true){
        blockChanged();
    }
}

void MainWindow::setBundle3DBackground(){
    int mx=worker->bundle->maxLength;//TEST length
    //int mx=20+worker->bundle->maxLength;//TEST length
    //int my=13*6;//worker->bundle->maxWidth;// TEST width 39 for reoraganization for eg7_39.txt
    int my=worker->bundle->maxWidth;// TEST width 39 for reoraganization for eg7_39.txt
    int tracknumber=worker->bundle->trackNum;
    int i=0;
    bundleImage=QImage(mx,my,QImage::Format_RGB16);
    bundleImage.fill(Qt::white);
    for (int tr = 0; tr<tracknumber; tr++){
        int xco=worker->bundle->getXCoordinate(tr);
        int width=worker->bundle->bundleState[tr].istate.size();
        int length=worker->bundle->bundleState[tr].istate[0].size();
       //uncomment below to reorganize the display of bundle for eg7_39.txt and NOTE change accoringly in else below
         /*
        if (tr==0) i=0;
        else if( tr==2) i=26;//13;  for eg7_26
        else if (tr==3) i=13;
        else if (tr==1) i=13;//26;
        else if (tr==4) i=0;
        else if (tr==5) i=13;
        else if (tr==6) i=0;*/
        //else if (tr==7) i=0; // */

     /*   if (tr==0) i=0;
        else if( tr==1) i=26;//13;  for eg7_21
        else if (tr==2) i=0;
        else if (tr==3) i=13;//26;
        else if (tr==4) i=26;//13;
        else if (tr==5) i=0;
        else if (tr==6) i=13;//0;*/
        //else if (tr==7) i=0; // */
        for (int la=0;la<width;la++){
            double c=255-255*sin((double(la)*2.0*4.0/13.0+6.5-4.0)/double(width)*3.1415);
            for(int j=0;j<mx;j++){
                if (j>=xco && j<xco+length){
                    if (threeDimensionEffctAct->isChecked()){
                       // uncomment the following for another background
                        /*double c=0;
                        if (2*(j-xco)<length){
                            c=1/(double(abs(la-width/2))*2.0/double(width))*(double(j-xco)/double(length));
                        }
                        else {
                            c=1/(double(abs(la-double(width)/2.0))*2.0/double(width))*(1.0-double(j-xco)/double(length));
                        }
                        c=min(c*length*3,100.0);*/
                        //bundleImage.setPixel(j+10,i+13,qRgb(int(c),int(c),int(c)));//+3 in order to add white space in image

                        bundleImage.setPixel(j,i,qRgb(int(c),int(c),int(c)));//+3 in order to add white space in image
                    }

                }
                else //if(tr==0 || tr==1 || tr==3) // else
                 //   bundleImage.setPixel(j,i,qRgb(255,255,255));// for white background
                    bundleImage.setPixel(j,i,qRgb(200,200,250));
            }
            i=i+1;
        }
    }
}


void MainWindow::setFlags(){
    if (worker->bundle->flagDen)
        startDensityAct->setChecked(true);
    else startDensityAct->setChecked(false);
    if (worker->bundle->flagCurrent) startCurrentAct->setChecked(true);
    else startCurrentAct->setChecked(false);
    //if (worker->bundle->flagInhomo) inhomoAct->setChecked(true);
    //else inhomoAct->setChecked(false);
    if (worker->bundle->flagRunLengthEnd) runLengthFromEndAct->setChecked(true);
    else runLengthFromEndAct->setChecked(false);
    if (worker->bundle->flagTipSize) startTipSizeAct->setChecked(true);
    else startTipSizeAct->setChecked(false);
    if (worker->bundle->unaccessibility)
        setUnaccessibleAct->setChecked(true);
    else
        setUnaccessibleAct->setChecked(false);
}

bool MainWindow::setLabelRegion(){
    int x1=worker->bundle->maxLength-10;
    int x2=worker->bundle->maxLength;
    if (worker->bundle->labelRegion.size()==2){
        x1=worker->bundle->labelRegion[0];
        x2=worker->bundle->labelRegion[1];
    }

    bool ok;
    QString str(QString::number(x1)+";"+QString::number(x2));
    QString text = QInputDialog::getText(this, tr("QInputDialog::getText()"),
                                         tr("e.g 100;120 (label particles in th region [100,120). )"), QLineEdit::Normal,
                                             str, &ok);
    if (ok){
        //double value = QInputDialog::getDouble(this, tr("QInputDialog::getDouble()"),tr("where label Region starts:"), initial, 0, 1, 3, &ok);
        QStringList list = text.split(";", QString::SkipEmptyParts);
        worker->bundle->labelRegion.clear();
        worker->bundle->labelRegion.push_back(list[0].toInt());
        worker->bundle->labelRegion.push_back(list[1].toInt());
        sort(worker->bundle->labelRegion.begin(),worker->bundle->labelRegion.end());
        //qDebug()<<worker->bundle->labelRegion[0];
    }
    return ok;
}

void MainWindow::setSpeed(){
    bool ok;
    double ini=worker->speed;
    double value = QInputDialog::getDouble(this, tr("QInputDialog::getSpeed()"),
                                 tr("set the speed of displaying kymograph: "), ini, 0.00001,100,4,&ok);//innitial, min, max ,decimals
    worker->speed=value;

}

void MainWindow::setOpenBoundary(){
    if (openBoundaryAct->isChecked()){
        for (int i=0;i<worker->bundle->trackNum;i++)
            worker->bundle->getTrack(i)->openBoundary=true;
    }
    else{
        for (int i=0;i<worker->bundle->trackNum;i++)
            worker->bundle->getTrack(i)->openBoundary=false;
    }
}




void MainWindow::setTipRegion(){
    int initial=worker->bundle->tipRegion;
    int Nmax=worker->bundle->maxLength;
    bool ok;
    int value = QInputDialog::getInt(this, tr("QInputDialog::getInteger()"),tr("where tip Region starts:"), initial, 0, Nmax, 1, &ok);
    worker->bundle->tipRegion=value;
}


void MainWindow::setTransitionRates(){
    setTransWindow->bindWorker(worker);
    setTransWindow->setOldRates();
    setTransWindow->show();
}

void MainWindow::setUnaccessible(){
    if (setUnaccessibleAct->isChecked()){
        worker->bundle->unaccessibility=true;
        worker->bundle->setUnaccessibleRegion();
        worker->bundle->setUnaccessibleBlockage();
        chooseInitialize();
        worker->bundle->initializeRateVector();
        blockChanged();
    }
    else{
        worker->bundle->unaccessibility=false;
        worker->bundle->removeUnaccessibleRegion();        
        chooseInitialize();
        worker->bundle->initializeRateVector();
        blockChanged();
    }
    showImage();
}

void MainWindow::setUpBundle(){
     std::vector<std::vector<int> > randBundle;
     int nb=worker->setRandomBundle(randBundle,pdum,4);// return number of MTs

     worker->flagKymograph=true;
     for (int i=0;i<nb;i++){
         plusLabel[i]->setGeometry(0,0,1,1);
     }
    showInitialDisplay();
}


void MainWindow::setPlusLabels(){
    int size=0;
    int tracknumber=worker->bundle->trackNum;
    for (int tr = 0; tr<tracknumber; tr++){
        int xco=worker->bundle->getXCoordinate(tr);
        if ( worker->bundle->bundleState[tr].getOrientation(xco)==1){
            plusLabel.at(size)->setScaledContents(true);
            plusLabel[size]->setGeometry(20,100+tr*78/tracknumber,10,10);//mx*3,my*3);
            // uncomment the following to show minus labels
            // minusLabel.at(tr)->setScaledContents(true);
           // minusLabel[tr]->setGeometry(20,100+tr*78/tracknumber,10,10);
            size++;
        }
        int le=worker->bundle->bundleState[tr].length;
        if( worker->bundle->bundleState[tr].getOrientation(xco+le-1)==0){
            // uncomment the following to show minus labels
            //minusLabel[tr]->setGeometry(960,100+tr*78/tracknumber,10,10);//mx*3,my*3);
            plusLabel.at(size)->setScaledContents(true);
            plusLabel[size]->setGeometry(960,100+tr*78/tracknumber,10,10);
            size++;
        }
    }
}


void MainWindow::showBlockages(){
    Bundle* temp_bundle=worker->bundle;

    int i=0;
    for (int tr = 0; tr<temp_bundle->trackNum; tr++){
      /*   i=0;
         else if( tr==1) i=26;//13;
         else if (tr==2) i=0;
         else if (tr==3) i=13;//26;
         else if (tr==4) i=26;//13;
         else if (tr==5) i=0;
         else if (tr==6) i=13;//0;*/
        int xco=temp_bundle->getXCoordinate(tr);
        for (int la=0;la<temp_bundle->bundleState[tr].width;la++){
            for(int j=xco;j<xco+temp_bundle->getTrack(tr)->length;j++){
                    if (worker->blocking==true){
                        if (temp_bundle->getTrack(tr)->istate[la][j-xco].unblocked==false)
                       //      bundleImage.setPixel(j+10,i+13,qRgb(100,100,100));
                        bundleImage.setPixel(j,i,qRgb(100,100,100));
                    }
                   // else bundleImage.setPixel(j+10,i+13,qRgb(0,0,0));
                    else bundleImage.setPixel(j,i,qRgb(0,0,0));
            }
            i=i+1;
        }
    }
}


void MainWindow::showCurrent(){
    int mx=worker->bundle->maxLength;
    QVector<double> x(mx);
    QVector<double> counts(mx,0);
    for (int i=0;i<mx;i++){
        x[i]=double(i);
    }
    worker->currentCalculate();

    for (int ty=0;ty<static_cast<int>(worker->current.size());ty++){
        QPen graphPen;
        graphPen.setWidthF(1);
        if (ty==0)
            graphPen.setColor(Qt::red);
        else if (ty==1)
            graphPen.setColor(Qt::green);
        for (int j=0;j<static_cast<int>(worker->current[0].size());j++){
            counts.fromStdVector(worker->current[ty][j]);
            plotWindow->plotState(x,counts,graphPen);
        }
    }

    plotWindow->customPlot->xAxis->setRange(0,mx);
    plotWindow->customPlot->yAxis->setRange(0,1);

    plotWindow->resize(621, 515);
    plotWindow->customPlot->setTitle("Current plot");
    plotWindow->customPlot->xAxis->setLabel(" x ");
    plotWindow->customPlot->yAxis->setLabel("Current");
    plotWindow->customPlot->replot();
    plotWindow->show();

}


void MainWindow::showDensity(){
    worker->densityCalculate();
    int mlength=worker->bundle->maxLength;
    int mx=300; // rectangle size of plot
    int my=mx;

    densityView = new QGraphicsView();
    QGraphicsScene * densityScene = new QGraphicsScene();
    densityView->setScene(densityScene);
    densityScene->setSceneRect(-100, -100, mx+200,my+200 );

    QPen pen(Qt::black);
    for (int ty=0;ty<static_cast<int>(worker->density.size());ty++){
        for (int la=0;la<static_cast<int>(worker->density[ty].size());la++){
            for (int si=0;si<static_cast<int>(worker->density[ty][la].size());si++){
                if (ty==0) pen.setColor(Qt::red);
                else if(ty==1) pen.setColor(Qt::green);
                densityScene->addEllipse(si*mx/mlength,worker->density[ty][la][si]*my, 1,1,pen);
            }
        }
    }
    densityView->scale(1, -1);
    // set axis and tick marks
    int tickNum_x=10;// num of ticks in x axis
    int tickNum_y=10;
    densityScene->addRect(0,0,mx,my);


    for (int i=0;i<=tickNum_x;i++){
        densityScene->addLine(i*mx/tickNum_x,0,i*mx/tickNum_x,-5);
        QGraphicsTextItem *ticks_x=new QGraphicsTextItem(QString::number(double(i)/double(tickNum_x)));
        ticks_x->setPos(i*mx/tickNum_x-8, -2);
        ticks_x->setTransform(densityView->transform().inverted());
       // x->setFont(QFont("Helvetica", 14, QFont::Bold));
        densityScene->addItem(ticks_x);
    }

    for (int j=0;j<=tickNum_y;j++){
        densityScene->addLine(-5,j*my/tickNum_y,0,j*my/tickNum_y);
        QGraphicsTextItem *ticks_y=new QGraphicsTextItem(QString::number(double(j)/double(tickNum_y)));
        ticks_y->setPos(-25,j*my/tickNum_y+10);
       // x->setFont(QFont("Helvetica", 14, QFont::Bold));
        densityScene->addItem(ticks_y);
        ticks_y->setTransform(densityView->transform().inverted());
    }

    // add xlable ylabels
    QString xlabel = QInputDialog::getText(this, "Enter label for x", "Enter label");
    QString ylabel = QInputDialog::getText(this, "Enter label for y", "Enter label");
    QGraphicsTextItem *x=new QGraphicsTextItem(xlabel);
    QGraphicsTextItem *y=new QGraphicsTextItem(ylabel);
    x->setPos(mx/2, -15);
    y->setPos(-55,my/2);
    x->setFont(QFont("Helvetica", 14, QFont::Bold));
    y->setFont(QFont("Helvetica", 14, QFont::Bold));
    densityScene->addItem(x);
    densityScene->addItem(y);

    x->setTransform(densityView->transform().inverted());
    y->setTransform(densityView->transform().inverted());
    y->setRotation(270);
    densityView->show();
}


void MainWindow::showDistribution(){
     int binSize=10;
     int size=worker->bundle->maxLength/binSize;
     QVector<double> bins(size+1);
     QVector<double> counts(size+1,0);
     for (int i=0;i<size;i++){
         bins[i]=double(i*binSize);
     }
     bins[size]=binSize*(size);
     for (int i=0;i<static_cast<int>(worker->bundle->runLenghValues.size());i++){
         int id=int(worker->bundle->runLenghValues[i]/binSize);
         if (id<size)  counts[id]++;
         else qDebug()<<"out of bins";

     }
     for (int i=0;i<size;i++){
         if (worker->bundle->runLenghValues.size()>0){
             counts[i]=counts[i]/double(worker->bundle->runLenghValues.size());
         }
     }
     plotWindow->plotDistribution(bins,counts);

     plotWindow->customPlot->xAxis->setRange(0,bins[size]);
     plotWindow->customPlot->yAxis->setRange(0,1);

     plotWindow->resize(621, 515);


     plotWindow->customPlot->setTitle("Distribution");
     plotWindow->customPlot->xAxis->setLabel("Run length");
     plotWindow->customPlot->yAxis->setLabel("Probability");
     plotWindow->customPlot->replot();
     plotWindow->show();
}


void MainWindow::showdocumentation()
{
    assistant->showDocumentation("content.htm");
}

void MainWindow::showDistributionScore(){
    int s=worker->score_t.size();
    if (s>0){
        QVector<double> bins;
        QVector<double> qdifference;
        for (int i=0;i<s;i++){
            bins.insert(bins.end(),worker->score_t[i][0]);
            qdifference.insert(qdifference.end(),worker->score_t[i][1]);
       }
        std::vector<double> difference=qdifference.toStdVector();

        plotWindow->plotDistribution(bins,qdifference);
        plotWindow->customPlot->xAxis->setRange(0,t);
        double m=*std::max_element(difference.begin(),difference.end());
        plotWindow->customPlot->yAxis->setRange(0,m+10);

        plotWindow->resize(621, 515);

        plotWindow->customPlot->setTitle("history of maximum difference to even distribution");
        plotWindow->customPlot->xAxis->setLabel("time (s)");
        plotWindow->customPlot->yAxis->setLabel("difference");
        plotWindow->customPlot->replot();
        plotWindow->show();
    }
}

void MainWindow::showEEdistance(){
    std::vector<int>region;
    worker->bundle->getBundleRegion(region);
    std::vector<int>dis;
    worker->calculateEEDistance(dis,region);


    int binSize=5;
    int size=worker->bundle->maxLength/binSize/10;
    QVector<double> bins(size+1);
    QVector<double> counts(size+1,0);
    for (int i=0;i<size;i++){
        bins[i]=double(i*binSize);
    }
    bins[size]=binSize*(size);
    for (int i=0;i<static_cast<int>(dis.size());i++){
        int id=int(dis[i]/binSize);
        if (id<size)  counts[id]++;
        else qDebug()<<"out of bins";

    }
    for (int i=0;i<size;i++){
        if (dis.size()>0){
            counts[i]=counts[i]/double(dis.size());
        }
    }
    plotWindow->plotDistribution(bins,counts);

    plotWindow->customPlot->xAxis->setRange(0,bins[size]);
    plotWindow->customPlot->yAxis->setRange(0,1);

    plotWindow->resize(621, 515);


    plotWindow->customPlot->setTitle("Distribution");
    plotWindow->customPlot->xAxis->setLabel("EE distance");
    plotWindow->customPlot->yAxis->setLabel("Probability");
    plotWindow->customPlot->replot();
    plotWindow->show();
}


void MainWindow::showEENum(){
    std::vector<int> state;
    worker->projectEE(state);
    double score=worker->calculateEvendistributionSore(state);

    int binSize=1;
    int size=worker->bundle->maxLength/binSize;
    QVector<double> bins(size+1);
    QVector<double> counts(size+1,0);
    for (int i=0;i<size;i++){
        bins[i]=double(i*binSize);
    }
    bins[size]=binSize*(size);
    for (int i=0;i<static_cast<int>(state.size());i++){

        int id=int(i/binSize);
        if (id<size)  counts[id]=counts[id]+state[i];
        else qDebug()<<"out of bins";
    }
    partial_sum(counts.begin(),counts.end(),counts.begin());

    plotWindow->plotDistribution(bins,counts);

    plotWindow->customPlot->xAxis->setRange(0,bins[size]);
    plotWindow->customPlot->yAxis->setRange(0,worker->bundle->particleNum);

    plotWindow->resize(621, 515);

    plotWindow->customPlot->setTitle("maximum difference to even distribution: "+QString::number(score));
    plotWindow->customPlot->xAxis->setLabel("distance to the left end of the bundle");
    plotWindow->customPlot->yAxis->setLabel("# EEs");
    plotWindow->customPlot->replot();
    plotWindow->show();
}

void MainWindow::showImage()
{
    if (showBundleAct->isChecked()) showLatticeBackup();    //showLattice();
    if (showMTAct->isChecked()) showShrink();
    if (showKymographAct->isChecked()) {
        showKymographTimeDown();
    }
}


void MainWindow::showInitialDisplay(){
    //setPlusLabels();
    setBundle3DBackground();
    initialKymographImage();
    if(worker->ending) showPMEnds();
    if(worker->blocking) showBlockages();
    if ((runLengthFromEndAct->isChecked()) & (worker->bundle->labelRegion.size()==2))
        worker->bundle->labelParticles(1,worker->bundle->labelRegion[0],worker->bundle->labelRegion[1]);
    showImage();
}


void MainWindow::showSetupWindow(){
    setWindow->bindWorker(worker);
    setWindow->setUp();
    setWindow->show();
}




void MainWindow::showKymograph(){
    QString s=QString::number(t,'f',4);
            s.append(" s ");
            s.append(" step: ");
            s.append(QString::number(step));
            s.append(" tip size: ");
            int tip=worker->bundle->calculateTipSize();
            s.append(QString::number(tip));
    timeLabel->setText(s);
    timeLabel->setGeometry(750,550,200,20);

    vector<double> num;
    worker->bundle->getKPercentage(num);
    QString per("#type 0 - Kinesin-carried EE (red): ");
    per.append(QString::number(num[0]));
    per.append("\n#type 1 - Dynein-carried EE (green): ");
    per.append(QString::number(num[1]));
    per.append("\nPercentage of Kinesin-carried EE: ");
    per.append(QString::number(num[2],'f',0));
    per.append("%");
    kPercentageLabel->setText(per);
    kPercentageLabel->setGeometry(50,550,400,60);

    bool show=false;// true- using memcpy
    if(show){
        if (tStep==1){
            worker->previousStep=tStep;
            initialKymographImage();
        }
        else if(tStep>worker->previousStep && (tStep-worker->previousStep)<kymSize){
           uchar *data = previousKymographImage.bits();
           memcpy(data+(tStep-worker->previousStep)*previousKymographImage.bytesPerLine(), data, previousKymographImage.byteCount()-(tStep-worker->previousStep)*previousKymographImage.bytesPerLine());
           for (int j=0;j<tStep-worker->previousStep;j++){
               for (int i=0;i<static_cast<int>(worker->bundle->vkymograph[0].size());i++){
                   if (worker->bundle->vkymograph[j][i]>0)
                   previousKymographImage.setPixel(i,j,qRgb(255.*(worker->bundle->vkymograph[j][i]%2==1),255.*(worker->bundle->vkymograph[j][i]>1),255.*0));
                   else  previousKymographImage.setPixel(i,j,qRgb(255,255,255));
               }
           }
            worker->previousStep=tStep;
        }
        else if((tStep-worker->previousStep)>=kymSize){
            previousKymographImage.fill(Qt::white);
              for (int i=0;i<previousKymographImage.width();i++){
                for (int j=0;j<previousKymographImage.height();j++){
                    // the pixel value can be 0(empty),1(kinesin only),2(dynein only),3(both)
                    if (worker->bundle->vkymograph[j][i]>0)
                    previousKymographImage.setPixel(i,j,qRgb(255.*(worker->bundle->vkymograph[j][i]%2==1),255.*(worker->bundle->vkymograph[j][i]==2),255.*(worker->bundle->vkymograph[j][i]>2)));
                }
              }
             worker->previousStep=tStep;
        }
        else{
            worker->previousStep=tStep;
        }

        kymographLabel->setScaledContents(true);
        kymographLabel->setGeometry(50,280,300*3,250);
        kymographLabel->setPixmap(QPixmap::fromImage(previousKymographImage));
    }
    else{
        int mx=worker->bundle->maxLength;
        if (tStep==1){
            worker->previousStep=tStep;
            initialKymographImage();
        }
        else if(tStep>worker->previousStep && (tStep-worker->previousStep)<kymSize){
            QImage previous=previousKymographImage.copy(0,0,mx,kymSize-(tStep-worker->previousStep));
            QImage newupdate(mx,tStep-worker->previousStep,previousKymographImage.format());
            newupdate.fill(Qt::white);
            for (int j=0;j<tStep-worker->previousStep;j++){
                for (int i=0;i<mx;i++){
                    int c=worker->bundle->vkymograph[j][i];
                    if (c>0){
                       // int rc=255*(c>0);//255.*(worker->bundle->vkymograph[j][i]%2==1)
                        //int gc=255*(c>0);//255.*(worker->bundle->vkymograph[j][i]>1)
                        newupdate.setPixel(i,j,qRgb(255.*(c%2==1),255.*((c%2==0)+(c%3==0)),150.*(c>=10)));
                    }
                }
            }
            QPainter painter(&previousKymographImage);
            painter.drawImage(0,0,newupdate);
            painter.drawImage(0,newupdate.height(),previous);
            worker->previousStep=tStep;
        }
        else if((tStep-worker->previousStep)>=kymSize){
            previousKymographImage.fill(Qt::white);
              for (int i=0;i<previousKymographImage.width();i++){
                for (int j=0;j<previousKymographImage.height();j++){
                    // the pixel value can be 0(empty),1(kinesin only),2(dynein only),3(both)
                    int c=worker->bundle->vkymograph[j][i];
                    if (c>0){
                        //int rc=255*(c>0)-(c>=10)*100;
                       // int gc=255*(c>0)-(c>=20)*100;
                        previousKymographImage.setPixel(i,j,qRgb(255.*(c%2==1),255.*((c%2==0)+(c%3==0)),150.*(c>=10)));
                    }
                }
              }
             worker->previousStep=tStep;
        }
        else{
            worker->previousStep=tStep;
        }
        kymographLabel->setScaledContents(true);
        kymographLabel->setGeometry(50,280,300*3,115);//250);
        kymographLabel->setPixmap(QPixmap::fromImage(previousKymographImage));
        if (worker->recordKymograph){
            QString s;
            s.append(dirname);
            s.append("/kymograph");
            s.append(QString::number(step));
            s.append(".jpg");
            previousKymographImage.scaled(kymographLabel->size(), Qt::IgnoreAspectRatio).save(s,"jpg");
        }

    }
    if (false){
    int mx=worker->bundle->maxLength;
    QImage image4(mx,worker->bundle->vkymograph.size(),QImage::Format_RGB16);
    image4.fill(Qt::white);
    if (worker->bundle->vkymograph.size()>0){
      for (int i=0;i<image4.width();i++){
        for (int j=0;j<image4.height();j++){
            if (worker->bundle->vkymograph[j][i]>0)
            image4.setPixel(i,j,qRgb(255.*(worker->bundle->vkymograph[j][i]%2==1),255.*(worker->bundle->vkymograph[j][i]>1),255.*0));
        }
      }
    }
    kymographLabel->setScaledContents(true);
    kymographLabel->setGeometry(50,280,300*3,250);
    kymographLabel->setPixmap(QPixmap::fromImage(image4));
    if (worker->recordKymograph){
        QString s;
        s.append(dirname);
        s.append("/kymograph");
        s.append(QString::number(step));
        s.append(".jpg");
        image4.scaled(kymographLabel->size(), Qt::IgnoreAspectRatio).save(s,"jpg");
    }
    }
}

void MainWindow::showKymographTimeDown(){
    QString s=QString::number(t,'f',4);
            s.append(" s ");
            s.append(" step: ");
            s.append(QString::number(step));
            s.append(" tip size: ");
            int tip=worker->bundle->calculateTipSize();
            s.append(QString::number(tip));
    timeLabel->setText(s);
    timeLabel->setGeometry(700,330,200,20);

    vector<double> num;
    worker->bundle->getKPercentage(num);
    QString per("#type 0 - Kinesin-carried EE (red): ");
    per.append(QString::number(num[0]));
    per.append("\n#type 1 - Dynein-carried EE (green): ");
    per.append(QString::number(num[1]));
    per.append("\nPercentage of Kinesin-carried EE: ");
    per.append(QString::number(num[2],'f',0));
    per.append("%");
    kPercentageLabel->setText(per);
    kPercentageLabel->setGeometry(50,330,400,60);

    bool show=true;// true- using memcpy
    if(show){
        int mx=worker->bundle->maxLength;
        if (tStep==1){
            worker->previousStep=tStep;
            initialKymographImageTimeDown();
        }
        else if(tStep>worker->previousStep && (tStep-worker->previousStep)<kymSize){
            QImage previous=previousKymographImage.copy(0,tStep-worker->previousStep,mx,kymSize-(tStep-worker->previousStep));
            QImage newupdate(mx,tStep-worker->previousStep,previousKymographImage.format());
            newupdate.fill(Qt::black);
            for (int j=tStep-worker->previousStep-1;j>=0;j--){
                for (int i=0;i<mx;i++){
                    int c=worker->bundle->vkymograph[j][i];
                    if (c>0){
                       // int rc=255*(c>0);//255.*(worker->bundle->vkymograph[j][i]%2==1)
                        //int gc=255*(c>0);//255.*(worker->bundle->vkymograph[j][i]>1)
                        newupdate.setPixel(i,j,qRgb(255.*(c>0),255.*(c>0),255.*(c>0))); // for whilte signals
                        //newupdate.setPixel(i,j,qRgb(255.*(c%2==1),255.*((c%2==0)+(c%3==0)),150.*(c>=10)));
                    }
                }
            }
            QPainter painter(&previousKymographImage);
            painter.drawImage(0,kymSize-newupdate.height(),newupdate);
            painter.drawImage(0,0,previous);
            worker->previousStep=tStep;
        }
        else if((tStep-worker->previousStep)>=kymSize){
            previousKymographImage.fill(Qt::white);
              for (int i=0;i<previousKymographImage.width();i++){
                for (int j=previousKymographImage.height()-1;j>0;j--){
                    // the pixel value can be 0(empty),1(kinesin only),2(dynein only),3(both)
                    int c=worker->bundle->vkymograph[j][i];
                    if (c>0){
                        //int rc=255*(c>0)-(c>=10)*100;
                       // int gc=255*(c>0)-(c>=20)*100;
                        previousKymographImage.setPixel(i,j,qRgb(255.*(c>0),255.*(c>0),255.*(c>0))); // for whilte signals
                        //previousKymographImage.setPixel(i,j,qRgb(255.*(c%2==1),255.*((c%2==0)+(c%3==0)),150.*(c>=10)));
                    }
                }
              }
             worker->previousStep=tStep;
        }
        else{
            worker->previousStep=tStep;
        }
        kymographLabel->setScaledContents(true);
        kymographLabel->setGeometry(50,400,300*3,115*2);//115*2);//250);
        kymographLabel->setPixmap(QPixmap::fromImage(previousKymographImage));
        if (worker->recordKymograph){
            QString s;
            s.append(dirname);
            s.append("/kymograph");
            s.append(QString::number(step));
            s.append(".tif");
            previousKymographImage.scaled(kymographLabel->size(), Qt::IgnoreAspectRatio).save(s,"tif");
        }

    }
 }


void MainWindow::showLattice()
{
    int mx=worker->bundle->maxLength;//TEST length
    int my=worker->bundle->maxWidth;// TEST width
    int tracknumber=worker->bundle->trackNum;
    QImage image(mx,my,QImage::Format_RGB16);//image width(length), image height(width)
    image.fill(Qt::black);
    int i=0;
    for (int tr = 0; tr<tracknumber; tr++){
        int xco=worker->bundle->getXCoordinate(tr);
        int width=worker->bundle->bundleState[tr].istate.size();
        int length=worker->bundle->bundleState[tr].istate[0].size();
        for (int la=0;la<width;la++){
            for(int j=0;j<mx;j++){
                if (j>=xco && j<xco+length){
                   double c=worker->bundle->inhomoWk(j,tr)/worker->bundle->transitionRates.turning[0];
                    int r=worker->bundle->bundleState[tr].istate[la][j-xco].occupiedNum[0];
                    int g=worker->bundle->bundleState[tr].istate[la][j-xco].occupiedNum[1];
                    int value=255;
                    if(r+g>0)
                        //image.setPixel(j,i,qRgb(value*(r>0),value*(g>0),0));
                       image.setPixel(j,i,qRgb(int(c)*(r==0)+value*(r>0),int(c)*(g==0)+value*(g>0),int(c)%value));
                    }
                else image.setPixel(j,i,qRgb(200,200,250));
            }
            i=i+1;
        }
    }

    imageLabel->setScaledContents(true);
    imageLabel->setGeometry(50,100,300*3,26*3);//mx*3,my*3);
    imageLabel->setPixmap(QPixmap::fromImage(image));

    if (worker->recordBundle){
        QString s;
        s.append(dirname);
        s.append("/bundle");
        s.append(QString::number(step));
        s.append(".jpg");
        image.scaled(imageLabel->size(), Qt::IgnoreAspectRatio).save(s,"jpg");//does not support gif
    }
}

void MainWindow::showLatticeBackup(){
    QString s=QString::number(t,'f',4);
            s.append(" s ");
            //s.append(" step: ");
            //s.append(QString::number(step));
            //s.append(" tip size: ");
            //int tip=worker->bundle->calculateTipSize();
            //s.append(QString::number(tip));
    //timeLabel1->setText(s);
    //timeLabel1->setGeometry(800,100,80,20);

  //  QPainter* painter = new QPainter();
   // painter->setPen(Qt::blue);
   // painter->setFont(QFont("Arial", 30));
   // painter->drawText(rect(), Qt::AlignCenter, "Text on Image");

    QImage image(bundleImage);
    int mx=worker->bundle->maxLength;//TEST length
    int tracknumber=worker->bundle->trackNum;
    int i=0;
    for (int tr = 0; tr<tracknumber; tr++){
        int xco=worker->bundle->getXCoordinate(tr);
        int width=worker->bundle->bundleState[tr].istate.size();
        int length=worker->bundle->getTrack(tr)->istate[0].size();
        for (int la=0;la<width;la++){
            for(int j=0;j<mx;j++){
                if (j>=xco && j<xco+length){
                    int r=worker->bundle->bundleState[tr].istate[la][j-xco].occupiedNum[0];
                    int g=worker->bundle->bundleState[tr].istate[la][j-xco].occupiedNum[1];
                    int value=255;
                    int temp=0;
                    if(r+g>0)  {
                        temp=255*(worker->bundle->getUnit(tr,la,j)->occupiedParticle[0].labeled);
                        image.setPixel(j,i,qRgb(value*(r>0),value*(g>0),temp));
                        //qDebug()<<"check out of region "<<tr<<" "<<i;
                    }
                }
            }
            i=i+1;
        }
    }
    //qDebug()<<image.pixel(2,1)<<" "<<image.pixel(3,1)<<" "<<image.pixel(4,1)<<image.pixel(1,1);
    imageLabel->setScaledContents(true);
    imageLabel->setGeometry(40,100,300*3,120);//(50,100,300*3,26*3);//mx*3,my*3);
    imageLabel->setPixmap(QPixmap::fromImage(image));

    if (worker->recordBundle){
        QString s;
        s.append(dirname);
        s.append("/bundle");
        s.append(QString::number(step));
        s.append(".tif");
        image.scaled(imageLabel->size(), Qt::IgnoreAspectRatio).save(s,"tif");//does not support gif
    }
}
/*
void MainWindow::showLatticeBackup(){
    QString s=QString::number(t,'f',3);
      //      s.append(" s ");
            //s.append(" step: ");
            //s.append(QString::number(step));
            //s.append(" tip size: ");
            //int tip=worker->bundle->calculateTipSize();
            //s.append(QString::number(tip));
    //timeLabel1->setText(s);
    //timeLabel1->setGeometry(800,100,80,20);

  //  QPainter* painter = new QPainter();
   // painter->setPen(Qt::blue);
   // painter->setFont(QFont("Arial", 30));
   // painter->drawText(rect(), Qt::AlignCenter, "Text on Image");

    QImage image(bundleImage);
    int mx=worker->bundle->maxLength;//TEST length
    int tracknumber=worker->bundle->trackNum;
    int i=0;
    for (int tr = 0; tr<tracknumber; tr++){
        int xco=worker->bundle->getXCoordinate(tr);
        int width=worker->bundle->bundleState[tr].istate.size();
        int length=worker->bundle->getTrack(tr)->istate[0].size();
        //uncomment below to reorganzie the display of bundle for eg67_39.txt
         /*
        if (tr==0) i=0;
        else if( tr==2) i=26;//13;  for eg7_26
        else if (tr==3) i=13;
        else if (tr==1) i=13;//26;
        else if (tr==4) i=0;
        else if (tr==5) i=13;
        else if (tr==6) i=0; */
        //else if (tr==7) i=0; */

  /*     if (tr==0) i=0;
        else if( tr==1) i=26;//13;  for eg7_21
        else if (tr==2) i=0;
        else if (tr==3) i=13;//26;
        else if (tr==4) i=26;//13;
        else if (tr==5) i=0;
        else if (tr==6) i=13;//0;
        //else if (tr==7) i=0; // */
    /*    for (int la=0;la<width;la++){
            for(int j=0;j<mx;j++){
                if (j>=xco && j<xco+length){
                    int r=worker->bundle->bundleState[tr].istate[la][j-xco].occupiedNum[0];
                    int g=worker->bundle->bundleState[tr].istate[la][j-xco].occupiedNum[1];
                    int value=255;
                    int temp=0;
                    if(r+g>0)  {
                        temp=255*(worker->bundle->getUnit(tr,la,j)->occupiedParticle[0].labeled);
                        image.setPixel(j+10,i+13,qRgb(value*(r>0),value*(g>0),temp));// i +3 for 3 line of white space
                       // qDebug()<<"check out of region "<<tr<<" "<<i;
                    }
                }
            }
            i=i+1;
        }
    }
  //  qDebug()<<image;
    //qDebug()<<image.pixel(2,1)<<" "<<image.pixel(3,1)<<" "<<image.pixel(4,1)<<image.pixel(1,1);
    imageLabel->setScaledContents(true);
    imageLabel->setGeometry(40-30,90,300*3+60,60*6);//(50,100,300*3,26*3);//mx*3,my*3);
    imageLabel->setPixmap(QPixmap::fromImage(image));
    // add flagmovie to record every constant time step
    extern double t_movie;
    extern bool flagmovie;
    if (flagmovie) t_movie=0;

    if (worker->recordBundle && (t>=t_movie || flagmovie)){//
        QString s;
        s.append(dirname);
        s.append("/bundle");
        s.append(QString::number(t_movie));
        s.append(".tif");
        image.scaled(imageLabel->size(), Qt::IgnoreAspectRatio).save(s,"tif");//does not support gif
        t_movie=t_movie+1;
        flagmovie=false;
    }
}
*/
void MainWindow::showPMEnds(){
    int tracknumber=worker->bundle->trackNum;
    int i=0;
    for (int tr = 0; tr<tracknumber; tr++){
        int xco=worker->bundle->getXCoordinate(tr);
        int width=worker->bundle->bundleState[tr].istate.size();
        // uncomment below to reorganize the display of bundle for eg7_39.txt // CARE of bundle size in the initial
     // /*
      /*  if (tr==0) i=0;
        else if( tr==2) i=26;//13;  for eg7_26
        else if (tr==3) i=13;
        else if (tr==1) i=13;//26;
        else if (tr==4) i=0;
        else if (tr==5) i=13;
        else if (tr==6) i=0;*/
        //else if (tr==7) i=0; // */

   /*     if (tr==0) i=0;
        else if( tr==1) i=26;//13;  for eg7_21
        else if (tr==2) i=0;
        else if (tr==3) i=13;//26;
        else if (tr==4) i=26;//13;
        else if (tr==5) i=0;
        else if (tr==6) i=13;//0;*/
        //else if (tr==7) i=0; // */
        for (int la=0;la<width;la++){
                int ps=worker->bundle->getTrack(tr)->plusEndLocations.size();
                int ms=worker->bundle->getTrack(tr)->minusEndLocations.size();
                for(int pend=0;pend<ps;pend++){
                    int site=worker->bundle->getTrack(tr)->plusEndLocations[pend];
                    if (site-xco<worker->bundle->getTrack(tr)->length){
                        //bundleImage.setPixel(site+10,i+13,qRgb(200,200,100));// NB 3 times elonged in horinzontal
                        bundleImage.setPixel(site,i,qRgb(200,200,100));
                    }
                    // include two sites for plus/minus ends marker  13_03_15
                  //  else bundleImage.setPixel(site-2,i,qRgb(250,200,100));
                    //else
                    if(site>xco){
                        //bundleImage.setPixel(site-1+10,i+13,qRgb(200,200,100));
                        bundleImage.setPixel(site-1,i,qRgb(200,200,100));
                    }
                 //   else bundleImage.setPixel(site+1,i,qRgb(250,200,100));
                }
                for(int mend=0;mend<ms;mend++){
                    int site=worker->bundle->getTrack(tr)->minusEndLocations[mend];
                    if (site-xco<worker->bundle->getTrack(tr)->length){
                        //bundleImage.setPixel(site+10,i+13,qRgb(100,200,200));
                        bundleImage.setPixel(site,i,qRgb(100,200,200));
                    }
                  //  else bundleImage.setPixel(site-2,i,qRgb(100,200,200));
                  /*  else
                    //if(site>xco){
                        bundleImage.setPixel(site+10-1,i+13,qRgb(100,200,200));*/
                   // }
                    else bundleImage.setPixel(site-1,i,qRgb(100,200,200));
                }
           i=i+1;
        }
    }
}


void MainWindow::showPMendChanged(){
    if(showPMendAct->isChecked()){
        worker->ending=true;
        showPMEnds();
        showLatticeBackup();
    }
    else{
        worker->ending=false;
        show3DEffect();
    }
}

// show just the MT view
void MainWindow::showShrink(){
    int mx=worker->bundle->maxLength;// TEST length
    int n_track=worker->bundle->trackNum;
    QImage image2(mx,n_track,QImage::Format_RGB16);//image width, image height
    image2.fill(Qt::white);
    for (int tr = 0; tr<n_track; tr++){
      for(int j = 0; j<mx;j++){
        int r=0,g=0;
        int xco=worker->bundle->getXCoordinate(tr);
        if (j<xco || j>=xco+static_cast<int>(worker->bundle->bundleState[tr].istate[0].size()))
            image2.setPixel(j,tr,qRgb(200,200,250));
        else{
            for (int la = 0; la<static_cast<int>(worker->bundle->bundleState[tr].istate.size()); la++){
                r=r+worker->bundle->bundleState[tr].istate[la][j-xco].occupiedNum[0];
                g=g+worker->bundle->bundleState[tr].istate[la][j-xco].occupiedNum[1];
            }
            r=(r>0);
            g=(g>0);
            if (r>0 || g>0) image2.setPixel(j,tr,qRgb(255*r,255*g,255.*0));
        }
      }
    }
    shrinkLabel->setScaledContents(true);
    shrinkLabel->setGeometry(50,300,300*3,10);
    shrinkLabel->setPixmap(QPixmap::fromImage(image2));

    if (worker->recordMT){
        QString s;
        s.append(dirname);
        s.append("/MT");
        s.append(QString::number(step));
        s.append(".jpg");
        image2.scaled(shrinkLabel->size(), Qt::IgnoreAspectRatio).save(s,"jpg");
    }
}



void MainWindow::show3DEffect(){
    if (threeDimensionEffctAct->isChecked()){
        worker->threeDEffect=true;
    }
    else{
        worker->threeDEffect=false;
    }
    setBundle3DBackground();
    if (worker->ending) showPMEnds();
    if (worker->blocking) showBlockages();
    showLatticeBackup();
}


void MainWindow::showTipSize(){
    int mx=worker->bundle->tipSize.size();
    if (mx>0){
    QVector<double> x(mx);
    QVector<double> counts(mx);
      for (int i=0;i<mx;i++){
        x[i]=worker->bundle->tipSize[i][0];
        counts[i]=worker->bundle->tipSize[i][1];
      }
    QVector<double> test(mx,1);
    QPen graphPen;
    graphPen.setWidthF(1);
    graphPen.setColor(Qt::red);
    plotWindow->plotState(x,counts,graphPen);

    plotWindow->customPlot->xAxis->setRange(x.front(),x.back());
    int ymax=(worker->bundle->maxLength-worker->bundle->tipRegion)*worker->bundle->trackNum*worker->bundle->laneNum;
    plotWindow->customPlot->yAxis->setRange(0,ymax);
    plotWindow->resize(621, 515);
    plotWindow->customPlot->setTitle("Tipsize Plot");
    plotWindow->customPlot->xAxis->setLabel(" steps ");
    plotWindow->customPlot->yAxis->setLabel("tip size");
    plotWindow->customPlot->replot();
    plotWindow->show();
    }
    else{
        QMessageBox msgBox;
        QString s("no data for plot");
        msgBox.setText(s);
        msgBox.exec();
    }
}

void MainWindow::startCurrentCalculate(){
    if (startCurrentAct->isChecked()){
        worker->bundle->flagCurrent=true;
        worker->initializeCurrentCalculation();
        qDebug()<<"current calculation starts";
    }
    else{
        worker->bundle->flagCurrent=false;
        worker->bundle->currentCounter=0;
    }
}


void MainWindow::startDensityCalculate(){
    if (startDensityAct->isChecked()){
        worker->bundle->flagDen=true;
        worker->initializeDensityCalculation();
    }
    else{
        worker->bundle->flagDen=false;
        worker->densityCounter=0;
    }
}

void MainWindow::startDistributionScore(){
    if (startEEDistriScoreAct->isChecked()){
        worker->bundle->flagEEDistriScore=true;
        worker->score_t.clear();
    }
    else{
        worker->bundle->flagEEDistriScore=false;
        worker->score_t.clear();
    }
}

void MainWindow::startRun()
{
    worker->running=true;
    int n=nStepEdit->text().toInt();
    while(worker->running==true){
        for (int i=0;i<n;i++){
            QString s("Running very ");
            s.append(QString::number(n)+" steps! ");
            s.append("Time elaps "+QString::number(t)+" s");
            statusBar()->showMessage(s);
            //worker->nsteps(n);
            worker->oneStep();
            if(step%50==0 & startEEDistriScoreAct->isChecked() ){//worker->bundle->flagEEDistriScore){
                  std:vector<int>state;
                  worker->projectEE(state);
                  double score=worker->calculateEvendistributionSore(state);
                 double a[2]={t,score};
                 std::vector<double>temp(a,a+2);
                worker->score_t.insert(worker->score_t.end(),temp);
            }
           /* if (worker->bundle->listOfParticles[0].type==0){
                    int temp=worker->bundle->listOfParticles[0].location[2];
                    std::cout<<temp<<std:endl;
            }*/
        }
        showImage();
        if (worker->saving==true) recordKymograph(f);
        QCoreApplication::processEvents();
    }
}


void MainWindow::startTipSize(){
    if (startTipSizeAct->isChecked()){
        worker->bundle->flagTipSize=true;
    }
    else{
        worker->bundle->flagTipSize=false;
        worker->bundle->tipSize.clear();
    }
}

void MainWindow::stopRun()
{
    worker->running=false;
    showImage();
    QCoreApplication::processEvents();
    statusBar()->showMessage("Ready");
}

void MainWindow::stopSave(){
    worker->saving=false;
    f->close();
}

void MainWindow::tipendInitialize(){


    defaultInitializerAct->setChecked(false);
    randomInitializerAct->setChecked(false);
    oneendInitializerAct->setChecked(true);
    statusBar()->showMessage("Reinitialize");
    worker->bundle->tipendInitializeState();
    worker->bundle->vkymograph.clear();
    worker->bundle->kymographRecord(kymSize,true);
    worker->bundle->initializeRateVector();
    worker->bundle->flagDen=false;
    tStep=1;
    step=0;
    t=0;
    showInitialDisplay();
    statusBar()->showMessage("Ready");
    extern bool flagmovie;
    flagmovie=true;
}

void MainWindow::test(){

}
