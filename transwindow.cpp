/*
 *this transwindow shows the urrent set of parameters used in the system and allows to change any these parameters used.
 *The parameters include physical quantities: velocities and run lengths
 * and transition rates: boundary rates, forward, lane change, turning rates
 *and inhomogeneous transition rates associated with minus/plus locations
 *
 */


#include <QMessageBox>
#include "parameter.h"
#include "transwindow.h"

TransWindow::TransWindow(QWidget *parent) :QWidget(parent){
    setWindowTitle(tr("transition rates"));
    infoLabel=new QLabel(tr("set up physical quantities to get new corresponding rates (/s)"));
    infoLabel->setFont(QFont("Helvetica", 14, QFont::Bold));
    infoLabel->setFrameStyle(QFrame::Panel | QFrame::Sunken);
    infoLabel->setAlignment(Qt::AlignCenter);
    QHBoxLayout *infoLayout=new QHBoxLayout;

    setLabel=new QLabel(tr("Click the buttons on the left to get new corresponding transition rates and display here the explainations!"));

    QString s1=QString("");//relationship between run length, velocity, forward rates, turning rates and lane changes rates:\n");
    s1.append("1. run length equals velocity over corresponding turning rate while velocity is forwad rate (including unblocked lane change rates) times the step size.\n");
    s1.append("2. default forward rate for type 0 equals velocity over step size and unblocked lane change rate 0.\n");
    s1.append("3. default unblocked lane change rate for type 1 is velocity*delta/2 (delta=");
    s1.append(QString::number(delta));
    s1.append(" times/um) and foward rate=velocity/stepsize -2*unblocked lane change rate.\n");

    s1.append("4. turning rate is velocity over corresponding runlength \n");
    s1.append("5. bundle region (unit of site): MTs coulbe be bundled to enable inhomogeneous turning or switching\n");
    s1.append("NB: lane change when blocked will be updated automatically to be half of the overall forward rate as pushing left buttons.");

    setLabel->setText(s1);
    //setEdit->setEnabled(false);

    boundaryBox=new QGroupBox(tr("boundary rate"));
    dataBox=new QGroupBox(tr("physical data"));
    singleRateBox=new QGroupBox(tr("single track rate grid"));
    inhomoBox=new QGroupBox(tr("inhomogeneous rate grid"));
    switchBox=new QGroupBox(tr("bundle-associate rate grid (continue values indicate the multiplication factor of single track rate)"));

    QGridLayout *velocityGrid=new QGridLayout;// for velocity and runlength parameters in physical units
    QGridLayout *boundaryGrid=new QGridLayout;
    QGridLayout *setGrid=new QGridLayout;// set forward and turning rates
    QGridLayout *paraGrid=new QGridLayout;// for turning , forward and lane change rates
   // QGridLayout *otherGrid=new QGridLayout;// other insinglegeneous rates particularly in multiple tracks
    QHBoxLayout *checkLayout=new QHBoxLayout;// for confirm, cancel buttons
    QGridLayout *switchGrid=new QGridLayout;//track associate
    QGridLayout *inhomoGrid=new QGridLayout;

    // bottons
    QPushButton *bgetForward=new QPushButton("get forward and turning rates",this);
    connect(bgetForward,SIGNAL(clicked()),this,SLOT(getForward()));
    //QPushButton *bgetTurn=new QPushButton("set turning rates",this);
    //connect(bgetTurn,SIGNAL(clicked()),this,SLOT(getTurning()));
    QPushButton *bgetVeloRunLen=new QPushButton("get velocities and run length",this);
    connect(bgetVeloRunLen,SIGNAL(clicked()),this,SLOT(getRunLength()));
    //QPushButton *getRunLength=new QPushButton("set run length",this);
    QPushButton *bconfirm=new QPushButton("OK", this);
    connect(bconfirm,SIGNAL(clicked()),this,SLOT(checkRates()));
    QPushButton *bcancel=new QPushButton("Cancel", this);
    connect(bcancel,SIGNAL(clicked()),this,SLOT(close()));
    QPushButton *bdefault=new QPushButton("Default",this);
    connect(bdefault,SIGNAL(clicked()),this,SLOT(setDefault()));
    //connect(bdefault,SIGNAL(clicked()),this,SLOT(close()));
    inhomoPlusCheck=new QCheckBox("load rate",this);

    connect(inhomoPlusCheck,SIGNAL(toggled(bool)),this,SLOT(inhomoPlusChanged()));

    checkLayout->addWidget(bconfirm);
    checkLayout->addWidget(bcancel);
    checkLayout->addWidget(bdefault);

    int ntype=num_of_type;
    QVector<QLabel *>typeLabel_1(ntype);//velocitygrid
    QVector<QLabel *>typeLabel_2(ntype);//for paraGrid
    QVector<QLabel *>typeLabel_4(ntype);//for boundaryGrid
    QVector<QLabel *>typeLabel_3(ntype*2);// for infLayout
    QVector<QLabel *>typeLabel_5(ntype);// for track switch grid

   // infLayout
    for (int i=0;i<4;++i){
        typeLabel_3[i]=new QLabel();
        infoLayout->addWidget(typeLabel_3[i]);
    }
    typeLabel_3[0]->setText("type 0: plus-directed particles");
    QImage image(QSize(10,10), QImage::Format_RGB32);
    image.fill(Qt::red);
    typeLabel_3[1]->setPixmap(QPixmap::fromImage(image));
    typeLabel_3[2]->setText("type 1: minus-directed particles");
    image.fill(Qt::green);
    typeLabel_3[3]->setPixmap(QPixmap::fromImage(image));



    injectEdit.resize(ntype);
    exitEdit.resize(ntype);
    velocityEdit.resize(ntype);
    runBundledEdit.resize(ntype);
    runUnbundledEdit.resize(ntype);
    turnuniEdit.resize(ntype);
    turnbiEdit.resize(ntype);
    turninhomoEdit.resize(3);
    forwardEdit.resize(ntype);
    laneChangeEdit.resize(ntype*2);
    forwardEndEdit.resize(4);
    switchInnerEdit.resize(4);
    switchEndEdit.resize(12);

    QLabel *stepSizeLabel=new QLabel("step size");

    QLabel *injectLabel=new QLabel("injection");
    QLabel *exitLabel=new QLabel("exit");
    QLabel *velocityLabel=new QLabel("velocities (um/s)");
    QLabel *runlengthunipolarLabel=new QLabel("run length(unbundled, um)");
    QLabel *runlengthbipolarLabel=new QLabel("run length(bundled, um)");
    QLabel *turnunbundledLabel=new QLabel("turning rates (unbundled)");
    QLabel *turnbundledLabel=new QLabel("turning rates (bundled)");
    QLabel *turninhomoLabel1=new QLabel("at its (-) end");
    QLabel *turninhomoLabel2=new QLabel("at other (-) ends");
    QLabel *turnloadLabel=new QLabel("at other (+) ends");
    QLabel *crossLabel=new QLabel("cross (-)");
    QLabel *bundleRegionLabel=new QLabel("bundle region:");
    QLabel *bundleRegiontoLabel = new QLabel(" - ");
    QLabel *forwardLabel=new QLabel("forward rates");
    QLabel *laneChangeUnblockedLable=new QLabel("lane change rates (unblocked)");
    QLabel *laneChangeBlockedLable=new QLabel("lane change rates (blocked)");
    QLabel *switchWithinBundleLabel=new QLabel("within bundle");
    QLabel *switchUnipolarJunctionLabel=new QLabel("unipolar bundle junction");
    QLabel *switchBipolarJunctionLabel=new QLabel("bipolar bundle junction");
    QVector<QLabel *>switchTypeChangeLabel(3);
    QVector<QLabel *>switchNoTypeChangeLabel(3);
    QVector<QLabel *>switchLabel(2);
    QVector<QLabel *>switchContinueLabel(2);
    QVector<QLabel *>switchIntoEndLabel(2);
    QVector<QLabel *>switchLeaveEndLabel(2);
    for (int i=0;i<2;i++){
        switchTypeChangeLabel[i]=new QLabel("switch & type change");
        switchNoTypeChangeLabel[i]=new QLabel("switch & no type change");
        switchLabel[i]=new QLabel("switch");
        switchContinueLabel[i]=new QLabel("continue");
        switchIntoEndLabel[i]=new QLabel("into ends");
        switchLeaveEndLabel[i]=new QLabel("leave ends");
    }
    switchTypeChangeLabel[2]=new QLabel("switch & type change");
    switchNoTypeChangeLabel[2]=new QLabel("switch & no type change");

    velocityGrid->addWidget(velocityLabel,1,3);
    velocityGrid->addWidget(runlengthunipolarLabel,1,4);
    velocityGrid->addWidget(runlengthbipolarLabel,1,5);


    boundaryGrid->addWidget(injectLabel,1,2);
    boundaryGrid->addWidget(exitLabel,1,3);
    boundaryBox->setLayout(boundaryGrid);

    stepSizeEdit=new QLineEdit("");
    stepSizeEdit->setText(QString::number(hs));
    stepSizeEdit->setEnabled(true);
    bundleRegionEdit.resize(2);
    for (int i=0;i<2;i++){
    bundleRegionEdit[i]=new QLineEdit("");
    bundleRegionEdit[i]->setEnabled(true);
    }
    bundleRegionEdit[0]->setText(QString::number(B1));
    bundleRegionEdit[1]->setText(QString::number(B2));
    paraGrid->addWidget(bgetForward,2,1,1,1);
    velocityGrid->addWidget(bgetVeloRunLen,2,1,1,1);

    setGrid->addWidget(stepSizeLabel,1,1);
    setGrid->addWidget(stepSizeEdit,1,3);
    setGrid->addWidget(bundleRegionLabel,2,1);
    setGrid->addWidget(bundleRegionEdit[0],3,1);
    setGrid->addWidget(bundleRegiontoLabel,3,2);
    setGrid->addWidget(bundleRegionEdit[1],3,3);

    setGrid->addWidget(setLabel,1,4,3,1);
    setGrid->addWidget(boundaryBox,1,5,3,3);
   // setGrid->setColumnStretch(2, 1);
   // setGrid->setColumnStretch(3, 4);


    paraGrid->addWidget(turnunbundledLabel,1,4);
    paraGrid->addWidget(turnbundledLabel,1,5);
    paraGrid->addWidget(forwardLabel,1,6);
    paraGrid->addWidget(laneChangeUnblockedLable,1,7);
    paraGrid->addWidget(laneChangeBlockedLable,1,8);

    paraGrid->addWidget(inhomoBox,5,1,1,8);


    ploadEdit=new QLineEdit("");
    turninhomoEdit[0]=new QLineEdit("");//for plus end of its own MTs
    turninhomoEdit[1]=new QLineEdit("");//for plus end of other MTs
    turninhomoEdit[2]=new QLineEdit("");// for plus end of other MTs

    QString s=QString("");
    s.append("1. inhomogeneous turning rate for type-1 particles at a - end of its own (left) or of another track with single minus in opposite orientation (right)\n");
    s.append("2. inhomogeneous rate at plus region is for type 0 particles and the rate is proportion to the density\n");
    s.append("which is assumed to be linear depending on the region size and dynein comet size.\n");
    s.append("3. tick box to enable the inhomogeneous rate for type-0 particles at (sub)apticals\n");
    s.append("4. inhomogeneous turning rate for type-0 particles when meeting subapical + ends of other MTs in the same orientation\n");
    s.append("5. cross (-): the rate of type-1 particles crossing the double minus end on the its track.\n");

    turnInfLabel=new QLabel(tr(""));
    turnInfLabel->setText(s);// setPlainText(s);
    crossEdit=new QLineEdit("");
    inhomoBox->setLayout(inhomoGrid);
    inhomoGrid->addWidget(turninhomoLabel1,1,1);
    inhomoGrid->addWidget(inhomoPlusCheck,2,1,1,1);
    inhomoGrid->addWidget(turninhomoEdit[0],1,2);//inhomogeneous turning at minus end for dynein
    inhomoGrid->addWidget(turninhomoLabel2,1,3);//inhomogeneous turning at minus end for dynein
    inhomoGrid->addWidget(turninhomoEdit[1],1,4);//inhomogeneous turning at minus end for dynein
    inhomoGrid->addWidget(ploadEdit,2,2);
    inhomoGrid->addWidget(turnloadLabel,2,3);
    inhomoGrid->addWidget(turninhomoEdit[2],2,4);//inhomogeneous turning at minus end for dynein
    inhomoGrid->addWidget(turnInfLabel,1,5,3,3);
    inhomoGrid->addWidget(crossLabel,3,1,1,2);//inhomogeneous turning at minus end for dynein
    inhomoGrid->addWidget(crossEdit,3,3,1,2);//inhomogeneous turning at minus end for dynein
  //  inhomoGrid->setColumnStretch(1, 1);
    //inhomoGrid->setColumnStretch(4, 5);

    for (int i=0;i<ntype;i++){
        injectEdit[i]=new QLineEdit("");
        exitEdit[i]=new QLineEdit("");
        velocityEdit[i]=new QLineEdit("");
        runBundledEdit[i]=new QLineEdit("");
        runUnbundledEdit[i]=new QLineEdit("");
        turnuniEdit[i]=new QLineEdit("");
        turnbiEdit[i]=new QLineEdit("");
        forwardEdit[i]=new QLineEdit("");
        laneChangeEdit[i*2]=new QLineEdit("");
        laneChangeEdit[i*2+1]=new QLineEdit("");
        QString s("type ");
        s.append(QString::number(i));
        typeLabel_1[i]=new QLabel(s);
        typeLabel_2[i]=new QLabel(s);
        typeLabel_4[i]=new QLabel(s);
        typeLabel_5[i]=new QLabel(s);
        boundaryGrid->addWidget(typeLabel_4[i],i+2,1);
        boundaryGrid->addWidget(injectEdit[i],i+2,2);
        boundaryGrid->addWidget(exitEdit[i],i+2,3);
        velocityGrid->addWidget(typeLabel_1[i],i+2,2);
        velocityGrid->addWidget(velocityEdit[i],i+2,3);
        velocityGrid->addWidget(runUnbundledEdit[i],i+2,4);
        velocityGrid->addWidget(runBundledEdit[i],i+2,5);
        paraGrid->addWidget(typeLabel_2[i],i+2,3);
        paraGrid->addWidget(turnuniEdit[i],i+2,4);
        paraGrid->addWidget(turnbiEdit[i],i+2,5);
        paraGrid->addWidget(forwardEdit[i],i+2,6);
        paraGrid->addWidget(laneChangeEdit[i*2],i+2,7);//unblocked
        paraGrid->addWidget(laneChangeEdit[i*2+1],i+2,8);//blocked
    }
   for (int i=0;i<forwardEndEdit.size();i++) forwardEndEdit[i]=new QLineEdit("");
    for (int i=0;i<switchInnerEdit.size();i++) switchInnerEdit[i]=new QLineEdit("");
    for (int i=0;i<switchEndEdit.size();i++) switchEndEdit[i]=new QLineEdit("");




    switchGrid->addWidget(switchWithinBundleLabel,1,2,1,2,Qt::AlignCenter);
    switchGrid->addWidget(switchUnipolarJunctionLabel,1,4,1,4,Qt::AlignCenter);
    switchGrid->addWidget(switchBipolarJunctionLabel,1,8,1,4,Qt::AlignCenter);
    switchGrid->addWidget(switchTypeChangeLabel[0],2,2,2,1,Qt::AlignCenter);
    switchGrid->addWidget(switchNoTypeChangeLabel[0],2,3,2,1,Qt::AlignCenter);
    for (int i=0;i<2;i++){
        switchGrid->addWidget(switchIntoEndLabel[i],2,i*4+4,1,2,Qt::AlignCenter);
        switchGrid->addWidget(switchLeaveEndLabel[i],2,i*4+6,1,2,Qt::AlignCenter);
        switchGrid->addWidget(switchContinueLabel[i],3,i*4+4);
        switchGrid->addWidget(switchLabel[i],3,i*4+5);
        switchGrid->addWidget(switchTypeChangeLabel[i+1],3,i*4+6);
        switchGrid->addWidget(switchNoTypeChangeLabel[i+1],3,i*4+7);
        switchGrid->addWidget(switchInnerEdit[i*2],4+i,2);
        switchGrid->addWidget(switchInnerEdit[i*2+1],4+i,3);
        switchGrid->addWidget(forwardEndEdit[i*2],4+i,4);
        switchGrid->addWidget(forwardEndEdit[i*2+1],4+i,8);
        for (int j=0;j<6;j++)
        switchGrid->addWidget(switchEndEdit[i*6+j],4+i,5*(j<3)+6*(j>=3)+j);
    }

    dataBox->setLayout(velocityGrid);
    singleRateBox->setLayout(paraGrid);
    //bundleRateBox->setLayout(otherGrid);
    switchBox->setLayout(switchGrid);
    switchGrid->addWidget(typeLabel_5[0],4,1);
    switchGrid->addWidget(typeLabel_5[1],5,1);
    QVBoxLayout *mainLayout=new QVBoxLayout;
    mainLayout->addWidget(infoLabel);
    mainLayout->addLayout(infoLayout);
    mainLayout->addWidget(dataBox);
    mainLayout->addLayout(setGrid);
    mainLayout->addWidget(singleRateBox);
   // mainLayout->addWidget(bundleRateBox);
    mainLayout->addWidget(switchBox);
    mainLayout->addLayout(checkLayout);
    setLayout(mainLayout);

}


void TransWindow::bindWorker(Worker *wk) {
    worker = wk;
}
void TransWindow::bundleLengthMessage(){
    QMessageBox messageBox(this);
    QString s("Current physical bundle length is ");
    double len=double(worker->bundle->maxLength)*stepSizeEdit->text().toDouble();
    s.append(QString::number(len)+"um!\n");
    s.append("you may wish to change physical bundle length by chaning number of site via changing tracks is the Geometry menu!");
    messageBox.setText(s);
    messageBox.exec();
}

void TransWindow::getNewRates(){
    worker->bundle->stepSize=stepSizeEdit->text().toDouble();
    worker->bundle->bundleRegion_1=bundleRegionEdit[0]->text().toDouble();
    worker->bundle->bundleRegion_2=bundleRegionEdit[1]->text().toDouble();

    for (int i=0;i<velocityEdit.size();i++){
        worker->bundle->velocity[i]=velocityEdit[i]->text().toDouble();
        worker->bundle->runLengthBundled[i]=runBundledEdit[i]->text().toDouble();
        worker->bundle->runLengthUnbundled[i]=runUnbundledEdit[i]->text().toDouble();
        worker->bundle->transitionRates.forward[i]=forwardEdit[i]->text().toDouble();
        worker->bundle->transitionRates.laneChange[i*2]=laneChangeEdit[i*2]->text().toDouble();
        worker->bundle->transitionRates.laneChange[i*2+1]=laneChangeEdit[i*2+1]->text().toDouble();
        worker->bundle->transitionRates.turning[i*2]=turnuniEdit[i]->text().toDouble();
        worker->bundle->transitionRates.turning[i*2+1]=turnbiEdit[i]->text().toDouble();
       for (int tr=0;tr<worker->bundle->trackNum;tr++){
           worker->bundle->transitionRates.boundaryIn[i][tr]=injectEdit[i]->text().toDouble();
           worker->bundle->transitionRates.boundaryOut[i][tr]=exitEdit[i]->text().toDouble();
       }
    }

    worker->bundle->transitionRates.turning[4]=turninhomoEdit[0]->text().toDouble();
    worker->bundle->transitionRates.turning[5]=turninhomoEdit[1]->text().toDouble();
    worker->bundle->transitionRates.turning[6]=turninhomoEdit[2]->text().toDouble();
    for (int i=0;i<forwardEndEdit.size();i++){
         worker->bundle->transitionRates.forwardEnd[i]=forwardEndEdit[i]->text().toDouble();
    }

    for (int i=0;i<switchEndEdit.size();i++){
         worker->bundle->transitionRates.trackSwitchAtEnd[i]=switchEndEdit[i]->text().toDouble();
    }
    for (int i=0;i<switchInnerEdit.size();i++){
         worker->bundle->transitionRates.trackSwitchAtInner[i]=switchInnerEdit[i]->text().toDouble();
    }
    if (inhomoPlusCheck->isChecked()) {
        //ploadEdit->setEnabled(true);
        worker->bundle->flagInhomo=true;
    }
    else {
        //ploadEdit->setEnabled(false);
        worker->bundle->flagInhomo=false;
    }

    worker->bundle->pLoad=ploadEdit->text().toDouble();
    worker->bundle->transitionRates.forward[2]=crossEdit->text().toDouble();
    //worker->bundle->minRate();
    worker->bundle->initializeRateVector();
}

void TransWindow::getForward(){

    double v=velocityEdit[0]->text().toDouble();
    double h=stepSizeEdit->text().toDouble();
    double newRate=v/h;
    forwardEdit[0]->setText(QString::number(newRate));
    laneChangeEdit[0]->setText("0");
    laneChangeEdit[1]->setText(QString::number(newRate/2.0));
    v=velocityEdit[1]->text().toDouble();
    double dance=delta;//times per nm;
    if (v/h-v*dance>=0){
        laneChangeEdit[2]->setText(QString::number(v*dance/2.0));  //lane_d=v*delta/2;//delta*v
        forwardEdit[1]->setText(QString::number(v/h-v*dance));//v/hs-v*delta
        laneChangeEdit[3]->setText(QString::number(v/h/2.0));
    }
    else{
        laneChangeEdit[2]->setText(QString::number(v/h/2.0));  //lane_d=v*delta/2;//delta*v
        forwardEdit[1]->setText("0");
        laneChangeEdit[3]->setText(QString::number(v/h/2.0));
    }
    getTurning();
}


void TransWindow::getTurning(){    
   // setEdit->setText("turning rate is velocity over corresponding runlength \n");
    for (int i=0;i<turnbiEdit.size();++i){
        turnbiEdit[i]->setText(QString::number(velocityEdit[i]->text().toDouble()/runBundledEdit[i]->text().toDouble()));
        turnuniEdit[i]->setText(QString::number(velocityEdit[i]->text().toDouble()/runUnbundledEdit[i]->text().toDouble()));
    }
}

void TransWindow::getRunLength(){

    double h=stepSizeEdit->text().toDouble();
    getVelocity();
    for (int i=0;i<turnbiEdit.size();++i){
        runBundledEdit[i]->setText(QString::number(velocityEdit[i]->text().toDouble()/turnbiEdit[i]->text().toDouble()));
        runUnbundledEdit[i]->setText(QString::number(velocityEdit[i]->text().toDouble()/turnuniEdit[i]->text().toDouble()));
        laneChangeEdit[i*2+1]->setText(QString::number(velocityEdit[i]->text().toDouble()/h/2.0));
    }
}

void TransWindow::getVelocity(){
    double h=stepSizeEdit->text().toDouble();
    double p=laneChangeEdit[0]->text().toDouble()*2.0+forwardEdit[0]->text().toDouble();
    double newRate=p*h;
    velocityEdit[0]->setText(QString::number(newRate));
    p=laneChangeEdit[2]->text().toDouble()*2.0+forwardEdit[1]->text().toDouble();
    newRate=p*h;
    velocityEdit[1]->setText(QString::number(newRate));
}

void TransWindow::setOldRates(){
    stepSizeEdit->setText(QString::number(worker->bundle->stepSize));
    bundleRegionEdit[0]->setText(QString::number(worker->bundle->bundleRegion_1));
    bundleRegionEdit[1]->setText(QString::number(worker->bundle->bundleRegion_2));
    for (int i=0;i<velocityEdit.size();i++){
        for (int tr=0;tr<worker->bundle->trackNum;tr++){
            injectEdit[i]->setText(QString::number(worker->bundle->transitionRates.boundaryIn[i][tr]));
            exitEdit[i]->setText(QString::number(worker->bundle->transitionRates.boundaryOut[i][tr]));
        }
        velocityEdit[i]->setText(QString::number(worker->bundle->velocity[i]));
        runBundledEdit[i]->setText(QString::number(worker->bundle->runLengthBundled[i]));
        runUnbundledEdit[i]->setText(QString::number(worker->bundle->runLengthUnbundled[i]));
        turnuniEdit[i]->setText(QString::number(worker->bundle->transitionRates.turning[i*2]));
        turnbiEdit[i]->setText(QString::number(worker->bundle->transitionRates.turning[i*2+1]));
        forwardEdit[i]->setText(QString::number(worker->bundle->transitionRates.forward[i]));
        laneChangeEdit[i*2]->setText(QString::number(worker->bundle->transitionRates.laneChange[i*2]));
        laneChangeEdit[i*2+1]->setText(QString::number(worker->bundle->transitionRates.laneChange[i*2+1]));
    }
    turninhomoEdit[0]->setText(QString::number(worker->bundle->transitionRates.turning[4]));
    turninhomoEdit[1]->setText(QString::number(worker->bundle->transitionRates.turning[5]));
    turninhomoEdit[2]->setText(QString::number(worker->bundle->transitionRates.turning[6]));
    for (int i=0;i<forwardEndEdit.size();i++){
         forwardEndEdit[i]->setText(QString::number(worker->bundle->transitionRates.forwardEnd[i]));
    }
    for (int i=0;i<switchEndEdit.size();i++){
         switchEndEdit[i]->setText(QString::number(worker->bundle->transitionRates.trackSwitchAtEnd[i]));
    }
    for (int i=0;i<switchInnerEdit.size();i++){
         switchInnerEdit[i]->setText(QString::number(worker->bundle->transitionRates.trackSwitchAtInner[i]));
    }
    ploadEdit->setText(QString::number(worker->bundle->pLoad));
    if (worker->bundle->flagInhomo) ploadEdit->setEnabled(true);
    else ploadEdit->setEnabled(false);

    if (worker->bundle->flagInhomo) inhomoPlusCheck->setChecked(true);
    else  inhomoPlusCheck->setChecked(false);

    crossEdit->setText(QString::number(worker->bundle->transitionRates.forward[2]));
    if (worker->bundle->trackNum<=1){

    }
    else{

    }

}

void TransWindow::setDefault(){
    worker->bundle->initializeTransitionRates();
    setOldRates();
    //setEdit->setText("Click the buttons on the left to get new corresponding transition rates and display here the explainations!");
}

void TransWindow::checkRates(){
    QMessageBox messageBox(this);
    QString s("Please check the inherit relation bewteen physical quantities and forward and turning rates is satisfied by clicking the buttons when you change the physical quantities!");
    messageBox.setText(s);
    messageBox.setInformativeText("Do you want to use the current set of rates?");
    QAbstractButton *confirmButton=messageBox.addButton(tr("Yes"), QMessageBox::AcceptRole);
    QAbstractButton *cancelButton=messageBox.addButton(tr("Reset"),QMessageBox::RejectRole);
    messageBox.exec();
    if (messageBox.clickedButton() == confirmButton) {
        getNewRates();
        close();
    }
    else if(messageBox.clickedButton()==cancelButton){
        close();
    }
}


void TransWindow::closeEvent(QCloseEvent *){
   bundleLengthMessage();
}

void TransWindow::inhomoPlusChanged(){
    if (inhomoPlusCheck->isChecked()){
        worker->bundle->flagInhomo=true;
        ploadEdit->setEnabled(true);
    }
    else{
        worker->bundle->flagInhomo=false;
        ploadEdit->setEnabled(false);
    }
}
