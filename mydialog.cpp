/*
 *this mydialog is used to add and change tracks in the bundle
 */
#include "MyDialog.h"
#include "ui_mydialog.h"
#include <QtGui>
#include <QLabel>
#include <QBoxLayout>

MyDialog::MyDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::MyDialog)
{
    ui->setupUi(this);
    okFlag=false;

    QLabel *trackLabel=new QLabel(tr("track index"));
    QLabel *trackwidthLabel = new QLabel(tr("width"));
    QLabel *trackMinusLabel = new QLabel(tr("minus locations(;)"));
    QLabel *trackPlusLabel = new QLabel(tr("plus locations(;)"));
    QLabel *trackDisLabel = new QLabel(tr("distance of dynein comet"));
    QLabel *trackDyneinLabel = new QLabel(tr("# dynein"));

    indexEdit = new QLineEdit("0");
    widthEdit =new QLineEdit(QString::number(num_of_lane));
    minusLocationEdit = new QLineEdit(QString::number(0));
    plusLocationEdit = new QLineEdit(QString::number(L));
    distanceEdit = new QLineEdit(QString::number(int(distancePlus/hs)));
    dyneinEdit = new QLineEdit(QString::number(num_of_dynein));
    confirmButton = new QPushButton(tr("confirm"));
    confirmButton->setDefault(true);
    confirmButton->setEnabled(true);
    closeButton = new QPushButton(tr("Cancel"));

    connect(confirmButton, SIGNAL(clicked()),this, SLOT(getPara()));
    connect(confirmButton,SIGNAL(clicked()),this,SLOT(confirm()));
    connect(confirmButton,SIGNAL(clicked()),this,SLOT(close()));
    connect(closeButton, SIGNAL(clicked()),this, SLOT(cancel()));
    connect(closeButton,SIGNAL(clicked()),this,SLOT(close()));
    QGridLayout *topLayout = new QGridLayout;
    topLayout->addWidget(trackLabel,1,1);
    topLayout->addWidget(trackwidthLabel,1,2);
    topLayout->addWidget(trackMinusLabel,1,3);
    topLayout->addWidget(trackPlusLabel,1,4);
    topLayout->addWidget(trackDisLabel,1,5);
    topLayout->addWidget(trackDyneinLabel,1,6);

    //QHBoxLayout *top1Layout = new QHBoxLayout;
    topLayout->addWidget(indexEdit,2,1);
    topLayout->addWidget(widthEdit,2,2);
    topLayout->addWidget(minusLocationEdit,2,3);
    topLayout->addWidget(plusLocationEdit,2,4);
    topLayout->addWidget(distanceEdit,2,5);
    topLayout->addWidget(dyneinEdit,2,6);


    QHBoxLayout *buttomLayout = new QHBoxLayout;
    buttomLayout->addWidget(confirmButton);
    buttomLayout->addWidget(closeButton);
    buttomLayout->addStretch();

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addLayout(topLayout);
    //mainLayout->addLayout(top1Layout);
    mainLayout->addLayout(buttomLayout);
    setLayout(mainLayout);
    setWindowTitle(tr("bundle geometry"));
    setFixedHeight(sizeHint().height());
}

MyDialog::~MyDialog()
{
    delete ui;
}

void MyDialog::cancel()
 {
     okFlag=false;
 }

void MyDialog::confirm(){
    okFlag=true;
}

void MyDialog::setPara(DialogPara & initial, bool & ok){
    ok=okFlag;
    qDebug()<<"set Parameter "<<ok;
    if (ok==true) initial=para;
}


void MyDialog::bindWorkerDialog(Worker *wk)
{
    workerDialog=wk;
}

void MyDialog::getPara(){
    para.trackindex=indexEdit->text().toInt();
    para.width=widthEdit->text().toInt();
    getMinus();
    getPlus();
    getDistance();
    getDynein();
}


void MyDialog::getMinus(){
    para.minus.clear();
    QString s=minusLocationEdit->text();
    while (s.size()>0){
        if (!s.contains(";")){
        para.minus.push_back(s.toInt());
        s="";
        }
        else{
           int m=s.indexOf(";",0);
           QString temp(s.left(m));
           para.minus.push_back(temp.toInt());
           s.remove(0,m+1);
        }
    }
}

void MyDialog::getPlus(){
    para.plus.clear();
    QString s=plusLocationEdit->text();
    while (s.size()>0){
        if (!s.contains(";")){
        para.plus.push_back(s.toInt());
        s="";
        }
        else{
           int m=s.indexOf(";",0);
           QString temp(s.left(m));
           para.plus.push_back(temp.toInt());
           s.remove(0,m+1);
        }
    }
}

void MyDialog::getDistance(){
    para.distance.clear();
    QString s=distanceEdit->text();
    while (s.size()>0){
        if (!s.contains(";")){
        para.distance.push_back(s.toInt());
        s="";
        }
        else{
           int m=s.indexOf(";",0);
           QString temp(s.left(m));
           para.distance.push_back(temp.toInt());
           s.remove(0,m+1);
        }
    }
}

void MyDialog::getDynein(){
    para.dynein.clear();
    QString s=dyneinEdit->text();
    while (s.size()>0){
        if (!s.contains(";")){
        para.dynein.push_back(s.toInt());
        s="";
        }
        else{
           int m=s.indexOf(";",0);
           QString temp(s.left(m));
           para.dynein.push_back(temp.toInt());
           s.remove(0,m+1);
        }
    }
}
