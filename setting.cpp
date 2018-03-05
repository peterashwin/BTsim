/*
 *this widget display the current geometry of the bundle set up in the worker object;
 *information includes number of tracks, plus\minus locations of each track
 */

#include <QLabel>
#include <QVector>
#include <QBoxLayout>
#include <QGridLayout>
#include "setting.h"
#include "worker.h"
#include <QDebug>
Setting::Setting(QWidget *parent) :
    QWidget(parent)
{
    editLayout=new QVBoxLayout;

    QLabel *trackLabel=new QLabel(tr("index"));
    QLabel *trackWidthLabel = new QLabel(tr("track width"));
    QLabel *minusLabel = new QLabel(tr("minus location"));
    QLabel *plusLabel = new QLabel(tr("plus locations"));
    QLabel *disLabel=new QLabel(tr("length of dynein comet"));
    QLabel *dyneinLabel = new QLabel(tr("# Dynein at plus ends"));

    topLeftLayout = new QHBoxLayout;
    topLeftLayout->addWidget(trackLabel);
    topLeftLayout->addWidget(trackWidthLabel);
    topLeftLayout->addWidget(minusLabel);
    topLeftLayout->addWidget(plusLabel);
    topLeftLayout->addWidget(disLabel);
    topLeftLayout->addWidget(dyneinLabel);



    mainLayout = new QGridLayout;
    mainLayout->addWidget(trackLabel,1,1);
    mainLayout->addWidget(trackWidthLabel,1,2);
    mainLayout->addWidget(minusLabel,1,3);
    mainLayout->addWidget(plusLabel,1,4);
    mainLayout->addWidget(disLabel,1,5);
    mainLayout->addWidget(dyneinLabel,1,6);
    for (int i=1;i<10;i++){
        QLabel *trackLabel=new QLabel("",this);
        QLineEdit *widthEdit = new QLineEdit("",this);
        QLineEdit *minusEdit = new QLineEdit("",this);
        QLineEdit *plusEdit = new QLineEdit("",this);
        QLineEdit *disEdit=new QLineEdit("",this);
        QLineEdit *dyneinEdit=new QLineEdit("",this);
        widthEdit->setEnabled(false);
        minusEdit->setEnabled(false);
        plusEdit->setEnabled(false);
        disEdit->setEnabled(false);
        dyneinEdit->setEnabled(false);
        vminusEdit.push_back(minusEdit);
        vplusEdit.push_back(plusEdit);
        vwidthEdit.push_back(widthEdit);
        vlabelEdit.push_back(trackLabel);
        vdisEdit.push_back(disEdit);
        vdyneinEdit.push_back(dyneinEdit);
        mainLayout->addWidget(trackLabel,i+1,1);
        mainLayout->addWidget(widthEdit,i+1,2);
        mainLayout->addWidget(minusEdit,i+1,3);
        mainLayout->addWidget(plusEdit,i+1,4);
        mainLayout->addWidget(disEdit,i+1,5);
        mainLayout->addWidget(dyneinEdit,i+1,6);
    }

    setLayout(mainLayout);
    setWindowTitle(tr("bundle geometry"));
}

void Setting::bindWorker(Worker *wk) {
    worker = wk;
}



void Setting::setUp(){
    int trackN=worker->bundle->trackNum;
    for (int i=0;i<trackN;i++){
        QString tstr("track ");
        tstr.append(QString::number(i));
        int width=worker->bundle->bundleState[i].istate.size();
        std::vector<int> minus=worker->bundle->bundleState[i].minusEndLocations;
        std::vector<int> plus=worker->bundle->bundleState[i].plusEndLocations;
        std::vector<int> dis=worker->bundle->bundleState[i].inhomoDistances;
        std::vector<int> dynein=worker->bundle->bundleState[i].dyneinNum;


        QString mstr,pstr,wstr,distr,dystr;
        wstr.append(QString::number(width));
        mstr.append(QString::number(minus[0]));
        pstr.append(QString::number(plus[0]));
        distr.append(QString::number(dis[0]));
        dystr.append(QString::number(dynein[0]));
        //qDebug()<<dystr;
        if (minus.size()>1){
            for (std::vector<int>::size_type j=1;j<minus.size();j++){
                mstr.append(";");
                mstr.append(QString::number(minus[j]));
            }
        }
        if (plus.size()>1){
            for (std::vector<int>::size_type j=1;j<plus.size();j++){
                pstr.append(";");
                pstr.append(QString::number(plus[j]));
                distr.append(";");
                distr.append(QString::number(dis[j]));
                dystr.append(";");
                dystr.append(QString::number(dynein[j]));
            }
        }
        vlabelEdit[i]->setText(tstr);
        vwidthEdit[i]->setText(wstr);
        vminusEdit[i]->setText(mstr);
        vplusEdit[i]->setText(pstr);
        vdisEdit[i]->setText(distr);
        vdyneinEdit[i]->setText(dystr);
    }

}
