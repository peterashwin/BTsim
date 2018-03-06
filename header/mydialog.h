#ifndef MYDIALOG_H
#define MYDIALOG_H

#include <QDialog>
#include <QSpinBox>
#include <QLineEdit>
#include <QPushButton>
#include "parameter.h"
#include "worker.h"

namespace Ui {
class MyDialog;
}

class MyDialog : public QDialog
{
    Q_OBJECT

public:
    explicit MyDialog(QWidget *parent = 0);
    ~MyDialog();
    void setPara(DialogPara &, bool &);
    void bindWorkerDialog(Worker *wk);

    DialogPara para;
protected:
    Worker *workerDialog;
    bool okFlag;

signals:


private slots:
     void cancel();
     void confirm();
     void getPara();// get track parameters include xcoordinate, minus/plus locations, width and track index
private:
     QPushButton *confirmButton;
     QPushButton *closeButton;

     QLineEdit *indexEdit;
     QLineEdit *widthEdit;
     QLineEdit *minusLocationEdit;
     QLineEdit *plusLocationEdit;
     QLineEdit *distanceEdit;
     QLineEdit *dyneinEdit;
public:
     void getMinus();//get minus locations from minusLocationEdit
     void getPlus();//get plus locations from plusLocationEdit
     void getDistance();
     void getDynein();
private:
    Ui::MyDialog *ui;
};

#endif // MYDIALOG_H
