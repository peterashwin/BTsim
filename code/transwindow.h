#ifndef TRANSWINDOW_H
#define TRANSWINDOW_H

#include <QPushButton>
#include <QLineEdit>
#include <QSpinBox>
#include <QVBoxLayout>
#include <QLabel>
#include <QVector>
#include <QWidget>
#include <QGroupBox>
#include <QTextEdit>
#include <QDebug>
 #include <QCheckBox>

#include "worker.h"
#include "parameter.h"

class TransWindow : public QWidget
{
    Q_OBJECT
public:
    explicit TransWindow(QWidget *parent = 0);
    void bindWorker(Worker*);
    Worker *worker;
    TransitionRates newRates;
    void setOldRates();
protected:
    void closeEvent(QCloseEvent *);
private:
    QCheckBox *inhomoPlusCheck;
    QLineEdit *ploadEdit;
    QLabel *infoLabel;
    QLineEdit *stepSizeEdit;
    QVector<QLineEdit*> bundleRegionEdit;
    QLineEdit *hyphalLengthEdit;
    QVector<QLineEdit*> velocityEdit;
    QVector<QLineEdit*> runUnbundledEdit;
    QVector<QLineEdit*> runBundledEdit;
    QVector<QLineEdit*> turnuniEdit;
    QVector<QLineEdit*> turnbiEdit;
    QVector<QLineEdit*> turninhomoEdit;
    QLineEdit* crossEdit;
    QVector<QLineEdit*> forwardEdit;
    QVector<QLineEdit*> laneChangeEdit;
    QVector<QLineEdit*> switchInnerEdit;
    QVector<QLineEdit*> switchEndEdit;
    QVector<QLineEdit*> forwardEndEdit;
    QVector<QLineEdit*> injectEdit;
    QVector<QLineEdit*> exitEdit;
    QLabel *setLabel;
    QLabel *turnInfLabel;
    QLabel *stepSizeLabel;
    QGroupBox *singleRateBox;
    QGroupBox *dataBox;
    QGroupBox *boundaryBox;
    QGroupBox *switchBox;
    QGroupBox *inhomoBox;

    void getNewRates();
    void bundleLengthMessage();

signals:
    
public slots:
    void checkRates();
    void getForward();
    void getTurning();
    void setDefault();// set default rates from the initializer
    void inhomoPlusChanged();
    void getRunLength();
    void getVelocity();
};

#endif // TRANSWINDOW_H
