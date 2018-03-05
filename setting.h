#ifndef SETTING_H
#define SETTING_H

#include <QWidget>
#include <QPushButton>
#include <QLineEdit>
#include <QSpinBox>
#include <QVBoxLayout>
#include <QLabel>
#include <QVector>
#include "worker.h"
class Setting : public QWidget
{
    Q_OBJECT
public:
    explicit Setting(QWidget *parent = 0);
    void bindWorker(Worker*);
    Worker *worker;
    void setUp();
private:
    QGridLayout *mainLayout;
    QHBoxLayout *topLeftLayout;
    QVBoxLayout *editLayout;
    QVector<QHBoxLayout*> vhboxLayout;
    QVector<QLineEdit*> vlengthEdit;
    QVector<QLineEdit*> vwidthEdit;
    QVector<QLineEdit*> vminusEdit;
    QVector<QLineEdit*> vplusEdit;
    QVector<QLineEdit*> vdisEdit;
    QVector<QLineEdit*> vdyneinEdit;
    QVector<QLabel*> vlabelEdit;
    QVector<QLineEdit*> vxcoorEdit;
    QVector<QLineEdit*> vorientationEdit;
    QSpinBox *trackBox;

};

#endif // SETTING_H
