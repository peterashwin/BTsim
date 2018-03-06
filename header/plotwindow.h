#ifndef PLOTWINDOW_H
#define PLOTWINDOW_H

#include "qcustomplot.h"
#include <QWidget>
#include <QFrame>
#include <QVBoxLayout>
#include <QVector>
#include <QPen>
class PlotWindow : public QWidget
{
    Q_OBJECT
public:
    explicit PlotWindow(QWidget *parent = 0);
    QFrame* frame;
    QVBoxLayout *verticalLayout;
    QCustomPlot *customPlot;

    void plotDistribution(QVector<double>& bins, QVector<double> & values);
    void plotState(QVector<double> & bins, QVector<double> & values, QPen & pen);
    void closeEvent(QCloseEvent *);
signals:
    
public slots:
    
};

#endif // PLOTWINDOW_H
