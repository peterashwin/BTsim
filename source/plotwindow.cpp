/*
 *this plot method is based on customplot source file using QVector for the data sets
 */
#include "plotwindow.h"

PlotWindow::PlotWindow(QWidget *parent) :
    QWidget(parent)
{
    setWindowTitle(tr("plot"));
    frame = new QFrame();
    customPlot = new QCustomPlot(frame);
    customPlot->setObjectName(QStringLiteral("customPlot"));
    QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::MinimumExpanding);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(customPlot->sizePolicy().hasHeightForWidth());
    customPlot->setSizePolicy(sizePolicy);
    verticalLayout = new QVBoxLayout(frame);
    verticalLayout->addWidget(customPlot);

    customPlot->setRangeDrag(Qt::Horizontal|Qt::Vertical);
    customPlot->setRangeZoom(Qt::Horizontal|Qt::Vertical);
    customPlot->setupFullAxesBox();
    setLayout(verticalLayout);
}

void PlotWindow::closeEvent(QCloseEvent *event){
    customPlot->removeGraph(0);
}

void PlotWindow::plotDistribution(QVector<double> & bins, QVector<double> & values){

    customPlot->addGraph();
    customPlot->graph()->setData(bins, values);
    customPlot->graph()->setLineStyle((QCPGraph::LineStyle)(2));
    customPlot->graph()->setScatterSize(4);
    QPen graphPen;
    graphPen.setColor(Qt::black);
    graphPen.setWidthF(2);
    customPlot->graph()->setPen(graphPen);
}

void PlotWindow::plotState(QVector<double> & bins, QVector<double> & values, QPen & pen){
    customPlot->addGraph();
    customPlot->graph()->setData(bins, values);
    customPlot->graph()->setLineStyle(QCPGraph::lsNone);
    customPlot->graph()->setScatterStyle(QCP::ssCrossCircle);
    customPlot->graph()->setScatterSize(4);
    customPlot->graph()->setPen(pen);
}
