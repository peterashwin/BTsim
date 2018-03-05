#ifndef MAINWINDOW_H
#define MAINWINDOW_H


#include <QMainWindow>
#include <QAction>
#include <QLabel>
#include <QMenu>
#include <QGraphicsView>
#include <QApplication>
#include <QFile>
#include <QDir>
#include <QVector>
#include <QImage>


#include "bundle.h"
#include "worker.h"
#include "setting.h"
#include "transwindow.h"
#include "assistant.h"
#include "plotwindow.h"
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    Setting *setWindow;
    TransWindow *setTransWindow;
    PlotWindow *plotWindow;

    void showImage();
    void bindWorker(Worker *wk);
    void setFlags();
    void setPlusLabels();
    void setBundle3DBackground();

    ~MainWindow();

protected:
    Worker *worker;
protected:
    void closeEvent(QCloseEvent *event);
    void paintEvent(QPaintEvent *event);

private:
    Ui::MainWindow *ui;
    QImage previousKymographImage;
    Assistant *assistant;
    QImage bundleImage;
    QPainter *painter;
    QLineEdit *nStepEdit;
    QGraphicsView * densityView;
    QVector<QLabel*> plusLabel;
    QLabel *imageLabel;
    QLabel *shrinkLabel;
    QLabel *kymographLabel;
    QLabel *timeLabel;
    QLabel *timeLabel1;
    QLabel *kPercentageLabel;
    QMenu *fileMenu;
    QMenu *recordMenu;
    QMenu *dataMenu;
    QMenu *imageMenu;
    QMenu *importMenu;
    QMenu *currentMenu;
    QMenu *densityMenu;
    QMenu *runlengthMenu;
    QMenu *colocalizeMenu;
    QMenu *pauseMenu;
    QMenu *tipSizeMenu;
    QMenu *projectMenu;
    QMenu *setUpMenu;
    QMenu *calculateMenu;
    QMenu *initializerMenu;
    QMenu *blockageMenu;
    QMenu *endColorMenu;
    QMenu *paraMenu;
    QMenu *displayMenu;
    QMenu *helpMenu;

    QAction *exitAct;
    QAction *saveAct;
    QAction *stopAct;
    QAction *exportAct;
    QAction *importGeomAct;
    QAction *importCompleAct;
    QAction *recordAct;
    QAction *recordBundleAct;
    QAction *recordMTAct;
    QAction *recordKymographAct;
    QAction *stopRecordAct;
    QAction *openBoundaryAct;
    QAction *startDensityAct;
    QAction *showDensityAct;
    QAction *stopDensityAct;
    QAction *savePlotsAct;
    QAction *saveDensityAct;
    QAction *startCurrentAct;
    QAction *showCurrentAct;
    QAction *saveCurrentAct;
    QAction *runLengthFromEndAct;
    QAction *runLengthFromBundleAct;
    QAction *showDistributionAct;
    QAction *saveRunLengthAct;
    QAction *colocalizeTipAct;
    QAction *saveColocalizeAct;
    QAction *pauseAct;
    QAction *savePauseAct;
    QAction *checkAct;
    QAction *tipRegionAct;
    QAction *startTipSizeAct;
    QAction *showTipsizeAct;
    QAction *saveTipSizeAct;
    QAction *eeDistanceAct;
    QAction *eeNumAct;
    QAction *startEEDistriScoreAct;
    QAction *showEEDistriScoreAct;
    QAction *saveEEStateAct;
    QAction *particleNumAct;
    QAction *setUpBundleAct;
    QAction *addTrackAct;
    QAction *removeTrackAct;
    QAction *changeTrackAct;
    QAction *breakTrackAct;
    QAction *geometryShowAct;
    QAction *checkBundleAccessAct;
    QAction *setUnaccessibleAct;
    QAction *showBlockageAct;
    QAction *activationAct;
    QAction *activateTipAct;
    QAction *blockNumAct;
    QAction *blockOutBundleAct;
    QAction *addBlockAct;
    QAction *showPMendAct;
    QAction *threeDimensionEffctAct;
    QAction *showBundleAct;
    QAction *showMTAct;
    QAction *showKymographAct;
    //QAction *inhomoAct;
    QAction *transitionRateAct;

    QAction *defaultInitializerAct;
    QAction *oneendInitializerAct;
    QAction *tipendInitializerAct;
    QAction *randomInitializerAct;
    QAction *kymographSpeedAct;

    QAction *aboutAct;
    QAction *helpAct;
    QFile *f;
    QDir *newfolder;
    QString dirname;

    void createActions();
    void createMenus();
    void createStatusBar();
    void recordKymograph(QFile*);
    void showInitialDisplay();
    void initialKymographImage();
    void initialKymographImageTimeDown();
    void creatDir();
    void chooseInitialize();
    void test();
    bool setLabelRegion();
    void setActChecked();

private slots:
    void breakTrack();
    void startDistributionScore();
    void showDistributionScore();
    void saveEEState();
    void showEEdistance();
    void showEENum();//per interval length; accumulated
    void setOpenBoundary();
    void randomInitialize();
    void oneendInitialize();
    void tipendInitialize();
    void defaultInitialize();
    void showLatticeBackup();
    void showLattice();
    void showShrink();
    void showKymograph();
    void showKymographTimeDown();
    void lotsUpdate();
    void startRun();
    void stopRun();
    void exportBundle();
    void exportBundleGeometry();
    void importBundle();
    void importGeometry();
    void quitMw();
    void saveKymograph();//save data for kymograph
    void stopSave();
    void recordBundle();
    void recordMT();
    void recordKymograph();
    void startDensityCalculate();
    void showDensity();
    void savePlots();
    void startCurrentCalculate();
    void showCurrent();
    void colocalize();
    void particleNumChange();
    void setUpBundle();
    void changeTrack();
    void addTrack();
    void removeTrack();
    void showSetupWindow();
    void setTransitionRates();
    void blockChanged();// display the blockage in the bundle view
    void blockOutBundle();
    void setBlockNum();
    void addBlock();
    void showPMendChanged();
    //void inhomoChanged();
    void setSpeed();// set speed in showing kymograph
    void about();
    void showdocumentation();
    void show3DEffect();
    void showDistribution();
    void runLengthFromEnd();
    //void runLengthFromBundle();
    void saveRunLength();
    void check();
    void showPMEnds();
    void showBlockages();
    void activation();// for single site
    void activateTip();
    void pauseCalculate();
    void savePause();
    void saveDensity();
    void saveCurrent();
    void saveColocalize();
    void setTipRegion();
    void startTipSize();
    void showTipSize();
    void saveTipSize();
    void checkBundleAccess();

    void setUnaccessible();
};



#endif // MAINWINDOW_H
