#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QVector>
#include <complex>
#include <QLabel>

using namespace std ;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;

 QString filename;

 QImage ImageLoaded,red_Img,green_Img,blue_Img,red_edges,green_edges,
         blue_edges,all_edges;

 QImage image_Grey;

 QImage GrayScale( QImage LoadedImage );

 QImage greyImage , imgAmp , IDFT_image;

 QVector< QVector <double> > AmpVector, AmpVectorNormalized, PhaseVector, RealVector ,ImagineryVector ;
 void DFT( QImage );
 void CalculatePhaseNAmp() ;
 void AmpPhasePlot();
 QImage Phase_image , Amplitude_image;
 QImage PhaseImage , AmplitudeImage;

 QPixmap z;                                     /*Loaded image pixmap*/
 int w , h ;                                    /*image height and width*/
 int red,green,blue;                            /*color components integers*/

 QVector< QVector < complex <double> > > ComplexDFT_array , IDFTComplex_array , DftNormalized_array;
 QVector< QVector<double> > inputSignal;
 QVector< QVector <double> > Amplitude_array, AmplitudeNormilizedArray, Phase_array ;
 void PlotIDFT();
 void PopulateComplexMatrix();
 void IDFT( );


private slots:
    void open_file();
    void on_DFT_clicked();
    void on_Grey_ScaleB_clicked();
    void on_IDFT_clicked();
};

#endif // MAINWINDOW_H
