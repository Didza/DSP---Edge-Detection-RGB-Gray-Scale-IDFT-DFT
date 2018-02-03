    #include "mainwindow.h"
    #include "ui_mainwindow.h"
    #include "qmath.h"
    #include "math.h"
    #include<QVector>
    #include<iostream>
    #include <QFileDialog>
    #include <QColor>
    #include <QtMath>
    #include <QtMath>
    #include<complex.h>


    using namespace std;

    MainWindow::MainWindow(QWidget *parent) :
        QMainWindow(parent),
        ui(new Ui::MainWindow) {
        ui->setupUi(this);

        /*================ Create Connections of Signals & Slots ======================*/
        connect(ui->actionOpen_file, &QAction::triggered, this, &MainWindow::open_file);

    }

    MainWindow::~MainWindow() {
        delete ui;
    }

    void MainWindow::open_file() 
{
        /*================ Loads Image when open file ======================*/

        QString imageLoaded = QFileDialog::getOpenFileName();
        QImage image = QImage(imageLoaded);
        ui->Picture_Box->setPixmap(z.fromImage(image));

        w = image.width();
        h = image.height();
        red_Img      = QImage(w, h, QImage::Format_RGB32);
        green_Img    = QImage(w, h, QImage::Format_RGB32);
        blue_Img     = QImage(w, h, QImage::Format_RGB32);
        red_edges    = QImage(w, h, QImage::Format_RGB32);
        green_edges  = QImage(w, h, QImage::Format_RGB32);
        blue_edges   = QImage(w, h, QImage::Format_RGB32);
        all_edges    = QImage(w, h, QImage::Format_RGB32);
        /*=================== Sobel Algorithm ========================*/
        int Gx[][3] = {{-1,0,1},
                       {-5,0,5},
                       {-1,0,1}};

        int Gy[][3] = {{1,5,1},
                       {0,0,0},
                       {-1,-5,-1}};

        for (int i = 0; i < w; i++) {
            for(int j = 0; j < h; j++) {

         /*================== Get Component from image ===============*/

                green = (image.pixel(i,j)&0x00FF00) >>8;
                red   = (image.pixel(i,j)&0xFF0000) >>16;
                blue  = (image.pixel(i,j)&0x0000FF);

                green_Img.setPixel(i,j,green<<8);
                red_Img.setPixel(i,j,red<<16);
                blue_Img.setPixel(i,j,blue);

                int green_Xcomp     = 0;
                int green_Ycomp     = 0;
                int red_Xcomp       = 0;
                int red_Ycomp       = 0;
                int blue_Xcomp      = 0;
                int blue_Ycomp      = 0;
                int green_magnitude = 0;
                int red_magnitude   = 0;
                int blue_magnitude  = 0;

                if (i == 0 || i == w-1)
                {
                    green_magnitude = 0;
                    blue_magnitude  = 0;
                    red_magnitude   = 0;
                }

                else if (j == 0 || j == h-1)
                {
                    blue_magnitude = 0;
                    red_magnitude = 0;
                    green_magnitude = 0;
                }

                else {
                    for(int k = -1; k <= 1; k++) {
                        for(int l = -1; l <= 1; l++) {

                            green = (image.pixel(k + i, l + j) & 0x00FF00) >>8;
                            red   = (image.pixel(k + i, l + j) & 0xFF0000) >>16;
                            blue  = (image.pixel(k + i, l + j) & 0x0000FF);

                            blue_Xcomp += blue*Gx[k+1][l+1];
                            blue_Ycomp += blue*Gy[k+1][l+1];

                            red_Xcomp += red*Gx[k+1][l+1];
                            red_Ycomp += red*Gy[k+1][l+1];

                            green_Xcomp += green*Gx[k+1][l+1];
                            green_Ycomp += green*Gy[k+1][l+1];


                        }
                    }
                    red_magnitude   = sqrt(red_Xcomp*red_Xcomp + red_Ycomp*red_Ycomp);
                    blue_magnitude  = sqrt(blue_Xcomp*blue_Xcomp + blue_Ycomp*blue_Ycomp);
                    green_magnitude = sqrt(green_Xcomp*green_Xcomp + green_Ycomp*green_Ycomp);

                }

                if (blue_magnitude  > 255)
                {
                    blue_magnitude = 255;
                }

                if (red_magnitude   > 255)
                {
                    red_magnitude = 255;
                }

                if (green_magnitude > 255)
                {
                    green_magnitude = 255;
                }


                green_edges.setPixel(i, j, green_magnitude<<8);
                red_edges.setPixel(i, j, red_magnitude<<16);
                blue_edges.setPixel(i, j, blue_magnitude);
                all_edges.setPixel(i, j, blue_magnitude + (red_magnitude<<16) + (green_magnitude<<8) );
            }
        }


        ui->Blue_Edges->setPixmap(z.fromImage(blue_edges));
        ui->Red_Edges->setPixmap(z.fromImage(red_edges));
        ui->Green_Edges->setPixmap(z.fromImage(green_edges));
        ui->Edges->setPixmap(z.fromImage(all_edges));
        ui->Blue_Image->setPixmap(z.fromImage(blue_Img));
        ui->Red_Image->setPixmap(z.fromImage(red_Img));
        ui->Green_Image->setPixmap(z.fromImage(green_Img));
    }

    void MainWindow::DFT(QImage )
    {
        double pi_double = 2.0 * M_PI;
        double Width_New = w;
        double Height_New = h;

        DftNormalized_array.resize(Height_New);
        DftNormalized_array[0].resize(Width_New);

        ComplexDFT_array.resize(Height_New);
        ComplexDFT_array[0].resize(Width_New);

        for (int j = 0; j < Width_New-1; ++j)
        {

            for (int k = 0; k < Height_New-1; ++k)
            {
                ComplexDFT_array[j][k] = 0;
                ComplexDFT_array[k].resize(w);

                DftNormalized_array[j][k] = 0;
                DftNormalized_array[k].resize(w);

                for (int r = 0; r < Width_New-1; ++r)
                {
                    for (int s = 0; s < Height_New-1; ++s)
                    {
                        ComplexDFT_array[j][k] += inputSignal[r][s] * exp(complex<double> (0,-pi_double*((j*r)/Width_New + (k*s)/Height_New)));
                        DftNormalized_array[j][k]+= inputSignal[r][s] * qPow(-1.0,r+s ) * exp(complex<double> (0,-pi_double*((j*r)/Width_New + (k*s)/Height_New)));
                    }
                }
            }
           ui->progressBar->setValue(100*j/Width_New+2);
        }
        ui->progressBar->setValue(0);
    }

    void MainWindow::CalculatePhaseNAmp()
    {
        Phase_array.resize(w);
        Amplitude_array.resize(w);
        AmplitudeNormilizedArray.resize(w);
        for( int r = 0 ; r < w-1 ; r++)
        {
            Phase_array[r].resize(h);
            Amplitude_array[r].resize(h);
            for( int s = 0 ; s < h-1 ; s++)
            {
                Phase_array[r][s] =arg(DftNormalized_array[r][s])*255/(2*3.14) ;
                Amplitude_array[r][s]= abs(DftNormalized_array[r][s]);
            }
        }
        int maximumValue = 0 ;
        for(int r = 0 ; r < w-1 ; r++)
        {
            for(int s = 0 ; s < h-1 ; s++)
            {
                if( Amplitude_array[r][s] > maximumValue)
                {
                    maximumValue = Amplitude_array[r][s] ;
                }

            }
        }
        double NormalizedValue = 255/(log10(maximumValue+1));
        for(int r = 0 ; r < w-1 ; r++)
        {
            AmplitudeNormilizedArray[r].resize(h);

            for(int s = 0 ; s < h-1 ; s++)
            {
                AmplitudeNormilizedArray[r][s] = NormalizedValue*log10(1+Amplitude_array[r][s]);
            }
        }
    }

    QImage MainWindow::GrayScale(QImage LoadedImage)
    {
        greyImage = QImage(w,h,QImage::Format_RGB32);

        inputSignal.resize(w);
        for ( int r = 0; r < w ; r++ )
        {
            inputSignal[r].resize(h);
            for ( int s = 0; s < h ; s++ )
            {
                QColor clearImage(LoadedImage.pixel(r,s));

                int greyImageConvert = (clearImage.red() + clearImage.green() + clearImage.blue() )/3 ;
                QColor changedImage(greyImageConvert,greyImageConvert,greyImageConvert,255);

                greyImage.setPixel(r,s,changedImage.rgba());
                inputSignal[r][s] = greyImageConvert ;
            }
        }

        return greyImage;
    }

    void MainWindow::on_Grey_ScaleB_clicked()
    {


        filename = QFileDialog::getOpenFileName();
        QImage LoadGScaleImg(filename);
        w = LoadGScaleImg.width();
        h = LoadGScaleImg.height();
        ui->label_10->setPixmap(QPixmap::fromImage(LoadGScaleImg));
        ImageLoaded = LoadGScaleImg ;
        greyImage = GrayScale(LoadGScaleImg);

        ui->label_9->setPixmap(QPixmap::fromImage(greyImage));
    }
    void MainWindow::AmpPhasePlot()
    {
        Phase_image = QImage(w,h,QImage::Format_RGB32);
        Amplitude_image = QImage(w,h,QImage::Format_RGB32);

            for ( int r = 1; r < w-1 ; r++ )
            {
                for ( int s = 1; s < h-1 ; s++ )
                {
                    int Phase_integer = Phase_array[r][s];
                    int Amplitude_integer = AmplitudeNormilizedArray[r][s];
                    Phase_image.setPixel(r,s,qRgb(Phase_integer,Phase_integer,Phase_integer));
                    Amplitude_image.setPixel(r,s,qRgb(Amplitude_integer,Amplitude_integer,Amplitude_integer));
                }
            }
                ui->label_12->setPixmap(QPixmap::fromImage(Phase_image));
                ui->label_11->setPixmap(QPixmap::fromImage(Amplitude_image));
    }

    void MainWindow::IDFT()
    {
        double M_width = w;
        double New_height = h;
        int pixel_Value;
        complex<double> Point_D(0,0);
        double Pi_Double = 2.0 * M_PI;

        IDFTComplex_array.resize(New_height);
        IDFTComplex_array[0].resize(w);
        IDFT_image = greyImage ;

        for (int x = 0; x < M_width-1; ++x)
        {
            for (int y = 0; y < New_height-1; ++y)
            {
                IDFTComplex_array[x][y] = 0;
                IDFTComplex_array[y].resize(w);

                for (int r = 0; r < M_width-1; ++r)
                {
                    for (int s = 0; s < New_height-1; ++s)
                    {
                        Point_D = ComplexDFT_array[r][s];
                        IDFTComplex_array[x][y] += Point_D * exp(complex<double> (0,Pi_Double*((r*x)/M_width + (s*y)/New_height)));

                    }
                }
                IDFTComplex_array[x][y] *= 1/(M_width*New_height);
                pixel_Value = (int)(std::abs(IDFTComplex_array[x][y]));
                IDFT_image.setPixel(x,y,qRgb(pixel_Value,pixel_Value,pixel_Value));
            }
            ui->progressBar->setValue(100*x/M_width+2);
        }
        ui->label_13->setPixmap(QPixmap::fromImage(IDFT_image));
    }



    void MainWindow::on_DFT_clicked()
    {
        DFT(greyImage);
        CalculatePhaseNAmp();
        AmpPhasePlot();
    }

    void MainWindow::on_IDFT_clicked()
    {
        IDFT();
    }






