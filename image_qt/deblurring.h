#ifndef DEBLURRING
#define DEBLURRING

#include <opencv.hpp>
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include<math.h>
#include<time.h>


using namespace cv;
using namespace std;

void FFT_1D(complex<double> * pCTData, complex<double> * pCFData, int nLevel);
void DIBFFT_2D(complex<double> * pCTData, int nWidth, int nHeight, complex<double> * pCFData);
void IFFT_2D(complex<double> * pCFData, complex<double> * pCTData, int nWidth, int nHeight);
Mat DIBBlindFilter(Mat pic);

#endif // DEBLURRING
