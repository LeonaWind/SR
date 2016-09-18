#ifndef SR
#define SR

#include <opencv.hpp>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include<vector>
#include <iomanip>

using namespace cv;
using namespace std;

#define  unchar unsigned char
#define D 9
#define mD 4
#define Threshold 110//运动估计筛选参数
#define unchar unsigned char
#define win_size 11
#define win_size2 21
#define sigh 2
#define diffTh 10

extern float final_x,final_y;
extern Mat result;

void registering(Mat image, Mat &temp, float x, float y, int mx, int my); //运动配准
void fitting(Mat &temp);//双三次样条插值
void splie2(int x1a[], int x2a[], double ya[][win_size], int m, int n, double y2a[][win_size],int b);
void splie2win2(int x1a[], int x2a[], double ya[][win_size2], int m, int n, double y2a[][win_size2],int b);
void splint(int xa[win_size],double ya[win_size],double y2a[win_size],int x,double &y,int n);
void spline(int x[], double y[], int n, double yp1,double  ypn, double y2[]);

void motion(Mat prev, Mat curr,float &x, float &y);//运动估计
Mat interpolation(int step,unchar* curr_block,int winSize,int expSize);//双线性插值
float iMAD(int curr_step,int step,unchar* prev_block,unchar* curr_block,int m);//计算两块的平均绝对帧差MAD
void find(int step,int search_step,unchar* prev_block,unchar* curr_block,int m);
void MySobel(IplImage* gray, IplImage* gradient, int i, int j, double &ux, double &uy);//sobel算子
#endif // SR

