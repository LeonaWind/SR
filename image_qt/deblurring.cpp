#include"deblurring.h"
/*************************************************************************
 *
 * \函数名称：
 *   FFT_1D()
 *
 * \输入参数:
 *   complex<double> * pCTData	- 指向时域数据的指针，输入的需要变换的数据
 *   complex<double> * pCFData	- 指向频域数据的指针，输出的经过变换的数据
 *   nLevel						－傅立叶变换蝶形算法的级数，2的幂数，
 *
 * \返回值:
 *   无
 *
 * \说明:
 *   一维快速傅立叶变换。
 *
 *************************************************************************
 */

void FFT_1D(complex<double> * pCTData, complex<double> * pCFData, int nLevel)
{
    // 循环控制变量
    int		i;
    int     j;
    int     k;

    double PI = 3.1415926;

    // 傅立叶变换点数
    int	nCount =0 ;

    // 计算傅立叶变换点数
    nCount =(int)pow(2,nLevel) ;

    // 某一级的长度
    int		nBtFlyLen;
    nBtFlyLen = 0 ;

    // 变换系数的角度 ＝2 * PI * i / nCount
    double	dAngle;

    complex<double> *pCW ;

    // 分配内存，存储傅立叶变化需要的系数表
    pCW  = new complex<double>[nCount / 2];

    // 计算傅立叶变换的系数
    for(i = 0; i < nCount / 2; i++)
    {
        dAngle = -2 * PI * i / nCount;
        pCW[i] = complex<double> ( cos(dAngle), sin(dAngle) );
    }

    // 变换需要的工作空间
    complex<double> *pCWork1,*pCWork2;

    // 分配工作空间
    pCWork1 = new complex<double>[nCount];

    pCWork2 = new complex<double>[nCount];


    // 临时变量
    complex<double> *pCTmp;

    // 初始化，写入数据
    memcpy(pCWork1, pCTData, sizeof(complex<double>) * nCount);

    // 临时变量
    int nInter;
    nInter = 0;

    // 蝶形算法进行快速傅立叶变换
    for(k = 0; k < nLevel; k++)
    {
        for(j = 0; j < (int)pow(2,k); j++)
        {
            //计算长度
            nBtFlyLen = (int)pow( 2,(nLevel-k) );

            //倒序重排，加权计算
            for(i = 0; i < nBtFlyLen/2; i++)
            {
                nInter = j * nBtFlyLen;
                pCWork2[i + nInter] =
                    pCWork1[i + nInter] + pCWork1[i + nInter + nBtFlyLen / 2];
                pCWork2[i + nInter + nBtFlyLen / 2] =
                    (pCWork1[i + nInter] - pCWork1[i + nInter + nBtFlyLen / 2])
                    * pCW[(int)(i * pow(2,k))];
            }
        }

        // 交换 pCWork1和pCWork2的数据
        pCTmp   = pCWork1	;
        pCWork1 = pCWork2	;
        pCWork2 = pCTmp		;
    }

    // 重新排序
    for(j = 0; j < nCount; j++)
    {
        nInter = 0;
        for(i = 0; i < nLevel; i++)
        {
            if ( j&(1<<i) )
            {
                nInter += 1<<(nLevel-i-1);
            }
        }
        pCFData[j]=pCWork1[nInter];
    }

    // 释放内存空间
    delete pCW;
    delete pCWork1;
    delete pCWork2;
    pCW		=	NULL	;
    pCWork1 =	NULL	;
    pCWork2 =	NULL	;

}

/*************************************************************************
 *
 * \函数名称：
 *   FFT_2D()
 *
 * \输入参数:
 *   complex<double> * pCTData	- 图像数据
 *   int    nWidth				- 数据宽度
 *   int    nHeight				- 数据高度
 *   complex<double> * pCFData	- 傅立叶变换后的结果
 *
 * \返回值:
 *   无
 *
 * \说明:
 *   二维傅立叶变换。
 *
 ************************************************************************
 */
void DIBFFT_2D(complex<double> * pCTData, int nWidth, int nHeight, complex<double> * pCFData)
{

    // 循环控制变量
    int	x;
    int	y;

    // 临时变量
    double	dTmpOne;
    double  dTmpTwo;

    // 进行傅立叶变换的宽度和高度，（2的整数次幂）
    // 图像的宽度和高度不一定为2的整数次幂
    int		nTransWidth;
    int 	nTransHeight;

    // 计算进行傅立叶变换的宽度	（2的整数次幂）
    dTmpOne = log(nWidth)/log(2);
    dTmpTwo = ceil(dTmpOne)		   ;
    dTmpTwo = pow(2,dTmpTwo)	   ;
    nTransWidth = (int) dTmpTwo	   ;

    // 计算进行傅立叶变换的高度 （2的整数次幂）
    dTmpOne = log(nHeight)/log(2);
    dTmpTwo = ceil(dTmpOne)		   ;
    dTmpTwo = pow(2,dTmpTwo)	   ;
    nTransHeight = (int) dTmpTwo	   ;

    // x，y（行列）方向上的迭代次数
    int		nXLev;
    int		nYLev;

    // 计算x，y（行列）方向上的迭代次数
    nXLev = (int) ( log(nTransWidth)/log(2) +  0.5 );
    nYLev = (int) ( log(nTransHeight)/log(2) + 0.5 );

    for(y = 0; y < nTransHeight; y++)
    {
        // x方向进行快速傅立叶变换
        FFT_1D(&pCTData[nTransWidth * y], &pCFData[nTransWidth * y], nXLev);
    }

    // pCFData中目前存储了pCTData经过行变换的结果
    // 为了直接利用FFT_1D，需要把pCFData的二维数据转置，再一次利用FFT_1D进行
    // 傅立叶行变换（实际上相当于对列进行傅立叶变换）
    for(y = 0; y < nTransHeight; y++)
    {
        for(x = 0; x < nTransWidth; x++)
        {
            pCTData[nTransHeight * x + y] = pCFData[nTransWidth * y + x];
        }
    }

    for(x = 0; x < nTransWidth; x++)
    {
        // 对x方向进行快速傅立叶变换，实际上相当于对原来的图象数据进行列方向的
        // 傅立叶变换
        FFT_1D(&pCTData[x * nTransHeight], &pCFData[x * nTransHeight], nYLev);
    }

    // pCFData中目前存储了pCTData经过二维傅立叶变换的结果，但是为了方便列方向
    // 的傅立叶变换，对其进行了转置，现在把结果转置回来
    for(y = 0; y < nTransHeight; y++)
    {
        for(x = 0; x < nTransWidth; x++)
        {
            pCTData[nTransWidth * y + x] = pCFData[nTransHeight * x + y];
        }
    }

    memcpy(pCTData, pCFData, sizeof(complex<double>) * nTransHeight * nTransWidth );
}

/*************************************************************************
 *
 * \函数名称：
 *   IFFT_2D()
 *
 * \输入参数:
 *   complex<double> * pCFData	- 频域数据
 *   complex<double> * pCTData	- 时域数据
 *   int    nWidth				- 图像数据宽度
 *   int    nHeight				- 图像数据高度
 *
 * \返回值:
 *   无
 *
 * \说明:
 *   二维傅立叶反变换。
 *
 ************************************************************************
 */
void IFFT_2D(complex<double> * pCFData, complex<double> * pCTData, int nWidth, int nHeight)
{
    // 循环控制变量
    int	x;
    int	y;

    // 临时变量
    double	dTmpOne;
    double  dTmpTwo;

    // 进行傅立叶变换的宽度和高度，（2的整数次幂）
    // 图像的宽度和高度不一定为2的整数次幂
    int		nTransWidth;
    int 	nTransHeight;

    // 计算进行傅立叶变换的宽度	（2的整数次幂）
    dTmpOne = log(nWidth)/log(2);
    dTmpTwo = ceil(dTmpOne)		   ;
    dTmpTwo = pow(2,dTmpTwo)	   ;
    nTransWidth = (int) dTmpTwo	   ;

    // 计算进行傅立叶变换的高度 （2的整数次幂）
    dTmpOne = log(nHeight)/log(2);
    dTmpTwo = ceil(dTmpOne)		   ;
    dTmpTwo = pow(2,dTmpTwo)	   ;
    nTransHeight = (int) dTmpTwo	   ;

    // 分配工作需要的内存空间
    complex<double> *pCWork= new complex<double>[nTransWidth * nTransHeight];

    //临时变量
    complex<double> *pCTmp ;

    // 为了利用傅立叶正变换,可以把傅立叶频域的数据取共轭
    // 然后直接利用正变换，输出结果就是傅立叶反变换结果的共轭
    for(y = 0; y < nTransHeight; y++)
    {
        for(x = 0; x < nTransWidth; x++)
        {
            pCTmp = &pCFData[nTransWidth * y + x] ;
            pCWork[nTransWidth * y + x] = complex<double>( pCTmp->real() , -pCTmp->imag() );
        }
    }

    // 调用傅立叶正变换
    ::DIBFFT_2D(pCWork, nWidth, nHeight, pCTData) ;

    // 求时域点的共轭，求得最终结果
    // 根据傅立叶变换原理，利用这样的方法求得的结果和实际的时域数据
    // 相差一个系数
    for(y = 0; y < nTransHeight; y++)
    {
        for(x = 0; x < nTransWidth; x++)
        {
            pCTmp = &pCTData[nTransWidth * y + x] ;
            pCTData[nTransWidth * y + x] =
                complex<double>( pCTmp->real()/(nTransWidth*nTransHeight),
                                 -pCTmp->imag()/(nTransWidth*nTransHeight) );
        }
    }
    delete pCWork ;
    pCWork = NULL ;
}

Mat DIBBlindFilter(Mat pic)
{
    unsigned char *IpSrc;//指向源图像的指针
    long lWidth;//图像的宽度
    long lHeight;//图像的高度
    long lLineBytes;//图像每行字节数
    long i,j;
    double temp,tempre,tempim,a,b,c,d,e,f,norm1,norm2;//临时变量
    //实际进行傅里叶变换的宽度和高度
    long lw = 1;
    long lh = 1;
    int wp = 0;
    int hp = 0;

    //得到图像的宽度和高度

    cout<<"DIBBlindFilter"<<endl;
    //Mat Gray1 = imread("result\\fitting.bmp",0);
    //cvNamedWindow("src");
    //cvShowImage("src",Gray);
    int lHeight1 = pic.rows;
    int lWidth1 = pic.cols;
    int lLineBytes1 = pic.step;

    cout<<"lWidth1 "<<lWidth1<<"	lLineBytes1 "<<lLineBytes1<<endl;

     //保证离散傅里叶变换的宽度和高度为2的整数幂
    while(lw*2 <= lWidth1)
    {
        lw = lw*2;
        wp++;
    }
    while(lh*2 <= lHeight1)
    {
        lh = lh*2;
        hp++;
    }

    int max1;
    if(wp<=hp)max1=lw;
    else max1=lh;

    Size czSize;              //目标图像尺寸
    czSize.width =(int)max1;
    czSize.height =(int)max1;

    Mat Gray;//=Mat(czSize1, CV_8UC1, cv::Scalar::all(0));
    resize(pic, Gray, czSize,CV_INTER_AREA);
    lHeight = Gray.rows;
    lWidth= Gray.cols;
    lLineBytes = Gray.step;

    cout<<"lWidth "<<lWidth<<"	lLineBytes "<<lLineBytes<<endl;

    //用来存储源图像和变换核的时域数据
    complex<double> *pCTSrc,*pCTH;
    //用来存储源图像和变换核的频域数据，*pCFHnew新的滤波器
    complex<double> *pCFSrc,*pCFH,*pGF,*pCFHnew,*f1,*f2;

    //为时域和频域的数组分配空间
    pCTSrc=new complex<double> [lHeight*lWidth];
    pCTH=new complex<double> [lHeight*lWidth];
    pCFSrc=new complex<double> [lHeight*lWidth];
    pCFH=new complex<double> [lHeight*lWidth];
    pGF=new complex<double> [lHeight*lWidth];
    pCFHnew=new complex<double> [lHeight*lWidth];
    f1=new complex<double> [lHeight*lWidth];
    f2=new complex<double> [lHeight*lWidth];
    //滤波器加权系数
    //double *pCFFiletr=new double [lHeight*lWidth];
    complex<double> *pCTNoise=new complex<double> [lHeight*lWidth];
    complex<double> *pCFNoise=new complex<double> [lHeight*lWidth];

    for(int g=0;g<3;g++)
    {
     for (int j = 0;j < lHeight ;j++)
    {
        for(int i = 0;i < lWidth ;i++)
        {
            // 指向退化图像倒数第j行，第i个象素的指针
            IpSrc = (unsigned char *)Gray.data + lLineBytes * j + i*3+g;
            //将像素值存储到时域数组中
            pCTSrc [lWidth * j + i]=complex<double> ((double)*IpSrc,0);
            //频域赋零值
            pCFSrc [lWidth * j + i]=complex<double> (0.0,0.0);

            /*变换核*/
            //退化系统时域及维纳滤波加权系数赋值
            if(i<5&&j<5)
            {
                pCTH [lWidth * j + i]=complex<double> (0.04,0.0);
                //pCFFiletr[lWidth * j + i]=0.05;
            }
            else
            {
                pCTH [lWidth * j + i]=complex<double> (0.0,0.0);
                //pCFFiletr[lWidth * j + i]=0.025;
            }
            //频域赋零值
            pCFH [lWidth * j + i]=complex<double> (0.0,0.0);

            /*噪声*/
            if(i == j)
            {
                pCTNoise [lWidth * j + i]=complex<double> (0.1,0.0);
            }
            else
            {
                pCTNoise [lWidth * j + i]=complex<double> (0.0,0.0);
            }

        }
     }

     DIBFFT_2D (pCTSrc,lWidth,lHeight,pCFSrc);//对退化图像进行FFT
     pGF=pCFSrc ;
     f1=pCFSrc ;
     DIBFFT_2D (pCTH,lWidth,lHeight,pCFH);//对变换核图像进行FFT
     DIBFFT_2D (pCTNoise,lWidth,lHeight,pCFNoise);//对过滤器进行FFT

     double alf=0;
     //计算噪声功率谱
      for (int i = 0;i < lHeight*lWidth ;i++)
      {
          a=pCFNoise[i].real();
          b=pCFNoise[i].imag();
          alf=alf+a*a+b*b;
      }
      alf=alf/(lWidth/lHeight);//!!!
      //cout<<alf<<endl;
      //设置最大迭代次数k，并通过迭代更新G,H求出复原图像
      for(int k=1;k<10;k++)
      {
          cout<<"第"<<k<<"次迭代"<<endl;
         for (int i = 0;i < lHeight*lWidth ;i++)
         {
             //赋值
             a=f1[i].real();
             b=f1[i].imag();
             c=pCFH [i].real();
             d=pCFH [i].imag();
             e=pGF [i].real();
             f=pGF [i].imag();

             norm1=a*a+b*b;
             norm2=c*c+d*d;
             temp=norm1+(1/norm2);
             {
                 tempre=(a*e+b*f)/temp;
                 tempim=(a*f-b*e)/temp;
                 pCFHnew [i]=complex<double> (tempre,tempim);
            }

        }

          for (int i = 0;i < lHeight*lWidth ;i++)
         {
             //赋值
             a=f1[i].real();
             b=f1[i].imag();
             c=pCFH [i].real();
             d=pCFH [i].imag();
             e=pGF [i].real();
             f=pGF [i].imag();

             norm1=a*a+b*b;
             norm2=c*c+d*d;
             temp=(1/norm1)*alf+norm2;
             {
                 tempre=(c*e+d*f)/temp;
                 tempim=(c*f-d*e)/temp;
                 pCFSrc[i]=complex<double> (tempre,tempim);
            }

        }
          f1=pCFSrc ;
          pCFH =pCFHnew ;
      }
 IFFT_2D (pCFSrc,pCTSrc,lWidth,lHeight);

  for (int j = 0;j < lHeight ;j++)
    {
        for(int i = 0;i < lWidth ;i++)
        {
            // 指向退化图像倒数第j行，第i个象素的指针
            IpSrc = (unsigned char *)Gray.data + lLineBytes * j + i*3+g;
            a=pCTSrc[lWidth * j + i].real();
            b=pCTSrc[lWidth * j + i].imag();
            norm2=a*a+b*b;
            norm2=sqrt(norm2);//图像亮度

            if(norm2>255)
                norm2=255.0;
            if(norm2<0)
                norm2=0;

            *IpSrc =(unsigned char)(norm2);//改变的位置

        }
     }
  }

    Size czSize1;              //目标图像尺寸
    czSize1.width =pic.cols;
    czSize1.height =pic.rows;

    resize(Gray, Gray, czSize1,CV_INTER_CUBIC);

  //释放存储空间
  delete pCTSrc ;
  delete pCTH ;
  delete pCFSrc ;
  delete pCFH ;
  delete pCTNoise ;
  delete pCFNoise ;
  return Gray;





}
