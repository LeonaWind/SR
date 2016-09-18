#include"sr.h"

float final_x=0,final_y=0;

void registering(Mat image, Mat &temp, float x, float y, int mx, int my)
{
    int dx, dy;
    int color;

    uchar *image_ptr = image.data;
    uchar *temp_ptr = temp.data;
    int image_step = image.step;
    int temp_step = temp.step;

    //把其他图像映射到网格中
    for (int i = abs(x); i<image.rows - abs(x); i++)
    {
        for (int j = abs(y); j<image.cols - abs(y); j++)
        {
            for (int h = 0; h<3; h++)
            {
                dx = (int)(i + x)*my;
                dy = (int)(j + y)*mx;
                if (dx < 0)dx = 0;
                if (dy < 0)dy = 0;
                if (dx>(temp.rows))dx = temp.rows;
                if (dy>(temp.cols))dy = temp.cols;
                color = *(temp_ptr + dx*temp_step + dy * 3 + h);//原来的像素点的值
                if (color != 0)
                {
                    *(temp_ptr + dx*temp_step + dy * 3 + h) = (int)((color + *(image_ptr + i*image_step + j * 3 + h)) / 2 + 0.5);//均分
                }
                else
                {
                    *(temp_ptr + dx*temp_step + dy * 3 + h) = *(image_ptr + i*image_step + j * 3 + h);//填入
                }
            }
        }
    }

}

void spline(int x[], double y[], int n, double yp1,double  ypn, double y2[])
    //yi=f(xi),插值函数在第0个点和n-1个点的一阶导数值yp1和ypn。本程序返回y2[]其值为插值函数在xi点的二阶导数。
    //若yp1或ypn大于或等于10^30,则本程序将相应的边界条件设为自然样条，即边界的二阶导数为0,n为有效个数
{
    double u[win_size],aaa,sig,p,bbb,ccc,qn,un;
    int i,k;
    if( yp1 > 9.9e+29 )
    {
        y2[0] = 0;
        u[0] = 0;
    }
    else
    {
        y2[0] = -0.5;
        aaa = (y[1] - y[0]) / (x[1] - x[0]);
        u[0] = (3.0 / (x[1] - x[0])) * (aaa - yp1);
    }
    for (i = 1;i<n-1;i++)
    {
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        p = sig * y2[i - 1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        aaa = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
        bbb = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        ccc = x[i + 1] - x[i - 1];
        u[i] = (6.0 * (aaa - bbb) / ccc - sig * u[i - 1]) / p;
    }
    if (ypn > 9.9e+29 )
    {
        qn = 0.0;
        un = 0.0;
    }
    else
    {
        qn = 0.5;
        aaa = ypn - (y[n-1] - y[n - 2]) / (x[n-1] - x[n - 2]);
        un = (3.0/ (x[n-1] - x[n - 2])) * aaa;
    }
    y2[n-1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
    for (k = n - 2;k>=0;k--)
        y2[k] = y2[k] * y2[k + 1] + u[k];
}

void splint(int xa[win_size],double ya[win_size],double y2a[win_size],int x,double &y,int n)
    //给定数组xa和ya作为列表函数，y2a由spline求出，给定x的值，返回x的三次样条插值结果y,n为有效数值个数
{
    int k;
    double h,b,a;


    int klo=0;
    int khi=n-1;
    while(khi-klo>1)
    {
        k=(khi+klo)>>1;
        if(xa[k]>x) khi=k;
        else klo=k;
    }
    h=xa[khi]-xa[klo];
    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;
    y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi]*(h*h))/6.0;
}

void splie2(int x1a[], int x2a[], double ya[][win_size], int m, int n, double y2a[][win_size],int b)
    //给定一个m*n阶的列表函数ya[0…m-1][0…n-1]及列表自变量x2a，本程序构造ya的一维自然三次样条，返回的二阶导数存于数组y2a中
    //（x1a包含在变量列表中主要是为了保持一致）
{
    int x3a[win_size];//非零有效x坐标
    double ytmp[win_size], y2tmp[win_size];//ytmp[i]=f(x3a[i]),y2tmp[i]是x3a[i]的二阶导数
    int p=0;//用于有效值计数
    for (int j = 0; j<m; j++)
    {
        //初始化
        p=0;
        for(int h=0;h<n;h++)
        {
            ytmp[h] =0;
            x3a[h]=0;
            y2tmp[h]=0;
        }

        for (int k = 0; k<n; k++)//每一行选取有效的数值
        {
            if(ya[j][k]!=sigh)
            {
                ytmp[p] = ya[j][k];//要求解的有效数据
                x3a[p]=x2a[k];//相应y坐标
                p++;
            }
        }
        if(p>1)//至少1个数据，才能求二阶导数
            spline(x3a, ytmp, p, 1e+30, 1e+30, y2tmp);//返回第j行有效数据的二阶导数,放在 y2tmp中
        for (int k = 0; k<n; k++)
        {
            int u=(int)x3a[k];//获取有效的x坐标
            if(u!=0)
                y2a[j][u-b] = y2tmp[k];//二阶导数存入y2a,有效数据存在二阶导数，无效数据为0
            //u-b转换为数组中的有效坐标
        }
    }
}

void splie2win2(int x1a[], int x2a[], double ya[][win_size2], int m, int n, double y2a[][win_size2],int b)
    //给定一个m*n阶的列表函数ya[0…m-1][0…n-1]及列表自变量x2a，本程序构造ya的一维自然三次样条，返回的二阶导数存于数组y2a中
    //（x1a包含在变量列表中主要是为了保持一致）
{
    int x3a[win_size2];//非零有效x坐标
    double ytmp[win_size2], y2tmp[win_size2];//ytmp[i]=f(x3a[i]),y2tmp[i]是x3a[i]的二阶导数
    int p=0;//用于有效值计数
    for (int j = 0; j<m; j++)
    {
        //初始化
        p=0;
        for(int h=0;h<n;h++)
        {
            ytmp[h] =0;
            x3a[h]=0;
            y2tmp[h]=0;
        }

        for (int k = 0; k<n; k++)//每一行选取有效的数值
        {
            if(ya[j][k]!=sigh)
            {
                ytmp[p] = ya[j][k];//要求解的有效数据
                x3a[p]=x2a[k];//相应y坐标
                p++;
            }
        }
        if(p>1)//至少1个数据，才能求二阶导数
            spline(x3a, ytmp, p, 1e+30, 1e+30, y2tmp);//返回第j行有效数据的二阶导数,放在 y2tmp中
        for (int k = 0; k<n; k++)
        {
            int u=(int)x3a[k];//获取有效的x坐标
            if(u!=0)
                y2a[j][u-b] = y2tmp[k];//二阶导数存入y2a,有效数据存在二阶导数，无效数据为0
            //u-b转换为数组中的有效坐标
        }
    }
}


void fitting(Mat &temp)
{
    uchar *temp_ptr=temp.data;
    int temp_step=temp.step;

    //x1a存放x坐标,x2a存放y坐标,x3a存放筛选后符合要求的(x,j),(j=0,…,y-1)坐标
    int x1a[win_size], x2a[win_size],x3a[win_size];
    //ya存放(x,y)对应的z值，y2a存放(x,y)对应的二阶导数值，不存在为0

    int m,n,i,j;
    int x1,x2;//要求的中心点坐标
    double y;//最后求得的插值值
    int q=0;//非零有效数值计数
    int c1=0,c2=0;

    for(int a=0;a<temp.rows-win_size2;a++)
    {
        for(int b=0;b<temp.cols-win_size2;b++)
        {
            m = win_size;//窗口行高
            n = win_size;//窗口列宽
            x1=a+(win_size-1)/2;
            x2=b+(win_size-1)/2;//中心点坐标
            if(*(temp_ptr+ x1*temp_step +x2*3)==sigh)
            {
                double  ya[win_size][win_size], y2a[win_size][win_size]={0};
                //初始化
                for (i = 0; i<m; i++)
                    x1a[i]=a+i;
                for (i = 0; i<n; i++)
                    x2a[i] =b+i;
                for(int d=0;d<3;d++)//通道
                {
                    //输入插值数据
                    for(int i=0;i<m;i++)
                    {
                        for(int j=0;j<n;j++)
                        {
                            ya[i][j]=*(temp_ptr+ (a+i)*temp_step +(b+j)*3+d);
                        }
                    }

                    //根据值判断是否需要进行插值
                    float valueSum=0,valueNum=0,valueAve=0,diffValue=0;
                    for(int ii=0;ii<m;ii++)
                    {
                        for(int jj=0;jj<n;jj++)
                        {
                            if(ya[ii][jj]!=sigh)
                            {
                                valueSum+=ya[ii][jj];
                                ++valueNum;
                            }
                        }
                    }
                    if(valueNum!=0)
                    {
                        valueAve=valueSum/valueNum;//均值
                        valueSum=0;
                        for(int ii=0;ii<win_size;ii++)
                        {
                            for(int jj=0;jj<win_size;jj++)
                            {
                                if(ya[ii][jj]!=sigh)
                                {
                                    valueSum+=(ya[ii][jj]-valueAve);
                                }
                            }
                        }
                        diffValue=valueSum/valueNum;//均值差
                        if(diffValue>diffTh)//此区域相差较大，插值
                        {
                            splie2(x1a, x2a, ya, m, n, y2a,b);//返回ya中非零值的二阶导数


                            double ya_t[win_size]={0},y2a_t[win_size]={0},yytmp[win_size]={-1};
                            //ya_t存放有效ya,y2a_t存放有效二阶导数
                            //yytmp存放行插值结果
                            double yytmp1[win_size]={0},ytmp[win_size]={0};
                            //yytmp1存放有效的行插值结果,ytmp是其二阶导数
                            int x1a1[win_size]={0};//有效y坐标

                            for(int j=0;j<win_size;j++)//行
                            {
                                q=0;
                                for(int h=0;h<win_size;h++)
                                {
                                    ya_t[h]=0;
                                    y2a_t[h]=0;
                                    x3a[h]=0;
                                }

                                for(int k=0;k<win_size;k++)//列
                                {
                                    if(ya[j][k]!=sigh)//选取有效数字
                                    {
                                        ya_t[q]=ya[j][k];//ya_t是函数值
                                        y2a_t[q]=y2a[j][k];//y2a_t是二阶导数值
                                        x3a[q]=x2a[k];
                                        q++;
                                    }
                                }
                                if(q>1)//第j行至少有两个点，则进行插值
                                    splint(x3a,ya_t,y2a_t,x2,yytmp[j],q);//插值，并返回第j行x2处的样条插值的值y，存放在yytmp[j]中
                                else yytmp[j]=-1;
                            }
                            yytmp[(win_size-1)/2]=-1;//要插值的点=-1

                            q=0;
                            for(int i=0;i<win_size;i++)
                            {
                                if(yytmp[i]!=-1)
                                {
                                    yytmp1[q]=yytmp[i];
                                    x1a1[q]=x1a[i];
                                    q++;
                                }
                            }


                            spline(x1a1, yytmp1,q, 1e+30, 1e+30, ytmp);//求yytmp的二阶导数值

                            if(q>1)
                            {
                                splint(x1a1,yytmp1,ytmp,x1,y,q);
                                if(y>=0&&y<=255)//改变temp值
                                    *(temp_ptr+ x1*temp_step +x2*3+d)=(int)(y+0.5);
                                else if(y<0)
                                {
                                    *(temp_ptr+ x1*temp_step +x2*3+d)=0;
                                }
                                else if(y>255)
                                {
                                    *(temp_ptr+ x1*temp_step +x2*3+d)=255;
                                }
                                c1++;
                            }
                            else
                            {
                                *(temp_ptr+ x1*temp_step +x2*3+d)=(*(temp_ptr+ x1*temp_step +(x2-1)*3+d)+*(temp_ptr+ (x1-1)*temp_step +x2*3+d))/2;
                                c2++;
                            }
                        }
                        else//此区域相差不大，填充
                        {
                            *(temp_ptr+ x1*temp_step +x2*3+d)=(int)(valueAve+0.5);
                        }
                    }
                    else//这个区域内没有值
                    {
                        double  ya[win_size2][win_size2], y2a[win_size2][win_size2]={0};
                        //初始化
                        m = win_size2;//窗口行高
                        n = win_size2;//窗口列宽
                        x1=a+(win_size2-1)/2;
                        x2=b+(win_size2-1)/2;//中心点坐标
                        for (int i = 0; i<m; i++)
                            x1a[i]=a+i;
                        for (int i = 0; i<n; i++)
                            x2a[i] =b+i;
                        for(int d=0;d<3;d++)//通道
                        {
                            //输入插值数据
                            for(int i=0;i<m;i++)
                            {
                                for(int j=0;j<n;j++)
                                {
                                    ya[i][j]=*(temp_ptr+ (a+i)*temp_step +(b+j)*3+d);
                                }
                            }

                            splie2win2(x1a, x2a, ya, m, n, y2a,b);//返回ya中非零值的二阶导数


                            double ya_t[win_size2]={0},y2a_t[win_size2]={0},yytmp[win_size2]={-1};
                            //ya_t存放有效ya,y2a_t存放有效二阶导数
                            //yytmp存放行插值结果
                            double yytmp1[win_size2]={0},ytmp[win_size2]={0};
                            //yytmp1存放有效的行插值结果,ytmp是其二阶导数
                            int x1a1[win_size2]={0};//有效y坐标

                            for(int j=0;j<win_size2;j++)//行
                            {
                                q=0;
                                for(int h=0;h<win_size2;h++)
                                {
                                    ya_t[h]=0;
                                    y2a_t[h]=0;
                                    x3a[h]=0;
                                }

                                for(int k=0;k<win_size2;k++)//列
                                {
                                    if(ya[j][k]!=sigh)//选取有效数字
                                    {
                                        ya_t[q]=ya[j][k];//ya_t是函数值
                                        y2a_t[q]=y2a[j][k];//y2a_t是二阶导数值
                                        x3a[q]=x2a[k];
                                        q++;
                                    }
                                }
                                if(q>1)//第j行至少有两个点，则进行插值
                                    splint(x3a,ya_t,y2a_t,x2,yytmp[j],q);//插值，并返回第j行x2处的样条插值的值y，存放在yytmp[j]中
                                else yytmp[j]=-1;
                            }
                            yytmp[(win_size2-1)/2]=-1;//要插值的点=-1

                            q=0;
                            for(int i=0;i<win_size2;i++)
                            {
                                if(yytmp[i]!=-1)
                                {
                                    yytmp1[q]=yytmp[i];
                                    x1a1[q]=x1a[i];
                                    q++;
                                }
                            }


                            spline(x1a1, yytmp1,q, 1e+30, 1e+30, ytmp);//求yytmp的二阶导数值

                            if(q>1)
                            {
                                splint(x1a1,yytmp1,ytmp,x1,y,q);
                                if(y>=0&&y<=255)//改变temp值
                                    *(temp_ptr+ x1*temp_step +x2*3+d)=(int)(y+0.5);
                                else if(y<0)
                                {
                                    *(temp_ptr+ x1*temp_step +x2*3+d)=0;
                                }
                                else if(y>255)
                                {
                                    *(temp_ptr+ x1*temp_step +x2*3+d)=255;
                                }
                                c1++;
                            }
                            else
                            {
                                *(temp_ptr+ x1*temp_step +x2*3+d)=(*(temp_ptr+ x1*temp_step +(x2-1)*3+d)+*(temp_ptr+ (x1-1)*temp_step +x2*3+d))/2;
                                c2++;
                            }
                        }
                    }
                }
            }
        }//end 行
    }//end 列

    //裁剪掉未拟合的黑色区域
    temp= temp(Rect((win_size-1)/2+10,(win_size-1)/2+10,temp.cols-win_size-20,temp.rows-win_size-20));
}

float iMAD(int step,int curr_step,unchar* prev_block,unchar* curr_block,int m)//计算两块的平均绝对帧差MAD
{
    int temp=0;
    int value=0;
    int sum=0;
    float mad_value=0;
    for(int i=-2;i<3;i++)
    {
        for(int j=-2;j<3;j++)
        {
            temp=*(curr_block+i*curr_step+j)-*(prev_block+i*step+j);//？
            sum+=abs(temp);//绝对值
        }
    }
    mad_value=(float)sum/(D*D);
    return mad_value;
}

void motion(Mat prev, Mat curr,float &motion_x, float &motion_y)//运动估计
{
    Size size = prev.size();//读取prev大小，int
    int x,y;
    float p_x,p_y;
    Mat tempMax,prevMax;
    IplImage* gray;
    IplImage* gradient;
    double ux = 0;
    double uy = 0;
    int tempstep;

    int record[120][120] = { 0 };
    int ii, jj, max = 0;
    float rx, ry;

    //动态申请二维数组存放数据
    float *data1=new float[2*prev.cols];
    float (*data)[2]= (float(*)[2])data1;

    for(int i=0;i<prev.cols;i++)
        for(int j=0;j<2;j++)
            data[i][j]=0;

    //result=Mat(Size(size.width,size.height), prev.type(), Scalar::all(0));

    //unchar *ipresult=(unchar*)(result.data);

    unchar *ipprev=(unchar*)(prev.data );//开始的位置
    unchar *ipcurr=(unchar*)(curr.data );

    //prev一行字节数/unchar
    int step= prev.step/sizeof(uchar);

    //初始化
    /*for (int i=0; i<result.rows ; i++)
        for (int j=0;j<result.cols ; j++)
             *(ipresult+i*result.step+j)=225;*/

    unchar* prev_block;
    unchar* curr_block;
    unchar* temp_block;
    int It;



    for (int m = D; m < (size.height - D); m++)//不考虑边界，行
    {
        for (int n = D; n < (size.width - D); n++)//列
        {
            IplImage temp1 = (IplImage)prev;
            IplImage *gray=&temp1;

            IplImage temp2 = (IplImage)curr;
            IplImage *gradient=&temp2;

            MySobel(gray, gradient, m, n, ux, uy);

            It = ((ipcurr + m*step + n) - (ipprev + m*step + n) +
                (ipcurr + (m + 1)*step + n) - (ipprev + (m + 1)*step + n) +
                (ipcurr + m*step + n + 1) - (ipprev + m*step + n + 1) +
                (ipcurr + (m + 1)*step + n + 1) - (ipprev + (m + 1)*step + n + 1)) / 4;

            if ((ux + uy) > Threshold&&It > Threshold){
                final_x = 0; final_y = 0;
                //第一次搜索
                prev_block = ipprev + m*step + n;//取出一个像素放到prev_block中去，起始+第几行*行块+第几列
                curr_block = ipcurr + m*step + n;
                temp_block = ipcurr + m*step + n;
                tempstep = step;
                find(step, tempstep, prev_block, temp_block, 1);//第一次搜索
                p_x = final_x; p_y = final_y;//第一次搜索找到的偏移
                x = m + p_x; y = n + p_y;
                temp_block = ipcurr + y + x*step;//相对原图像，第一次搜索找到的位置

                //找到偏移,第二次搜索
                tempMax = interpolation(step, temp_block - mD*step - mD, D, 2 * D);//以第一次搜索找到的位置为中心插值
                prevMax = interpolation(step, prev_block - mD*step - mD, D, 2 * D);//原图像插值,(n,m)为中心点
                prev_block = prevMax.data + (prevMax.rows / 2 - 1)*prevMax.step + prevMax.cols / 2 - 1;//插值后图像中心
                curr_block = tempMax.data + (tempMax.rows / 2 - 1)*tempMax.step + tempMax.cols / 2 - 1;
                find(prevMax.step, tempMax.step, prev_block, curr_block, 2);//第二次搜索
                p_x += final_x / 2; p_y += final_y / 2;//第一次搜索偏移+第二次搜索偏移
                x = final_x + tempMax.rows / 2 - 1; y = final_y + tempMax.cols / 2 - 1;
                temp_block = curr_block + y + x*tempMax.step;//相对第一次插值后图片，搜索找到的位置

                //找到偏移,第三次搜索
                tempstep = tempMax.step;
                tempMax = interpolation(tempstep, temp_block - mD* tempstep - mD, D, 2 * D);//以第二次搜索找到的位置为中心插值
                prevMax = interpolation(prevMax.step, prev_block - mD*prevMax.step - mD, D, 2 * D);//原图像插值,4倍
                prev_block = prevMax.data + (prevMax.rows / 2 - 1)*prevMax.step + prevMax.cols / 2 - 1;//插值后图像中心
                curr_block = tempMax.data + (tempMax.rows / 2 - 1)*tempMax.step + tempMax.cols / 2 - 1;
                find(prevMax.step, tempMax.step, prev_block, curr_block, 4);//第三次搜索
                p_x += final_x / 4; p_y += final_y / 4;//第一次搜索偏移+第二次搜索偏移+第三次搜索偏移
                x = final_x + tempMax.rows / 2 - 1; y = final_y + tempMax.cols / 2 - 1;
                temp_block = curr_block + y + x*tempMax.step;//第三次搜索找到的位置
                /*
                if (p_x != 0 && p_y != 0 && m % 5 == 0 && n % 20 == 0)//屏幕打印移动轨迹
                {
                    line(result, cvPoint(n, m), cvPoint(n + p_x, m + p_y), CV_RGB(0, 0, 0), 1, 8, 0);

                }*/
                rx = p_x * 4;
                ry = p_y * 4;
                ii = (int)(rx + 60);
                jj = (int)(ry + 60);
                ++record[ii][jj];
            }//列
        }//行
    }

    for (int i = 0; i < 120; i++)
    {
        for (int j = 0; j < 120; j++)
        {
            if ((i != 60 || j != 60) && record[i][j] > max)
            {
                rx = (float)i - 60.0;
                ry = (float)j - 60.0;
                motion_x = rx / 4.0;
                motion_y = ry / 4.0;
                max = record[i][j];
            }
        }
    }
}

void MySobel(IplImage* gray, IplImage* gradient, int i, int j, double &ux, double &uy)
{
    unsigned char a00, a01, a02;
    unsigned char a10, a11, a12;
    unsigned char a20, a21, a22;

    CvScalar color;

    a00 = cvGet2D(gray, i - 1, j - 1).val[0];
    a01 = cvGet2D(gray, i - 1, j).val[0];
    a02 = cvGet2D(gray, i - 1, j + 1).val[0];
    a10 = cvGet2D(gray, i, j - 1).val[0];
    a11 = cvGet2D(gray, i, j).val[0];
    a12 = cvGet2D(gray, i, j + 1).val[0];
    a20 = cvGet2D(gray, i + 1, j - 1).val[0];
    a21 = cvGet2D(gray, i + 1, j).val[0];
    a22 = cvGet2D(gray, i + 1, j + 1).val[0];
    // x方向上的近似导数
    ux = a20 * (1) + a21 * (2) + a22 * (1)
        + (a00 * (-1) + a01 * (-2) + a02 * (-1));
    // y方向上的近似导数
    uy = a02 * (1) + a12 * (2) + a22 * (1)
        + a00 * (-1) + a10 * (-2) + a20 * (-1);

}

void find(int step,int search_step,unchar* prev_block,unchar* curr_block,int m)
{

    unchar* temp_block;
    int x=mD,y=mD;
    float temp,min;
    unchar* ip=curr_block-mD*search_step-mD;
    float t1=0,t2=0;

    temp_block=ip+x+y*search_step;
    min=iMAD(step,search_step,prev_block,temp_block,m);

      for(int i=0;i<D;++i)
        for(int j=0;j<D;++j)
        {
            temp_block=ip+i+j*search_step;
            temp=iMAD(step,search_step,prev_block,temp_block,m);
            if(temp<min)
            {
                x=j;y=i;
                min=temp;
            }
        }


        final_x=x-mD;final_y=y-mD;
}

Mat interpolation(int step,unchar* curr_block,int winSize,int expSize)
{

    Mat matDst1 = Mat(cv::Size(expSize,expSize), CV_8UC1, cv::Scalar::all(0));
    /*双线性插值*/
    unchar* dataDst=matDst1.data;
    int stepDst=matDst1. step/sizeof(unchar);
    unchar* dataSrc=curr_block;
    int stepSrc=step;
    int iWidthSrc=winSize;
    int iHeightSrc=winSize;
    int cols=matDst1.cols;
    int rows=matDst1.rows;
    double scale_x = (double) iWidthSrc / cols;
    double scale_y = (double) iHeightSrc / rows;

    for(int j=0;j<matDst1.rows;++j)//行
    {
        float fy=(float)((j+0.5)*scale_y - 0.5);

        int sy=cvFloor(fy);//sy取不超过fy的整数
        fy-=sy;//fy取小数部分
        sy=min(sy,iHeightSrc -2);
        sy=max(0,sy);

        short cbufy[2];
        cbufy[0]=saturate_cast<short>((1.f-fy)*2048);//1.f C++中的浮点1
        cbufy[1]=2048-cbufy[0];

        for(int i=0;i<matDst1.cols;++i)
        {
            float fx=(float)((i+0.5)*scale_x - 0.5);
            int sx=cvFloor(fx);
            fx-=sx;

            if(sx<0)
            {
                 fx = 0;
                 sx = 0;
            }
            if (sx >= iWidthSrc - 1)
            {
                 fx = 0;
                 sx = iWidthSrc - 2;
             }
             short cbufx[2];
             cbufx[0] = cv::saturate_cast<short>((1.f - fx) * 2048);
             cbufx[1] = 2048 - cbufx[0];

                *(dataDst+ j*stepDst + i ) =
                (*(dataSrc + sy*stepSrc + sx ) * cbufx[0] * cbufy[0] +
                *(dataSrc + (sy+1)*stepSrc + sx ) * cbufx[0] * cbufy[1] +
                *(dataSrc + sy*stepSrc + sx+1 ) * cbufx[1] * cbufy[0] +
                *(dataSrc + (sy+1)*stepSrc + sx+1 ) * cbufx[1] * cbufy[1]) >>22;

                double x1,y1;
                x1=i%(int)(1/scale_x);
                y1=j%(int)(1/scale_y);

                if(x1==0&&y1==0)
                {
                    *(dataDst+ j*stepDst + i ) =	*(dataSrc + (j/(int)(1/scale_y))*stepSrc +(i/(int)(1/scale_x)) ) ;
                }
        }
    }
    return matDst1 ;
}
