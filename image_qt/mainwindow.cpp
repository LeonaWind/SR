#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "sr.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    select_picture=new QWidget(this);
    setCentralWidget(select_picture);
    select_label=new QLabel();
    select_btn=new QPushButton("选择图片");
    begin_btn = new QPushButton("开始");

    QGridLayout *v = new QGridLayout();

    v->addWidget(select_label,0,0,1,2);
    v->addWidget(select_btn,1,0);
    v->addWidget(begin_btn,1,1);

    select_picture->setLayout(v);
    QObject::connect(select_btn, SIGNAL(clicked()), SLOT(open()));
    QObject::connect(begin_btn, SIGNAL(clicked()), SLOT(sr_start()));

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::open()
{
    QString path = QFileDialog::getOpenFileName(this, tr("Open Image"), ".", tr("Image Files(*.jpg *.png)"));
    QString show="选择的图片:\n";
    if(path.length() == 0)
    {
        select_label->setText("You didn't select any files.");
    }
    else
    {
        pathVector.append(path);

        QVector<QString>::iterator i;
        for (i = pathVector.begin(); i != pathVector.end(); ++i) {
          show+= (*i)+"\n"; // 使用 * 运算符获取遍历器所指的元素

        }
        select_label->setText(show);
    }
}

void MainWindow::sr_start()
{
    select_label->setText("请等待");
    //执行程序
    Mat prev, curr, prevc, currc;//curr指向对比图片，prev指向模板图片
    QVector<QString>::iterator iter=pathVector.begin();
    int i=0,j=0,k=0;
    string filename;//图片名
    int mx = 2;
    int my = 2;

    prev = imread((*iter).toStdString(), 0);
    prevc = imread((*iter).toStdString());

    temp = Mat(cv::Size((prev.cols * my), (prev.rows * mx)), CV_8UC3, cv::Scalar(2, 2, 2));

    uchar *image_ptr = prev.data;//指向第一个位置的指针
    uchar *temp_ptr = temp.data;
    int image_step = prev.step;//一行的像素数
    int temp_step = temp.step;
    uchar *imagec_ptr = prevc.data;//指向第一个位置的指针
    int imagec_step = prevc.step;//一行的像素数
    float x, y;//位移信息



    for (i = 0; i<prev.rows; i++)
        {
            for (j = 0; j<prev.cols; j++)
            {
                for (k = 0; k<3; k++)
                    *(temp_ptr + i * mx * temp_step + j * my * 3 + k) = *(imagec_ptr + i*imagec_step + j * 3 + k);
            }
        }

      for (iter=pathVector.begin();iter!= pathVector.end(); ++iter)
        {
           filename =(*iter).toStdString();
          //filename ="G:/image_qt/build-image-Desktop_Qt_5_4_0_MinGW_32bit-Debug/1.jpg";
            curr = imread(filename.c_str(), 0);//读取当前处理图片
            currc = imread(filename.c_str());
            /*运动估计*/
            motion(prev, curr, x, y);
            if(x<5&&x>-5&&y<5&&y>-5)
            {
                /*图像配准*/
                registering(currc, temp, x/2 , y/2 , mx, my);
            }
        }

        /*插值*/
        fitting(temp);
        imwrite("result.jpg", temp);//temp中存放配准后的矩阵
    select_label->setText("完成!");

    //显示结果
    QWidget *results_picture=new QWidget();
    QLabel *result_label=new QLabel();

    //改变图像大小
    QPixmap result_image("./result.jpg");
    QSize pixmapSize(500,500);
    result_image = result_image.scaled(pixmapSize,Qt::KeepAspectRatio);

    result_label->setPixmap(result_image);

    QPushButton *denoise_btn=new QPushButton("去噪");
    QPushButton *deblurring_btn = new QPushButton("去模糊");
    QPushButton *save_btn = new QPushButton("保存");

    QGridLayout *u = new QGridLayout();

    u->addWidget(result_label,0,0,1,3);
    u->addWidget(denoise_btn,1,0);
    u->addWidget(deblurring_btn,1,1);
    u->addWidget(save_btn,1,2);

    results_picture->setLayout(u);
    results_picture->show();
    //保存或其他处理
    QObject::connect(denoise_btn, SIGNAL(clicked()), SLOT(denoise()));
    QObject::connect(deblurring_btn, SIGNAL(clicked()), SLOT(deblurring()));
    QObject::connect(save_btn, SIGNAL(clicked()), SLOT(save()));
}

void MainWindow:: denoise()
{
        Mat noise = Mat(cv::Size(temp.cols , temp.rows ), CV_8UC1, cv::Scalar::all(1));//噪声检测矩阵
        double ave[4],t[4];//灰度平均值，判定阈值
        int tb=200;//阈值基值
        int noise_win=10;
        int half_win=(int)noise_win/2;
        uchar *temp_ptr=temp.data;
        int temp_step=temp.step;
        uchar *noise_ptr=noise.data;
        int noise_step=noise.step;
        int N;//非噪声点计数


        //计算噪声矩阵
        for(int a=noise_win+1;a<temp.rows-noise_win-1;a++)//行
        {
            for(int b=noise_win+1;b<temp.cols-noise_win-1;b++)//列
            {
                int x=*(temp_ptr+a*temp_step+b);

                N=0;
                for(int k=-1-noise_win;k<=-1;k++)//上邻域
                {
                    for(int l=-half_win;l<=half_win;l++)
                    {
                        ave[0]+=*(temp_ptr+(a+k)*temp_step+(b+l))**(noise_ptr+(a+k)*noise_step+(b+l));
                        N+=*(noise_ptr+(a+k)*noise_step+(b+l));
                    }
                }
                ave[0]=ave[0]/N;
                t[0]=tb*(1-abs(ave[0]-127.5)/127.5);



                if(x>ave[0]+t[0]||x<ave[0]-t[0])
                {
                    N=0;
                    for(int k=1;k<=1+noise_win;k++)//下邻域
                    {
                        for(int l=-half_win;l<=half_win;l++)
                        {
                            ave[1]+=*(temp_ptr+(a+k)*temp_step+(b+l))**(noise_ptr+(a+k)*noise_step+(b+l));
                            N+=*(noise_ptr+(a+k)*noise_step+(b+l));
                        }
                    }
                    ave[1]=ave[1]/N;
                    t[1]=tb*(1-abs(ave[1]-127.5)/127.5);


                    if(x>ave[1]+t[1]||x<ave[1]-t[1])
                    {
                        N=0;
                        for(int k=-half_win;k<=half_win;k++)//左邻域
                        {
                            for(int l=-1-noise_win;l<=-1;l++)
                            {
                                ave[2]+=*(temp_ptr+(a+k)*temp_step+(b+l))**(noise_ptr+(a+k)*noise_step+(b+l));
                                N+=*(noise_ptr+(a+k)*noise_step+(b+l));
                            }
                        }
                        ave[2]=ave[2]/N;
                        t[2]=tb*(1-abs(ave[2]-127.5)/127.5);


                        if(x>ave[2]+t[2]||x<ave[2]-t[2])
                        {
                            N=0;
                            for(int k=-half_win;k<=half_win;k++)//右邻域
                            {
                                for(int l=1;l<=1+noise_win;l++)
                                {
                                    ave[3]+=*(temp_ptr+(a+k)*temp_step+(b+l))**(noise_ptr+(a+k)*noise_step+(b+l));
                                    N+=*(noise_ptr+(a+k)*noise_step+(b+l));
                                }
                            }
                            ave[3]=ave[3]/N;
                            t[3]=tb*(1-abs(ave[3]-127.5)/127.5);


                            if(x>ave[3]+t[3]||x<ave[3]-t[3])
                            {
                                *(noise_ptr+a*noise_step+b)=0;//同时满足上述条件，为脉冲噪声
                            }//右邻域
                        }//左邻域
                    }//下邻域
                }//上邻域


            }
        }

        int c1=0,c2=0,c3=0,c4=0;
        //脉冲噪声滤除
        for(int a=noise_win+1;a<temp.rows-noise_win-1;a++)//行
        {
            for(int b=noise_win+1;b<temp.cols-noise_win-1;b++)//列
            {
                //cout<<(int)*(noise_ptr+a*noise_step+b)<<" ";
                vector<int> y;//非噪声集合
                if(*(noise_ptr+a*noise_step+b)==0)
                {
                    //system("pause");
                    for(int k=-half_win;k<=half_win;k++)//采集非噪声点
                    {
                        for(int l=-half_win;l<=half_win;l++)
                        {
                            if(*(noise_ptr+(a+k)*noise_step+(b+l))!=0)
                            {
                                y.push_back(*(temp_ptr+(a+k)*temp_step+(b+l)));
                            }
                        }
                    }

                    //排序
                    sort(y.begin(),y.end());


                    int m=y.size();
                    if(m==0)
                    {
                        ++c1;
                    }
                    else if(m%2==1)//奇数
                    {
                        ++c2;
                        *(temp_ptr+a*temp_step+b)=y[(m-1)/2];
                    }
                    else if(m%2==0)//偶数
                    {
                        ++c3;
                        *(temp_ptr+a*temp_step+b)=(y[m/2-1]+y[m/2])/2;
                    }
                }
                else
                {
                    c4++;
                }
            }
        }

        imwrite("./denoise.jpg", temp);//temp中存放配准后的矩阵
        imshow("./denoise", temp);
}

void MainWindow:: deblurring()
{
    temp = DIBBlindFilter(temp);
    imwrite("result\\DIBBlindFilter.bmp", temp);//temp中存放配准后的矩阵
    imshow("result\\DIBBlindFilter", temp);
}

void MainWindow:: save()
{

}
