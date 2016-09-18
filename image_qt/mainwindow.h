#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QLabel>
#include<QMessageBox>
#include<QPushButton>
#include<QGridLayout>
#include<QFileDialog>
#include<QObject>
#include<iostream>
#include<QLabel>
#include<QImage>
#include<QPainter>
#include<QVector>
#include"sr.h"
#include"deblurring.h"


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
    QWidget *select_picture=0;
    QPushButton *select_btn;
    QPushButton *begin_btn;
    QLabel *select_label;
    QVector<QString> pathVector;
    Mat temp;
public slots:
    void open();
    void sr_start();
    void denoise();
    void deblurring();
    void save();
};

#endif // MAINWINDOW_H
