#-------------------------------------------------
#
# Project created by QtCreator 2016-03-17T21:05:01
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = image
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    sr.cpp \
    deblurring.cpp

HEADERS  += mainwindow.h \
    sr.h \
    deblurring.h

FORMS    += mainwindow.ui

INCLUDEPATH+=G:\image_qt\opencv\include\opencv\
             G:\image_qt\opencv\include\opencv2\
             G:\image_qt\opencv\include
LIBS+=G:\image_qt\opencv\lib\libopencv_calib3d249.dll.a\
        G:\image_qt\opencv\lib\libopencv_contrib249.dll.a\
        G:\image_qt\opencv\lib\libopencv_core249.dll.a\
        G:\image_qt\opencv\lib\libopencv_features2d249.dll.a\
        G:\image_qt\opencv\lib\libopencv_flann249.dll.a\
        G:\image_qt\opencv\lib\libopencv_gpu249.dll.a\
        G:\image_qt\opencv\lib\libopencv_highgui249.dll.a\
        G:\image_qt\opencv\lib\libopencv_imgproc249.dll.a\
        G:\image_qt\opencv\lib\libopencv_legacy249.dll.a\
        G:\image_qt\opencv\lib\libopencv_ml249.dll.a\
        G:\image_qt\opencv\lib\libopencv_objdetect249.dll.a\
        G:\image_qt\opencv\lib\libopencv_video249.dll.a
