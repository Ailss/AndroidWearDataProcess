%校准陀螺仪的代码
close all;
clear;
clc;

%1.读入数据
whichAxes='Z';
circleCount=10;
path=['/Users/ly/Documents/MATLAB/AndroidWearDataProcess/CalibrationGyro/data/',whichAxes,int2str(circleCount),'/'];
file=dir(path);

%2.提取gyro数据
readGyroData;

%3.画出gyro数据
plotGyroData

%4.积分
%找出积分区间，积分区间应该均大于0.3
integralData

