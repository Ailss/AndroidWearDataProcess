% 校准加速度的主程序

close all;
clear;
clc;
%%
% 将data文件夹下的所有加速度数据合并，并进行采样
path='/Users/ly/Documents/MATLAB/watchDataProcess/data/calibrationData';
file=dir(path);
samplingRate=10;
jointAccData;

%%
AutoAccCalibration;