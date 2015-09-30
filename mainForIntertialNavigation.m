% 使用惯导的主程序，调用Project:armEstimation中INS文件下的Staged_SINS.m文件

clear
clc
close all
load /Users/ly/Documents/MATLAB/armEstimation/INS/data/dataAccAndGyro.mat

%%
% 整理数据，对接接口
initailPosition=[0 1 0];
initialVelocity=[0 0 0];
timeLine8s=timeLine8s-8;
timeLine8s(1)=0;
acc.time=timeLine8s;
acc.signals.values=acc8s(:,1:3);
gyro.time=timeLine8s;
gyro.signals.values=gyro8s;
clear acc8s gyro8s time8s timeLine8s

%%
% 调用惯导
sins=Staged_SINS(initialAttitude,initailPosition,initialVelocity,gyro,acc);