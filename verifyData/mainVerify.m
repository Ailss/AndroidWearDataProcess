close all;
clear
clc
whichaction='1';
path=['action',whichaction,'/data/'];
file=dir(path);
%whichOne='Acc';
whichOne='Gyro';

%将动作的ACC和GYRO读入
readAccGyro
clearvars path
load(['action',whichaction,'/action.mat'])

%画图
plotRes
