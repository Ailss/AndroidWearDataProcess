close all;
clear
clc

%% 1.define basic varibales
path='txtdata/';
file=dir(path);
%whichOne='Acc';
whichOne='Gyro';
load('action.mat');

%% 2.read the txt data
readAccGyro
clearvars path

%% 3.format the data
addpath('../INS/')
for i=1:1:length(file)-2
    ii=num2str(i);
    eval(['Acc',ii,'=SimulinkTimeSeriesWrapper(Acc',ii,'(:,1),Acc',ii,'(:,2:1:4));']);
    eval(['Gyro',ii,'=SimulinkTimeSeriesWrapper(Gyro',ii,'(:,1),Gyro',ii,'(:,2:1:4));']);
end
clearvars i ii

%% 4.plot
plotRes
