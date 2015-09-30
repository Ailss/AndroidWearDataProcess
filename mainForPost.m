% 主功能：第1-5部分获得相对姿态
% 辅助功能：第6部分产生惯导所需的数据，供mainForIntertialNavigation使用
%         第7部分可以进行一些画图验证
%         第8部分可以为Project:inmed2014-master产生原始数据

clear;
clc;
close all;

%%
%1. 写路径
% path='/Users/ly/Documents/MATLAB/watchDataProcess/';
path='/Users/ly/Documents/MATLAB/watchDataProcess/data/';
% which=1;
% switch(which)
%     case 1%Z
%         data='14_07_25_';
%         time='18_22_42_';
%     case 2%X
%         data='14_07_25_';
%         time='18_33_10_';
%     case 3%X
%         data='14_07_25_';
%         time='18_33_42_';
%     case 4%Y
%         data='14_07_25_';
%         time='18_58_37_';
%     case 5%Y
%         data='14_07_25_';
%         time='18_48_21_';
% end
% suffix='sensor_GWatch.txt';
% filepath=[path,data,time,suffix];
file=dir(path);
which=1;
filepath=[path,file(which+2).name];
%%
% 2. 对源数据分析
sensorTextReadFromDisk;
%%
% 3. 得到初始post
getAngle;
%%
% 4. 对数据处理，校准
dataAlignment2;
%%
% 5. 提取正式数据（8s后的数据，前8s是准备阶段）
timeLine8s=timeLine2(timeLine2>8);
time8s=timeLine(timeLine2>8);
post8s=post(timeLine2>8,:);
post8s=post8s-kron([postStart(:,1),0,0],ones(length(post8s),1));
acc8s=acc(timeLine2>8,:);
gyro8s=gyro(timeLine2>8,:);
% figure
% grid on
% hold on
% plot(time8s,post8s);
% % plot(timeLine2,post);
% 
% % result=angle2quat(post8s(:,1),post8s(:,2),post8s(:,3));
% % plot(time8s,result);
% hold off
save data.mat
%%
% % 6. 产生惯导所需的数据，供mainForIntertialNavigation使用
% % 存储惯导的测试数据
% initialAttitude=angle2quat(postStart(1),postStart(2),postStart(3));
% save /Users/ly/Documents/MATLAB/armEstimation/INS/data/dataAccAndGyro.mat acc8s gyro8s time8s timeLine8s initialAttitude;

%%
% 7. 测试部分，画图
% figure
% executeWhichPara='post';
% hold on
% eval(['plot(timeLine2,',executeWhichPara,'(:,1),''r.'')']);
% eval(['plot(timeLine2,',executeWhichPara,'(:,2),''b.'')']);
% eval(['plot(timeLine2,',executeWhichPara,'(:,3),''k.'')']);
% grid
% axis([timeLine2(1),timeLine2(length(timeLine2)),-6,6]);
% plot(ones(1000,1).*3,-6:12/999:6);
% plot(ones(1000,1).*8,-6:12/999:6);
% hold off
% 
% figure
% executeWhichPara='gamepost';
% hold on
% eval(['plot(timeLine2,',executeWhichPara,'(:,1))']);
% eval(['plot(timeLine2,',executeWhichPara,'(:,2))']);
% eval(['plot(timeLine2,',executeWhichPara,'(:,3))']);
% grid
% hold off

%%
% 8. 为获得准确姿态做准备
% close all;
% Acc=acc;
% Gyro=gyro;
% Mag=mag;
% TimeLine=timeLine2;
% clear firstSensorEventTime postStart timeLine earth gamepost old ori post post8s rawData time8s
% clear acc gyro mag timeLine2
% 
% save /Users/ly/Documents/MATLAB/inmed2014-master/Scripts/data.mat