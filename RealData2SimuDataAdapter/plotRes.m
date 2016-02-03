close all
l=length(data.timeAcc);
if strcmp(whichOne,'Acc')
    eval(['plot(data.time',whichOne,'(30:1:l),data.data',whichOne,'(30:1:l,1),''r'',data.time',whichOne,'(30:1:l),data.data',whichOne,'(30:1:l,2),''g'',data.time',whichOne,'(30:1:l),data.data',whichOne,'(30:1:l,3),''k'');']);
else
    eval(['plot(data.time',whichOne,',data.data',whichOne,'(:,1),''r'',data.time',whichOne,',data.data',whichOne,'(:,2),''g'',data.time',whichOne,',data.data',whichOne,'(:,3),''k'');']);
    figure
    data.dataGyro=data.dataGyro/pi*180;
    integralRes=cumtrapz(data.timeGyro,data.dataGyro);
    plot(data.timeGyro,integralRes(:,1),'r',data.timeGyro,integralRes(:,2),'g',data.timeGyro,integralRes(:,3),'b');
end
legend('x','y','z');
title('simulation');

load('offsetGyro.mat')
for i=1:1:length(file)-2
    figure
    ii=num2str(i);
    if strcmp(whichOne,'Acc')
        eval(['plot(',whichOne,ii,'.time,',whichOne,ii,'.signals.values(:,1),''r'',',whichOne,ii,'.time,',whichOne,ii,'.signals.values(:,2),''g'',',whichOne,ii,'.time,',whichOne,ii,'.signals.values(:,3),''k'');']);
    else
        eval(['plot(',whichOne,ii,'.time,',whichOne,ii,'.signals.values(:,1),''r'',',whichOne,ii,'.time,',whichOne,ii,'.signals.values(:,2),''g'',',whichOne,ii,'.time,',whichOne,ii,'.signals.values(:,3),''k'');']);
        figure
        eval(['dataTemp=Gyro',ii,'.signals.values;']);
        eval(['timeTemp=Gyro',ii,'.time;']);
        dataTemp=dataTemp./pi.*180;
        integralRes=cumtrapz(timeTemp,dataTemp(:,:));
        %计算误差
        offsetGyroData=[timeTemp.*offset(1),timeTemp.*offset(2),timeTemp.*offset(3)];
        integralRes=integralRes+offsetGyroData;
        plot(timeTemp,integralRes(:,1),'r',timeTemp,integralRes(:,2),'g',timeTemp,integralRes(:,3),'b');
    end
    legend('x','y','z');
    title([whichOne,ii]);
end
clearvars dataTemp file i ii integralRes 