close all
if strcmp(whichOne,'Acc')
    eval(['plot(data.time',whichOne,',data.data',whichOne,'(:,1),''r'',data.time',whichOne,',data.data',whichOne,'(:,2),''g'',data.time',whichOne,',data.data',whichOne,'(:,3),''k'');']);
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
        eval(['plot(',whichOne,ii,'(:,1),',whichOne,ii,'(:,2),''r'',',whichOne,ii,'(:,1),',whichOne,ii,'(:,3),''g'',',whichOne,ii,'(:,1),',whichOne,ii,'(:,4),''k'');']);
    else
        eval(['plot(',whichOne,ii,'(:,1),',whichOne,ii,'(:,2),''r'',',whichOne,ii,'(:,1),',whichOne,ii,'(:,3),''g'',',whichOne,ii,'(:,1),',whichOne,ii,'(:,4),''k'');']);
        figure
        eval(['dataTemp=Gyro',ii,'(:,1:1:4);']);
        dataTemp(:,2:1:4)=dataTemp(:,2:1:4)/pi*180;
        integralRes=cumtrapz(dataTemp(:,1),dataTemp(:,2:1:4));
        %计算误差
        offsetGyroData=[dataTemp(:,1).*offset(1),dataTemp(:,1).*offset(2),dataTemp(:,1).*offset(3)];
        integralRes=integralRes+offsetGyroData;
        plot(dataTemp(:,1),integralRes(:,1),'r',dataTemp(:,1),integralRes(:,2),'g',dataTemp(:,1),integralRes(:,3),'b');
    end
    legend('x','y','z');
    title([whichOne,ii]);
end
clearvars dataTemp file i ii integralRes 