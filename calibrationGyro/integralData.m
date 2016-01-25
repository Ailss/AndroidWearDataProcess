whichOne=1;
w=num2str(whichOne);
if whichAxes=='X'
    whichAxesStr='2';
    color='r';
else
    if whichAxes=='Y'
        whichAxesStr='3';
        color='g';
    else
        whichAxesStr='4';
        color='k';
    end
end
dataAxes=[Gyro1(:,1),Gyro1(:,2)];
figure
eval(['dataAxes=[Gyro',w,'(:,1),Gyro',w,'(:,',whichAxesStr,')];']);
plot(dataAxes(:,1),dataAxes(:,2),color);
dataAxesAngle=dataAxes(:,2)/pi*180;
figure
plot(dataAxes(:,1),dataAxesAngle,color);
integralRes=cumtrapz(dataAxes(:,1),dataAxesAngle);
figure
plot(dataAxes(:,1),integralRes);