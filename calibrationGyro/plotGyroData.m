whichOne=1;
w=num2str(whichOne);
figure
eval(['plot(Gyro',w,'(:,1),Gyro',w,'(:,2),''r'',Gyro',w,'(:,1),Gyro',w,'(:,3),''g'',Gyro',w,'(:,1),Gyro',w,'(:,4),''k'');'])
legend('red    --X','green--Y','balck --Z');