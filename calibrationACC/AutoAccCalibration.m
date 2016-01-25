tic
load acc.mat
% close all;
% clearvars -except old;
% acc=old.acc;
% clear old;
% 
% accfilt= filtfilt(ones(1,20)/20,1,acc(:,1:3));
% accdiff = diff(accfilt);
% accdiff = abs(accdiff)*ones(3,1);
% [pks,locs] = findpeaks(accdiff,'MINPEAKHEIGHT',0.1,'MINPEAKDISTANCE',50);
% keyPoints = [];
% for i = 1:length(pks)-1
%     [c,index] = min(accdiff(locs(i):locs(i+1)));
%     if c<=0.005
%         keyPoints = [ keyPoints index+locs(i)];
%     end
% end

% keyPointsData = acc(keyPoints,1:3);
temp=norm(acc);
keyPoints=1:1:length(acc);
keyPointsData=acc;


% keyPointsData(:,1) = keyPointsData(:,1)./keyPointsData(:,4);
% keyPointsData(:,2) = keyPointsData(:,2)./keyPointsData(:,4);
% keyPointsData(:,3) = keyPointsData(:,3)./keyPointsData(:,4);

% gcf = figure;
% axes1 = axes('Parent',gcf,'PlotBoxAspectRatio',[1 1 1],...
%     'DataAspectRatio',[1 1 1]);
% hold on;
% scatter3(keyPointsData(:,1),keyPointsData(:,2),keyPointsData(:,3),'filled');
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% for i = 1:length(keyPointsData)
%     text(keyPointsData(i,1)+0.3,keyPointsData(i,2)+0.3,keyPointsData(i,3)+0.3,num2str(i));
% end
% hold off;

figure
hold on;
plot(acc);
plot(keyPoints,acc(keyPoints,4),'.','MarkerSize',20);
hold off;

E = 9.8*ones(length(keyPointsData),1);

x0 = [1 1 1 0 0 0 0 0 0];
[x , resnorm] = lsqcurvefit(@accelerometerError,x0,keyPointsData,E);

scaleMatrix = [x(1) x(4) x(5);
               x(4) x(2) x(6);
               x(5) x(6) x(3)];
offsetVector = [x(7) x(8) x(9)];
acc(:,1) = acc(:,1) - x(7);
acc(:,2) = acc(:,2) - x(8);
acc(:,3) = acc(:,3) - x(9);
acc = scaleMatrix * acc(:,1:3)';
acc = acc';
acc(:,4) = sqrt(sum((acc.^2)')');
save CalibrationData.mat x;
clearvars -except x
toc