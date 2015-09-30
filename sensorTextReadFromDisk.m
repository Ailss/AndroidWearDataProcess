% 原始数据前3秒是标定位置
% 接下来5秒是过度时间，5秒结束开始正式数据的采集
format long g;
digits(64);
tic;
format long;

fid = fopen(filepath,'rt');
DS = textscan(fid,'%s %n %n %n %n %n');
fclose(fid);
accLength = 0;
magLength = 0;
gyroLength =0;
oriLength = 0;
quatLength = 0;
earthLength = 0;
gquatLength = 0;
other = 0;

% magStartTime=DS{2}(1);
% accStartTime=DS{2}(2);
% oriStartTime=DS{2}(3);
% quatStartTime=DS{2}(4);
% gquatStartTime=DS{2}(5);
% earthStartTime=DS{2}(6);
% gyroStartTime=DS{2}(7);

for i = 1:length(DS{1})
    type = DS{1}(i);
    switch type{1}(1:2)
        case 'ac'
            accLength = accLength +1;
        case 'gy'
            gyroLength = gyroLength +1;
        case 'ma'
            magLength = magLength+1;
        case 'or'
            oriLength = oriLength+1;
        case 'qu'
            quatLength = quatLength+1;
        case 'ea'
            earthLength = earthLength+1;
        case 'gq'
            gquatLength = gquatLength+1;
%         case 't'
%             timeString= type{1}(5:end);
%             sensorStartTime=str2num(timeString)/1000;
        otherwise
            other = other+1;
    end
end

acc = zeros(accLength,4);
ori = zeros(oriLength,3);
mag = zeros(magLength,4);
gyro = zeros(gyroLength,3);
quat = zeros(quatLength,4);
gquat = zeros(gquatLength,4);
earth = zeros(earthLength,3);

gyroTimeLine= zeros(gyroLength,1);
accTimeLine = zeros(accLength,1);
oriTimeLine = zeros(oriLength,1);
magTimeLine = zeros(magLength,1);
quatTimeLine = zeros(quatLength,1);
earthTimeLine = zeros(earthLength,1);
gquatTimeLine = zeros(gquatLength,1);

accCount = 1;magCount = 1; gyroCount = 1; oriCount = 1; earthCount = 1;
quatCount = 1;
gquatCount = 1;

matrix = cell2mat(DS(2:6));
for i = 1:length(matrix)
    type = DS{1}(i);
    time = matrix(i,1);
    
        switch type{1}(1:2)
        case 'ac'
            data = matrix(i,2:4);
%             time=time-accStartTime;
            all = norm(data);
            acc(accCount,:) = [data all];
            accTimeLine(accCount)=time;
            accCount = accCount+1;
        case 'gy'
%             time=time-gyroStartTime;
            data = matrix(i,2:4);
            gyro(gyroCount,:) = data;
            gyroTimeLine(gyroCount)=time;
            gyroCount = gyroCount +1;
        case 'ma'
%             time=time-magStartTime;
            data = matrix(i,2:4);
            all = norm(data);
            mag(magCount,:) = [data all];
            magTimeLine(magCount)=time;
            magCount = magCount+1;
        case 'or'
%             time=time-oriStartTime;
            data = matrix(i,2:4);
            ori(oriCount,:) = data;
            oriTimeLine(oriCount)=time;
            oriCount = oriCount+1;           
        case 'qu'
%             time=time-quatStartTime;
            data = matrix(i,2:5);
            quat(quatCount,:) = data;
            quatTimeLine(quatCount)=time;
            quatCount = quatCount+1; 
        case 'ea'
%             time=time-earthStartTime;
            data = matrix(i,2:4);
            earth(earthCount,:) = data;
            earthTimeLine(earthCount)=time;
            earthCount = earthCount+1; 
        case 'gq'
%             time=time-gquatStartTime;
            data = matrix(i,2:5);
            gquat(gquatCount,:) = data;
            gquatTimeLine(gquatCount)=time;
            gquatCount = gquatCount+1; 
        end
end

rawData.acc = acc;
rawData.gyro = gyro;
rawData.mag = mag;
rawData.ori = ori;
rawData.quat = quat;
rawData.earth = earth;
rawData.gquat = gquat;
%%
%处理漂移
dealwithdrift;
%%
rawData.accTimeLine = accTimeLine;
rawData.oriTimeLine = oriTimeLine;
rawData.magTimeLine = magTimeLine;
rawData.gyroTimeLine = gyroTimeLine;
rawData.quatTimeLine = quatTimeLine;
rawData.gquatTimeLine = gquatTimeLine;
rawData.earthTimeLine = earthTimeLine;
% rawData.sensorStartTime = sensorStartTime;
[postx,posty,postz]=quat2angle(quat);
post=[postx,posty,postz];
postTimeLine=quatTimeLine;
[gamepostx,gameposty,gamepostz]=quat2angle(gquat);
gamepost=[gamepostx,gameposty,gamepostz];
gamepostTimeLine=gquatTimeLine;


clear *Length *Count;
clear *x *y *z;
clearvars -except sensorStart* gyro* acc* mag* ori* earth* post* gamepost* operationTime rawData filename;

old.acc=acc;
old.earth=earth;
old.post=post;
old.gamepost=gamepost;
old.gyro=gyro;
old.ori=ori;
old.mag=mag;


clock;
% fileName = ['zzsensor @' num2str(ans(2)) '.' num2str(ans(3)) '.' num2str(ans(4)) '.' num2str(ans(5)) '.mat'];

% save(fileName);
% movefile(fileName,'./data/','f');
clear fileName;
clear ans;

disp(['Sensor.txt Read done.']);
toc;


