
tic;

firstSensorEventTime = min([accTimeLine(1) magTimeLine(1) gyroTimeLine(1) oriTimeLine(1) postTimeLine(1) earthTimeLine(1) gamepostTimeLine(1)]);

[accTimeLine,timeLineIndex] = sort(accTimeLine);
acc(:,:) = acc(timeLineIndex,:);
[gyroTimeLine,timeLineIndex] = sort(gyroTimeLine);
gyro(:,:) = gyro(timeLineIndex,:);
[magTimeLine,timeLineIndex] = sort(magTimeLine);
mag(:,:) = mag(timeLineIndex,:);
[oriTimeLine,timeLineIndex] = sort(oriTimeLine);
ori(:,:) = ori(timeLineIndex,:);
[earthTimeLine,timeLineIndex] = sort(earthTimeLine);
earth(:,:) = earth(timeLineIndex,:);
[postTimeLine,timeLineIndex] = sort(postTimeLine);
post(:,:) = post(timeLineIndex,:);
[gamepostTimeLine,timeLineIndex] = sort(gamepostTimeLine);
gamepost(:,:) = gamepost(timeLineIndex,:);

% for i = 2:length(gyro)-1
%     if abs(gyro(i,1))>0.06 && abs(gyro(i,1))<0.09 && abs(gyro(i-1,1))<0.02 && abs(gyro(i+1,1))<0.02
%         gyro(i,1) = gyro(i-1,1) + (gyroTimeLine(i)-gyroTimeLine(i-1))/(gyroTimeLine(i)-gyroTimeLine(i+1))*(gyro(i+1,1)-gyro(i-1,1));
%     end
%     if abs(gyro(i,3))>0.06 && abs(gyro(i,3))<0.09 && abs(gyro(i-1,3))<0.02 && abs(gyro(i+1,3))<0.02
%         gyro(i,3) = gyro(i-1,3) + (gyroTimeLine(i)-gyroTimeLine(i-1))/(gyroTimeLine(i)-gyroTimeLine(i+1))*(gyro(i+1,3)-gyro(i-1,3));
%     end
% end


accxts = timeseries(acc(:,1),accTimeLine,'name','accx');
accxts.TimeInfo.Units='nanoseconds';
accyts = timeseries(acc(:,2),accTimeLine,'name','accy');
accyts.TimeInfo.Units='nanoseconds';
acczts = timeseries(acc(:,3),accTimeLine,'name','accz');
acczts.TimeInfo.Units='nanoseconds';
accTC = tscollection({accxts accyts acczts});

magxts = timeseries(mag(:,1),magTimeLine,'name','magx');
magxts.TimeInfo.Units='nanoseconds';
magyts = timeseries(mag(:,2),magTimeLine,'name','magy');
magyts.TimeInfo.Units='nanoseconds';
magzts = timeseries(mag(:,3),magTimeLine,'name','magz');
magzts.TimeInfo.Units='nanoseconds';
magTC = tscollection({magxts magyts magzts});

gyroxts = timeseries(gyro(:,1),gyroTimeLine,'name','gyrox');
gyroxts.TimeInfo.Units='nanoseconds';
gyroyts = timeseries(gyro(:,2),gyroTimeLine,'name','gyroy');
gyroyts.TimeInfo.Units='nanoseconds';
gyrozts = timeseries(gyro(:,3),gyroTimeLine,'name','gyroz');
gyrozts.TimeInfo.Units='nanoseconds';
gyroTC = tscollection({gyroxts gyroyts gyrozts});

orixts = timeseries(ori(:,1),oriTimeLine,'name','orix');
orixts.TimeInfo.Units='nanoseconds';
oriyts = timeseries(ori(:,2),oriTimeLine,'name','oriy');
oriyts.TimeInfo.Units='nanoseconds';
orizts = timeseries(ori(:,3),oriTimeLine,'name','oriz');
orizts.TimeInfo.Units='nanoseconds';
oriTC = tscollection({orixts oriyts orizts});

earthxts = timeseries(earth(:,1),earthTimeLine,'name','earthx');
earthxts.TimeInfo.Units='nanoseconds';
earthyts = timeseries(earth(:,2),earthTimeLine,'name','earthy');
earthyts.TimeInfo.Units='nanoseconds';
earthzts = timeseries(earth(:,3),earthTimeLine,'name','earthz');
earthzts.TimeInfo.Units='nanoseconds';
earthTC = tscollection({earthxts earthyts earthzts});

postxts = timeseries(post(:,1),postTimeLine,'name','postx');
postxts.TimeInfo.Units='nanoseconds';
postyts = timeseries(post(:,2),postTimeLine,'name','posty');
postyts.TimeInfo.Units='nanoseconds';
postzts = timeseries(post(:,3),postTimeLine,'name','postz');
postzts.TimeInfo.Units='nanoseconds';
postTC = tscollection({postxts postyts postzts});

gamepostxts = timeseries(gamepost(:,1),gamepostTimeLine,'name','gamepostx');
gamepostxts.TimeInfo.Units='nanoseconds';
gamepostyts = timeseries(gamepost(:,2),gamepostTimeLine,'name','gameposty');
gamepostyts.TimeInfo.Units='nanoseconds';
gamepostzts = timeseries(gamepost(:,3),gamepostTimeLine,'name','gamepostz');
gamepostzts.TimeInfo.Units='nanoseconds';
gamepostTC = tscollection({gamepostxts gamepostyts gamepostzts});

startTime = [accTimeLine(1) gyroTimeLine(1) magTimeLine(1) oriTimeLine(1) earthTimeLine(1) postTimeLine(1) gamepostTimeLine(1)];
startTime  = sort(startTime);
startTime = startTime(end);

endTime = [accTimeLine(end) gyroTimeLine(end) magTimeLine(end) oriTimeLine(end) earthTimeLine(end) postTimeLine(end) gamepostTimeLine(end)];
endTime  = sort(endTime);
endTime = endTime(1);

%  200Hz resample.
accTC = resample(accTC,startTime:5*10e5:endTime);
gyroTC = resample(gyroTC,startTime:5*10e5:endTime,'linear');
magTC = resample(magTC,startTime:5*10e5:endTime);
oriTC = resample(oriTC,startTime:5*10e5:endTime);
earthTC = resample(earthTC,startTime:5*10e5:endTime);
postTC = resample(postTC,startTime:5*10e5:endTime);
gamepostTC = resample(gamepostTC,startTime:5*10e5:endTime);


%  back to normal data series
acc = [accTC.accx.Data accTC.accy.Data accTC.accz.Data];
gyro = [gyroTC.gyrox.Data gyroTC.gyroy.Data gyroTC.gyroz.Data];
mag = [magTC.magx.Data magTC.magy.Data magTC.magz.Data];
earth = [earthTC.earthx.Data earthTC.earthy.Data earthTC.earthz.Data];
ori = [oriTC.orix.Data oriTC.oriy.Data oriTC.oriz.Data];
post = [postTC.postx.Data postTC.posty.Data postTC.postz.Data];
gamepost = [gamepostTC.gamepostx.Data gamepostTC.gameposty.Data gamepostTC.gamepostz.Data];

cutStart = 50; cutEnd = length(acc);
acc = acc(cutStart:cutEnd,:);
gyro = gyro(cutStart:cutEnd,:);
mag = mag(cutStart:cutEnd,:);
ori = ori(cutStart:cutEnd,:);
earth=earth(cutStart:cutEnd,:);
post=post(cutStart:cutEnd,:);
gamepost=gamepost(cutStart:cutEnd,:);
timeLine = accTC.Time;
timeLine = timeLine(cutStart:cutEnd);

timeLine2 = (timeLine - firstSensorEventTime)/1e9;
% timeLine2 = timeLine2 + sensorStartTime;

for i = 1:length(acc)
acc(i,4) = norm(acc(i,1:3));
% mag(i,4) = norm(mag(i,1:3));
end


disp(['Data Alignment Done. ']);
toc;

clear i;
clear *ts;
clear endTime startTime ;
clear cut*;
clear *TC *TimeLine *Index;
