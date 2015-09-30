accTimeLine2=(accTimeLine-accTimeLine(1))/1e9;
magTimeLine2=(magTimeLine-magTimeLine(1))/1e9;
oriTimeLine2=(oriTimeLine-oriTimeLine(1))/1e9;
gyroTimeLine2=(gyroTimeLine-gyroTimeLine(1))/1e9;
earthTimeLine2=(earthTimeLine-earthTimeLine(1))/1e9;
postTimeLine2=(postTimeLine-postTimeLine(1))/1e9;
gamepostTimeLine2=(gamepostTimeLine-gamepostTimeLine(1))/1e9;
% % 测试
% executeWhichPara='post';
% hold on
% eval(['plot(',executeWhichPara,'TimeLine2,',executeWhichPara,'(:,1))']);
% eval(['plot(',executeWhichPara,'TimeLine2,',executeWhichPara,'(:,2))']);
% eval(['plot(',executeWhichPara,'TimeLine2,',executeWhichPara,'(:,3))']);
% grid
% plot(ones(1000,1).*3,-6:12/999:6);
% plot(ones(1000,1).*8,-6:12/999:6);
% hold off

% 前3秒是用来标定人的朝向的
timetemp1=find(postTimeLine2>1.5);
timetemp2=find(postTimeLine2<2.5);
timetemp=intersect(timetemp1,timetemp2);
% time=postTimeLine2(timetemp);
posttemp=post(timetemp,:);
% hold on
% plot(time,posttemp(:,1));
% plot(time,posttemp(:,2));
% plot(time,posttemp(:,3));
% hold off
postStart=mean(posttemp);
postStart=postStart-[pi/2,0,0];

% % 8秒后为正式采集到的数据
% timetemp=find(postTimeLine2>8);
% time=postTimeLine2(timetemp);
clear *2;
clear timetemp timetemp1 posttemp;