% string=['acc','mag','earth','post','gamepost','gyro'];
% for n=1:1:6
%     eval(['for i=2:1:length(',string(n),'TimeLine)']);
%     eval(['if(',string(n),'TimeLine(i)-',string(n),'TimeLine(i-1))>1e9']);
%     eval('tempPos=i;');
%     eval('break;');
%     eval('end;');
%     eval('end;');
%     eval([string(n),'TimeLine(tempPos:length(',string(n),'TimeLine))=',string(n),'TimeLine(tempPos:length(',string(n),'TimeLine))-temp;']);
% end


for i=2:1:length(accTimeLine)
    if(accTimeLine(i)-accTimeLine(i-1))>1e9
        tempPos=i;
        break;
    end
end
temp=accTimeLine(tempPos)-accTimeLine(tempPos-1)-(accTimeLine(tempPos-1)-accTimeLine(tempPos-2));
accTimeLine(tempPos:length(accTimeLine))=accTimeLine(tempPos:length(accTimeLine))-temp;
for i=2:1:length(magTimeLine)
    if(magTimeLine(i)-magTimeLine(i-1))>1e9
        tempPos=i;
        break;
    end
end
temp=magTimeLine(tempPos)-magTimeLine(tempPos-1)-(magTimeLine(tempPos-1)-magTimeLine(tempPos-2));
magTimeLine(tempPos:length(magTimeLine))=magTimeLine(tempPos:length(magTimeLine))-temp;
for i=2:1:length(gyroTimeLine)
    if(gyroTimeLine(i)-gyroTimeLine(i-1))>1e9
        tempPos=i;
        break;
    end
end
temp=gyroTimeLine(tempPos)-gyroTimeLine(tempPos-1)-(gyroTimeLine(tempPos-1)-gyroTimeLine(tempPos-2));
gyroTimeLine(tempPos:length(gyroTimeLine))=gyroTimeLine(tempPos:length(gyroTimeLine))-temp;
for i=2:1:length(earthTimeLine)
    if(earthTimeLine(i)-earthTimeLine(i-1))>1e9
        tempPos=i;
        break;
    end
end
temp=earthTimeLine(tempPos)-earthTimeLine(tempPos-1)-(earthTimeLine(tempPos-1)-earthTimeLine(tempPos-2));
earthTimeLine(tempPos:length(earthTimeLine))=earthTimeLine(tempPos:length(earthTimeLine))-temp;
for i=2:1:length(oriTimeLine)
    if(oriTimeLine(i)-oriTimeLine(i-1))>1e9
        tempPos=i;
        break;
    end
end
temp=oriTimeLine(tempPos)-oriTimeLine(tempPos-1)-(oriTimeLine(tempPos-1)-oriTimeLine(tempPos-2));
oriTimeLine(tempPos:length(oriTimeLine))=oriTimeLine(tempPos:length(oriTimeLine))-temp;
for i=2:1:length(quatTimeLine)
    if(quatTimeLine(i)-quatTimeLine(i-1))>1e9
        tempPos=i;
        break;
    end
end
temp=quatTimeLine(tempPos)-quatTimeLine(tempPos-1)-(quatTimeLine(tempPos-1)-quatTimeLine(tempPos-2));
quatTimeLine(tempPos:length(quatTimeLine))=quatTimeLine(tempPos:length(quatTimeLine))-temp;
for i=2:1:length(gquatTimeLine)
    if(gquatTimeLine(i)-gquatTimeLine(i-1))>1e9
        tempPos=i;
        break;
    end
end
temp=gquatTimeLine(tempPos)-gquatTimeLine(tempPos-1)-(gquatTimeLine(tempPos-1)-gquatTimeLine(tempPos-2));
gquatTimeLine(tempPos:length(gquatTimeLine))=gquatTimeLine(tempPos:length(gquatTimeLine))-temp;

