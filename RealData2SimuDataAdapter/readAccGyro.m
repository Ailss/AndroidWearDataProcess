tic
for i=3:1:length(file)
    fid=fopen([path,file(i).name]);
    if(fid<=0)
        error(['error: The file named "',file(i).name,'" was not open!']);
    end
    DS = textscan(fid,'%s %n %n %n %n %n');
    fclose(fid);
    GyroTemp=[];
    AccTemp=[];
    for j=1:1:length(DS{1})
        if(strcmp(DS{1}(j),'gyro'))
            GyroTemp=[GyroTemp;(DS{2}(j)-DS{2}(1))/1000000000,DS{3}(j),DS{4}(j),DS{5}(j),norm([DS{3}(j),DS{4}(j),DS{5}(j)])];
        %加速度都要全反
        elseif (strcmp(DS{1}(j),'acc'))
            AccTemp=[AccTemp;(DS{2}(j)-DS{2}(1))/1000000000,-DS{3}(j),-DS{4}(j),-DS{5}(j),norm([-DS{3}(j),-DS{4}(j),-DS{5}(j)])];
        end
    end
    eval(['Gyro',num2str(i-2),'=GyroTemp;']);
    eval(['Acc',num2str(i-2),'=AccTemp;']);
end
%去掉前3秒的准备时间
for i=1:1:length(file)-2
    ii=num2str(i);
    eval(['Acc',ii,'=Acc',ii,'(find(Acc',ii,'(:,1)>3),:);']);
    eval(['Acc',ii,'(:,1)=Acc',ii,'(:,1)-3;']);
    eval(['Gyro',ii,'=Gyro',ii,'(find(Gyro',ii,'(:,1)>3),:);']);
    eval(['Gyro',ii,'(:,1)=Gyro',ii,'(:,1)-3;']);
end
toc
clearvars i ii fid DS ans GyroTemp AccTemp j