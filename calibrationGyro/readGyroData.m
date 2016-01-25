tic
for i=3:1:length(file)
    fid=fopen([path,file(i).name]);
    if(fid<=0)
        error(['error: The file named "',file(i).name,'" was not open!']);
    end
    DS = textscan(fid,'%s %n %n %n %n %n');
    fclose(fid);
    GyroTemp=[];
    for j=1:1:length(DS{1})
        if(strcmp(DS{1}(j),'gyro'))
            GyroTemp=[GyroTemp;(DS{2}(j)-DS{2}(1))/1000000000,DS{3}(j),DS{4}(j),DS{5}(j),norm([DS{3}(j),DS{4}(j),DS{5}(j)])];
        end
    end
    eval(['Gyro',num2str(i-2),'=GyroTemp;']);
end
toc