tic

accOrigin=[];
%1. 将所有文件的acc合并
for i=3:1:length(file)
% i=3;
    fid=fopen([path,file(i).name]);
    if(fid<=0)
        error(['error: The file named "',file(i).name,'" was not open!']);
    end
    DS = textscan(fid,'%s %n %n %n %n %n');
    fclose(fid);
    accTemp=[];
    for j=1:1:length(DS{1})
        if(strcmp(DS{1}(j),'acc'))
            accTemp=[accTemp;DS{3}(j),DS{4}(j),DS{5}(j),norm([DS{3}(j),DS{4}(j),DS{5}(j)])];
        end
    end
    accOrigin=[accOrigin;accTemp];
end

%2. 采样
sampling=1:samplingRate:length(accOrigin);
acc=accOrigin(sampling,:);



clearvars -except acc
save acc.mat
toc