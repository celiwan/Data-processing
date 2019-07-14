filesDetect1 = dir('tracks/*.txt'); %list of the data files in parent directory
filesDetect2 = dir('tracksR/*.txt'); %list of the data files in parent directory
dataDetect1 = []; %main data matrix declaration
meanVel1 = [];
errorMeanVel1 = [];
dataDetect2 = []; %main data matrix declaration
meanVel2 = [];
errorMeanVel2 = [];

flowRates = [60:20:160];
%DATA IMPORT PART
k=0;
for i=1:length(filesDetect1)
    %eval(['load ' filesDetect(i).name ' -ascii']);
    k=k+1;
end

dataCell1 = cell(k);
cd tracks
for i=1:length(filesDetect1)
    dataDetect1=importdata(filesDetect1(i).name);
    dataCell1(i) = {dataDetect1};
end
cd ../tracksR

k=0;
for i=1:length(filesDetect2)
    %eval(['load ' filesDetect(i).name ' -ascii']);
    k=k+1;
end

dataCell2 = cell(k);
for i=1:length(filesDetect2)
    dataDetect2=importdata(filesDetect2(i).name);
    dataCell2(i) = {dataDetect2};
end

%COLORMAPS AND PLOTTING TOOLS
cmap1=winter(4);
cmap2=summer(4);
cmap3=winter(4);

for i=7:length(filesDetect1)
    meanVel1 = [meanVel1 mean(dataCell1{i,1}.data(:,17))];
    errorMeanVel1 = [errorMeanVel1 ...
        0.5*(mean(dataCell1{i,1}.data(:,21)))];
        %/sqrt(length(dataCell{i,1}.data(:,17))))];
    %disp(min(dataCell1{i,1}.data(:,17)));
end

for i=7:length(filesDetect2)
    meanVel2 = [meanVel2 mean(dataCell2{i,1}.data(:,17))];
    errorMeanVel2 = [errorMeanVel2 ...
        0.5*(mean(dataCell2{i,1}.data(:,21)))];
        %/sqrt(length(dataCell{i,1}.data(:,17))))];
    %disp(min(dataCell1{i,1}.data(:,17)));
end

figure(1)


errorbar(flowRates,meanVel1,errorMeanVel1,'-o','lineWidth',1,'color','red',...
    'DisplayName','Injection zone')
hold on
errorbar(flowRates,meanVel2,errorMeanVel2,'-o','lineWidth',1,'color','blue',...
    'DisplayName','Reduction zone')

ylabel('Velocity (mm/sec)')
grid on
ylim([0 30])

xlabel('Flow rate (µl/min)')
title('Flowing capsules in rectangular channel')
ylim([0 30])
grid on
ylabel('Mean velocity (mm/sec)')
