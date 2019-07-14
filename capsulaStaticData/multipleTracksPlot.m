filesDetect = dir('*.txt'); %list of the data files in parent directory
dataDetect = []; %main data matrix declaration
beamsSpacing = 874;%beams spacing in µm declaration
minTimeBetwPeaks = 0.01;%minimum time between two peaks in seconds declaration
zeroPeak2 = 0.3; %peak 2 zero value measured previously to data processing declaration
zeroTrigger = 0.6;%trigger peak value measured previously to data processing declaration
totalIntensity = zeroPeak2+zeroTrigger;%total signal level calculation (2 beams summed)
zeroPeak2Norm = zeroPeak2/totalIntensity;
zeroTriggerNorm = zeroTrigger/totalIntensity;
Nsampling = 1000;%multiplicative coefficient for getting index from time values
spheroVelo = []; %spheroid velocity array declaration
finalFile = []; %final measurement array declaration
caps = 0; %peak counter declaration
realCaps = 0;
deltaTime=[];
j=0;
minTabD3=[];
fwhmTabD3=[];
minTabD6=[];
fwhmTabD6=[];
minTabD9=[];
fwhmTabD9=[];



%DATA IMPORT PART
k=0;
for i=1:length(filesDetect)
    eval(['load ' filesDetect(i).name ' -ascii']);
    k=k+1;
end
cmap1=winter(4);
cmap2=summer(4);
cmap3=winter(4);
dataCell = cell(k);
for i=1:length(filesDetect)
    dataDetect=importdata(filesDetect(i).name);
    dataCell(i) = {dataDetect};
end

figure(1)

for i=1:length(filesDetect)
    disp(i);
    dataPlot=dataCell{i,1};
    dataYsmoothed = smooth(dataPlot(:,1),dataPlot(:,2),1);
    y = dataYsmoothed;
    y(y>1)=1;
    t = dataPlot(:,1);
    
    [pn,ln,width] = findpeaks(-y,t,'MinPeakProminence',0.2,'Annotate','extents');
    [W,INITCROSS,FINALCROSS] = pulsewidth(-y,t,'Tolerance',15);
    
    title('D3')
    %ylim([0 1])
    %xlim([-300 300])
    xlabel('Translation (µm)')
    ylabel('Transmittance (norm.)')
    plot((t-(0.5*(FINALCROSS-INITCROSS)+INITCROSS))*100,y,'-','lineWidth',1,'color',cmap1(i,:));
    hold on
    legend('0.0116 mm^2','0.0205 mm^2','0.0142 mm^2','0.0077 mm^2','Location','southoutside')
    grid on
    
end
