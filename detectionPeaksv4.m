%DetectionPeak v.4.1 - Code by Nelson Hélaine
%This code detects peaks in data collected by a capsula system, providing distribution of FWHM and absorption minimum for
%detected organoids. Thresholds of detection are tunable in variables initialization.


%data files reading and variable declaration


filesDetect = dir('*.txt'); %list of the data files in parent directory
dataDetect = []; %main data matrix declaration
beamsSpacing = 874;%beams spacing in µm declaration
minTimeBetwPeaks = 0.01;%minimum time between two peaks in seconds declaration
zeroPeak2 = 0.337; %peak 2 zero value measured previously to data processing declaration
zeroTrigger = 0.50;%trigger peak value measured previously to data processing declaration
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
%DATA IMPORT PART

for i=1:length(filesDetect)
    eval(['load ' filesDetect(i).name ' -ascii']);
end
for i=1:length(filesDetect)
    dataDetect(:,:,i) = importdata(filesDetect(i).name);
end
%PRINCIPAL LOOP ON THE DATA ARRAY (3D ARRAY)
for i=1:size(dataDetect,3) %loop on data files
    tic
    ii=0;%looping variable declaration
    disp('File #')%ID file prompting
    disp(filesDetect(i).name)
    sizeFile=size(dataDetect(: ,:,i)); %facultatif
    subDataDetect=dataDetect(:,:,i); %stockage des données du fichier dans une matrice 2D (pour plus de facilité de manip)
    subDataDetect(:,2)=subDataDetect(:,2)./max(subDataDetect(:,2)); %normalisation des données
    YsubDataDetect=smooth(subDataDetect(:,1),subDataDetect(:,2),100); %smoothing sur les données (facultatif) et stockage dans un vecteur Y
    y = YsubDataDetect(:,1);
    t = subDataDetect(:,1);% idem (stockage données X)
    %figure(i+10000)
    %findpeaks(-y,t,'MinPeakDistance',minTimeBetwPeaks,'MinPeakProminence',0.05,'Annotate','extents')             
    [pn ln width] = findpeaks(-y,t,'MinPeakDistance',minTimeBetwPeaks,'MinPeakProminence',0.05,'Annotate','extents'); % recherche des pics et leur loca avec une distance minimale de 300 ms entre chaque
    pn(pn<-0.8)=NaN; %manipulation matrice pour ne garder que les pics rentrant dans les critères de seuil
    peakData = [1+pn ln width]; %stockage dans un tableau des pics à étudier
    peakData(any(isnan(peakData), 2), :) = []; %manip matrice pour remplacer les valeurs NaN par des cases vides
    for k=1:size(peakData,1)%boucle d'exploitation des pics
        %findpeaks(-y,t,'MinPeakDistance',minTimeBetwPeaks,'MinPeakProminence',0.3,'Annotate','extents')
        if (peakData(k,1)<zeroTriggerNorm) && (size(peakData,1)>1) && (k>1) %&& (peakData(k,3)<25e-3) && (peakData(k,3)>5e-3)%si le pic est inférieur au zéro du trigger, c'est un pic à analyser et non pas un trigger
            caps = caps+1;%on augmente le compteur de capsule
            idxMinPeak2 = find(YsubDataDetect(:,1)==peakData(k,2)); %on trouve l'indice correspondant au minimum du pic à analyser
            windowDetect = 10e3*[peakData(k,2)-4*peakData(k,3) peakData(k,2)+4*peakData(k,3)];
            if (windowDetect(1)<0)
                windowDetect(1) = 1;
            end
            if (windowDetect(2)<40000) && (peakData(k,1)>zeroPeak2Norm) && (peakData(k-1,1) < peakData(k,1)) && (peakData(k,1)<zeroTriggerNorm) % condition sur la fenêtre de détection pour éviter les effets de bord et mêmes condition que pour le min
                ii = ii +1;
                deltaMin = peakData(k,2)-peakData(k-1,2);%on calcule le delta t entre le pic trigger et le pic qu'on analyse
                peak2 = subDataDetect(round(windowDetect(1)):round(windowDetect(2)),1:2); %on définit une fenêtre de données autour du minimum pour retrouver la HWHM
                peak2(:,2) = YsubDataDetect(round(windowDetect(1)):round(windowDetect(2)));
                [minMaxPeak2(:,caps)] = statelevels(peak2(:,2),1000,'mode');
                sizePeak2A = size(pulsewidth(peak2(:,2),peak2(:,1),'Polarity','Negative','Tolerance',2,'StateLevels',(minMaxPeak2(:,caps)).'),1);
                sizePeak2B = size(pulsewidth(peak2(:,2),peak2(:,1),'Polarity','Negative','Tolerance',2,'StateLevels',(minMaxPeak2(:,caps)).'),2);
                if (sizePeak2A~=0) && (sizePeak2A<2) && (sizePeak2B~=0) && (sizePeak2B<2)
                    j=j+1;
                    [W(ii,i),InitCross(ii,i),FinalCross(ii,i),MeadVal(ii,i)]=pulsewidth(peak2(:,2),peak2(:,1),'Polarity','Negative','Tolerance',2,'StateLevels',(minMaxPeak2(:,caps)).');
                    %figure((i*100)+k)
                    %pulsewidth(peak2(:,2),peak2(:,1),'Polarity','Negative','Tolerance',2,'StateLevels',(minMaxPeak2(:,caps)).');
                    timePeak2 = InitCross(ii,i)+((FinalCross(ii,i)-InitCross(ii,i))/2);
                    timePeakTrigger = peakData(k-1,2);
                    deltaTime(j) = timePeak2-timePeakTrigger;                   
                    spheroVelo(caps,i) = beamsSpacing/deltaTime(j); %on calcule la vitesse de la capsule analysée
                    if spheroVelo(caps,i)>1 && ((spheroVelo(caps,i))<250000) && (W(ii,i)*spheroVelo(caps,i)<600)
                        minPeak2(caps,i)=minMaxPeak2(1,caps); % on le stocke
                        HWHM(caps,i)=W(ii,i)*spheroVelo(caps,i);
                        figure((i*100)+k)
                        findpeaks(-y,t,'MinPeakDistance',minTimeBetwPeaks,'MinPeakProminence',0.05,'Annotate','extents')
                        realCaps =  realCaps+1;
                    else
                        minPeak2(caps,i)=NaN; % on le stocke
                        HWHM(caps,i)=NaN;
                        spheroVelo(caps,i)=NaN;
                    end         
                else
                    HWHM(caps,i)=NaN;
                    minPeak2(caps,i)=NaN;
                    spheroVelo(caps,i) = NaN; %on calcule la vitesse de la capsule analysée
                end
                windowDetect = [];
            else
                HWHM(caps,i)=NaN;
                minPeak2(caps,i) = NaN; %sinon on ne stocke que du vide
                spheroVelo(caps,i) = NaN;
            end
            minPeak2(minPeak2==0)=NaN;
            HWHM(HWHM==0)=NaN;
            spheroVelo(spheroVelo==0)=NaN;
            %finalFile = [finalFile;spheroVelo(caps,i),minPeak2(caps,i),HWHM(caps,i)]; %on stocke tout dans un fichier de données final
        end      
    end
    peakData = [];
    filename = ['datatest', num2str(i)];  
    savefig(filename);
    toc  
end
disp(finalFile);
binCount = round(1+3.22*log(realCaps));%Sturge's rule for bin count
histoFilesName = ['Chip T ,150/450 µl/mn, 2 identical beams - 874 µm separation - N = ',num2str(realCaps)];
%on trace les distributions
figure()
m=spheroVelo;
m = (m(~isnan(m)));
histogram(m,binCount);
grid on
ax = gca;
ax.FontSize = 14; 
title(histoFilesName,'FontSize',14)
xlabel('Spheroids velocity (µm/s)','FontSize',14)
ylabel('Counts','FontSize',14);
figure()
m=minPeak2;
m = (m(~isnan(m)));
histogram(m,binCount);
grid on
ax = gca;
ax.FontSize = 14;
title(histoFilesName,'FontSize',14)
xlabel('Minimum voltage detected (V)','FontSize',14);
ylabel('Counts','FontSize',14);
figure()
m=HWHM;
m = (m(~isnan(m)));
histogram(m,binCount);
grid on
ax = gca;
ax.FontSize = 14;
title(histoFilesName,'FontSize',14)
xlabel('HWHM (µm)','FontSize',14)
ylabel('Counts','FontSize',14);
