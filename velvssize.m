cd Links
filesDetect=dir('*.csv');
Link=struct;
%dataDetect = struct;
dataCell = {};
k=0;
cd ..
radiuses =importdata('diameters.xlsx','\t',1); 
cd Links
cmap = jet(8);
radTab = [7*10^-5:0.1*10^-5:11*10^-5];
dynaVisc = 10^-3;
g = 9.81;
deau = 1000;
rhopart = 1.0311*10^3;
a = [];

for i=1:length(filesDetect)
    k=k+1;
    if contains(filesDetect(i).name,'._Link')
        disp('o');
    else
      dataDetect(i).data=importdata(filesDetect(i).name,';',1);
      a = importdata(filesDetect(i).name,';',1);
      disp(a);
      dataCell(i) = {dataDetect(i).data};
      disp(dataCell(i));
    end

end
meanVel = [];
errorVel = [];
figure(1)
for i=9:length(filesDetect)
    temp = [];
    temp=dataCell{1,i}.data(:,9);
    disp(length(temp))
    k=0;
    temp2=temp;
    for ii=1:length(temp)
        k=k+1;
        if temp(ii)==0 && ii<round(0.3*length(temp)) && ii >round(0.6*length(temp))
            temp2(k,:)=[];
            k=k-1;
        end
        
    end
    temp(temp==0)=NaN;
    %disp(temp);
    meanVel = [meanVel mean(temp2)];
    errorVel = [errorVel 0.5*(abs(max(temp2)-mean(temp2))+...
        abs(min(temp2)-mean(temp2)))./sqrt(length(temp2))]
    %disp('file')
    %disp(temp)
    %disp('mean');
    %disp(mean(temp2));
    plot(dataCell{1,i}.data(:,7),smooth(dataCell{1,i}.data(:,7),...
        temp,10),'.-','color',cmap(i-8,:),'DisplayName',...
        ['R = ' num2str(radiuses.data(i-8,3)) ' µm'])
    
    hold on
end
finalTab = [meanVel' radiuses.data(:,3) errorVel'];
%disp(finalTab)
finalTab = sortrows(finalTab,2);
finalTab(2,:)=[];

veloSph = 700e-6;
rsph = 60e-6;
dpart = [];

for i=1:length(finalTab)
    dpart(i) = (4.5*dynaVisc*(finalTab(i,1).*10^-3)/(g*(finalTab(i,2).*10^-6).^2))+deau;
end
dpart = dpart';
dpart(2,:)=[];
meandPart = mean(dpart);


%disp(finalTab)
xlim([0.25 2])
ylim([0 2])
grid on
ylabel('Velocity (m/s)')
bar = [finalTab(:,2)];
%finalTab(4,:)=[];
disp(bar);
%colorbar('TicksMode','manual',);
figure(2)
plot(radTab(1,:),(2/(9*dynaVisc))*(rhopart-deau)*9.81*((radTab(1,:)).^2))
hold on

errorbar(finalTab(:,2).*10^-6,finalTab(:,1).*10^-3,finalTab(:,3).*10^-3,'o')
xlabel('Spheroid radius (m)','Interpreter','latex')
ylabel('mean velocity $(m.s^{-1})$','Interpreter','latex')
grid on
leg1 = legend('$v = \frac{2r^2\Delta\rho g}{9\eta}$','Experimental points');
set(leg1,'Interpreter','latex');



%figure(3)

%plot(finalTab(:,1),dpart(:))
    
