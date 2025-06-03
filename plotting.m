clear;
addpath(fullfile(pwd,'utils'));
subfolder='test';
folderpath=fullfile(pwd,'output',subfolder);
time_index_max=size(dir(fullfile(folderpath,'*.mat')),1)-1;

iTP = 2;

%load input
filename='input.xlsx';
    data=xlsread(filename,'data');
    input_solid=readtable(filename).solid';
    input_liquid=readtable(filename).liquid';
    N=readtable(filename).N';
    Kd_cpx=readtable(filename).Kd_Cpx';

TE_list = split('La Ce Pr Nd Sm Eu Gd Tb Dy Ho Er Tm Yb Lu'); nTE = length(TE_list);
TE_list_full = split('La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu');

%normalization
solid_N_aux=input_solid./N;
liquid_N_aux=input_liquid./N;
eq_solid_N_aux=Kd_cpx.*input_liquid./N;
data_N_aux=data./N;

%Pm interpolation
Pm_pos=find(contains(TE_list_full,"Pm"));
solid_N=[solid_N_aux(1:Pm_pos-1) sqrt(solid_N_aux(Pm_pos-1)*solid_N_aux(Pm_pos)) solid_N_aux(Pm_pos:end)] ; 
liquid_N=[liquid_N_aux(1:Pm_pos-1) sqrt(liquid_N_aux(Pm_pos-1)*liquid_N_aux(Pm_pos)) liquid_N_aux(Pm_pos:end)] ;
eq_solid_N=[eq_solid_N_aux(1:Pm_pos-1) sqrt(eq_solid_N_aux(Pm_pos-1)*eq_solid_N_aux(Pm_pos)) eq_solid_N_aux(Pm_pos:end)] ;
data_N=[data_N_aux(:,1:Pm_pos-1) sqrt(data_N_aux(:,Pm_pos-1).*data_N_aux(:,Pm_pos)) data_N_aux(:,Pm_pos:end)] ;

%Eu anomaly
Sm_pos=find(contains(TE_list_full,"Sm"));
Eu_pos=find(contains(TE_list_full,"Eu"));
Gd_pos=find(contains(TE_list_full,"Gd"));
data_Eu_anomaly=data_N(:,Eu_pos)./(sqrt(data_N(:,Sm_pos).*data_N(:,Gd_pos)));   

%% REE patterns (Fig. 1)

figure(1); clf(1)
ymin = 0.1; ymax = 1000;

%data
    subplot(1,2,1);
    
    for i = 1:length(data_N(:,1))
        semilogy(0:nTE,data_N(i,:),'k');
        hold on
    end
    
    xticks(0:nTE); xticklabels(TE_list_full)
    ylabel('Normalized REE concentrations')  
    ylim([ymin ymax])

%model
    subplot(1,2,2); 
    zoning_index=[1,5,10,15,20,25,30,35,40,45,51]; %positions along core-to-rim profile
    cmap=parula(length(zoning_index));
    
    % semilogy(0:nTE,solid_N(:),'k'); hold on %plot input liquid
    % semilogy(0:nTE,liquid_N(:),'--r'); %plot input cpx
    % semilogy(0:nTE,eq_solid_N(:),'r'); %plot equilibrated cpx
    
    % bottom cpx at time_index_max (e.g. 20kyr)
    depth_index = 1; %out of 100 (from 1-bottom to 100-top)
    time_index=time_index_max;
    baseFileName = sprintf('savetime_%d', time_index);
    load(fullfile(folderpath,baseFileName));       
    
    for i = 1:length(zoning_index)
        TE_pattern_aux(zoning_index(i),:) = ME_save_tp{iTP}(depth_index*nTE-nTE+1:depth_index*nTE,zoning_index(i))./N'; %normalization
        Pm(zoning_index(i)) = sqrt(TE_pattern_aux(zoning_index(i),Pm_pos-1)*TE_pattern_aux(zoning_index(i),Pm_pos)); %Pm interpolation
        TE_pattern(zoning_index(i),:) = [TE_pattern_aux(zoning_index(i),1:Pm_pos-1) Pm(zoning_index(i)) TE_pattern_aux(zoning_index(i),Pm_pos:end)];
        semilogy(0:nTE,TE_pattern(zoning_index(i),:),"LineStyle","-","Color",cmap(i,:));
        hold on
    end
    
    % 30m cpx at 25% time_index_max (e.g. 5kyr)
    depth_index = 3; %out of 100 (from 1-bottom to 100-top)
    time_index=round(0.25*time_index_max);
    baseFileName = sprintf('savetime_%d', time_index);
    load(fullfile(folderpath,baseFileName));    
    
    for i = 1:length(zoning_index)
        TE_pattern_aux(zoning_index(i),:) = ME_save_tp{iTP}(depth_index*nTE-nTE+1:depth_index*nTE,zoning_index(i))./N'; %normalization
        Pm(zoning_index(i)) = sqrt(TE_pattern_aux(zoning_index(i),Pm_pos-1)*TE_pattern_aux(zoning_index(i),Pm_pos)); %Pm interpolation
        TE_pattern(zoning_index(i),:) = [TE_pattern_aux(zoning_index(i),1:Pm_pos-1) Pm(zoning_index(i)) TE_pattern_aux(zoning_index(i),Pm_pos:end)];
        model1=semilogy(0:nTE,TE_pattern(zoning_index(i),:),"LineStyle","--","Color",cmap(i,:));
        hold on
    end
    
    xticks(0:nTE); xticklabels(TE_list_full)
    ylabel('Normalized REE concentrations')  
    ylim([ymin ymax])

%% Covariation (Eu/Eu*)N vs Eu (Fig. 2)

figure(2); clf(2)
xmin = 0.1; xmax = 10;
ymin = 0.4; ymax = 2;

%model
    depth_index=1;
    time_interval=10; %plot every 1/n of time_index_max

    i=0;
    for time_index=1:time_index_max
        if rem(time_index,round(time_index_max/time_interval,0))==0
            baseFileName = sprintf('savetime_%d', time_index);
            load(fullfile(folderpath,baseFileName));    
            i=i+1;
            Eu(i,:) = ME_save_tp{iTP}(depth_index*nTE-nTE+Eu_pos-1,:);
            Sm_N(i,:) = ME_save_tp{iTP}(depth_index*nTE-nTE+Sm_pos-1,:)./N(Sm_pos-1);
            Eu_N(i,:) = ME_save_tp{iTP}(depth_index*nTE-nTE+Eu_pos-1,:)./N(Eu_pos-1);
            Gd_N(i,:) = ME_save_tp{iTP}(depth_index*nTE-nTE+Gd_pos-1,:)./N(Gd_pos-1);
            Eu_anomaly_N(i,:) = Eu_N(i,:)./(sqrt(Sm_N(i,:).*Gd_N(i,:)));
            plot(Eu(i,:), Eu_anomaly_N(i,:));
            hold on
        end
    end
    set(gca,'colororder',parula(i))

%data   
    scatter(data(:,Eu_pos-1),data_Eu_anomaly(:),10,'filled','k');                                               
    hold on

xlabel('Eu (ppm)')
ylabel('(Eu/Eu*)_N')
xlim([xmin xmax]);
ylim([ymin ymax]);  
set(gca,'xscale','log')
