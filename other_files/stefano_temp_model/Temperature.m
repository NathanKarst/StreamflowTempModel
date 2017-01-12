
%The streamflow divided channel by channel is Qout_net
load Qout_net %Matrix of water discharge exiting each channel at each time step
load Qin_net %Matrix of water discharge entering each channel at each time step
load V % Matrix of water volumes in each channel at each time step
load d %Matrix of water depths in each channel at each time step
load Hill_discharge %Matrix of water discharge entering the channel from the hillslope at each time step
load elder30m.asc % Digital elevation model raster file
load channels150.asc %Channel network raster file
load watersheds150.asc %Watersheds raster file
load Ta % Atmospheric temperature at each time step
load EV %Evaporation from open water at each time step
load daydlhour %Number of daylight hours at the daily scale
load awv %Actual water vapor in Pascal
load watershed_area %Area of each subwatershed of the basin [m^2] 

Qout_net=Qout_net';
Qin_net=Qin_net';
V=V';
d=d';
Hill_discharge=Hill_discharge';

if exist('day_start')
else
     day_start=datenum(2010,8,1);
     day_end=datenum(2011,8,1);
end

%Light data are weekly so I need to convert the discharge into a weekly
%time scale.
if exist('dtk')
    dt=dtk;
else
    dt=6; %Original Time step [min]
end
dt_T=8; % time step for the temperature model [min]
N=length(Qout_net(1,:));
N_T=ceil(length(Qout_net(1,:))/(dt_T/dt));

Hill_dischargeT=zeros(length(link_no),N_T);
 Qout_netT=zeros(length(link_no),N_T);
 Qin_netT=zeros(length(link_no),N_T);
 V_T=zeros(length(link_no),N_T);
 d_T=zeros(length(link_no),N_T);
 Ta_T=zeros(1,N_T);
 EV_T=zeros(1,N_T);
 awv_T=zeros(1,N_T);
for i=1:length(link_no)
   Qout=reshape(Qout_net(i,:),dt_T/dt,N_T);
   Qin=reshape(Qin_net(i,:),dt_T/dt,N_T);
   hillQ=reshape(Hill_discharge(i,:),dt_T/dt,N_T);
   V1=reshape(V(i,:),dt_T/dt,N_T);
   d1=reshape(d(i,:),dt_T/dt,N_T);  
   for j=1:N_T
   %mean_new_Q is the discharge produced or received in the new time period by each channel, considering everything
   %that happens upstream.
   Qout_netT(i,j)=sum(Qout(:,j));
   Qin_netT(i,j)=sum(Qin(:,j));
   %mean_new_hill_Q is the discharge produced in the new time period by each sub-basin hillslope only.
   Hill_dischargeT(i,j)=sum(hillQ(:,j));
   %mean_new_V is the volume in the new time period for each channel.
   V_T(i,j)=mean(V1(:,j));
    %mean_new_d is the depth in the new time period for each channel.
   d_T(i,j)=mean(d1(:,j));
   end
end
 Ta1=reshape(Ta,dt_T/dt,N_T);
 EV1=reshape(EV,dt_T/dt,N_T);
 awv1=reshape(awv,dt_T/dt,N_T);
for j=1:N_T
    Ta_T(j)=mean(Ta1(:,j));
    EV_T(j)=sum(EV1(:,j)); % mm in the time step dt (note that dt is note the same as dt in "model_network_new").
    awv_T(j)=mean(awv1(:,j));
end

% This is to extract the light in x channel pixel. I need to do it only
% if I need to built a new matrix with other light data. Otherwise it
% should be saved as M.mat in the folder. M.mat has a number of rows equal
% to the number of channel pixels and a number of columns equal to 365. The
% index of the rows can be univoccally traced back to the position in the
% channel using chanlightind1 and chanlightind2.

chanlightind1=[];
chanlightind2=[];
for in=1:length(link_no)
    tmp=find(channels150==1 & watersheds150==link_no(in));
    for i=1:length(find(channels150==1 & watersheds150==link_no(in)))
    chanlightind1=[chanlightind1 link_no(in)];
    chanlightind2=[chanlightind2 tmp(i)];
    end
end
tmp=[];

if 0
j=1;
lweek=51; % number of weeks for he light dataset
chanpix=length(find(channels150==1));
M=zeros(chanpix,lweek);

for i=5:7:362
    if i<10
f= strcat('C:\Users\stefano\Desktop\ssr\elder30myr4d00',num2str(i),'.asc');
    else
        if i<100
            f= strcat('C:\Users\stefano\Desktop\ssr\elder30myr4d0',num2str(i),'.asc');
        else
            f= strcat('C:\Users\stefano\Desktop\ssr\elder30myr4d',num2str(i),'.asc');
        end
    end
MM=dlmread(f, '',6,0);
for in=1:length(link_no)
M(chanlightind1==link_no(in),j)=MM(channels150==1 & watersheds150==link_no(in));
end
j=j+1;
end
%questo spalma la matrice settimanale della luce in scala giornaliera. La
%matrice che viene fuori ha un numero di colonnne pari a 364 (i.e., 62*7) e un numero di
%righe pari al numero di channel pixels.
M=reshape(repmat(M,7,1),chanpix,7*length(M(1,:)));
%Devo aggiungere un'altra colonna per arrivare a 365, la aggiungo alla fine
%e la metto uguale al primo giorno:
M=[M M(:,1)];
if 0
if length(M(1,:))<length(mean_new_hill_Q(1,:))
    for i=length(M(1,:)):length(mean_new_hill_Q(1,:))-1
    M=horzcat(M,M(:,end));
    end
end
end
save('M.mat','M');
else
    load M.mat;
end

% Light matrix goes from January 1st to December 31st. The following is to
% make it start from the first day of the hydrologic year based on date_start.
j=0;
while j<ceil((day_end-day_start)/365)
    M=[M M];
    j=j+1;
end
j=[];
tmp33=datevec(day_start);
tmp331=day_start-datenum(tmp33(1),1,1);
newM=[M(:,tmp331+1:365) M];
%newM(:,day_end-day_start+1:end)=[];
newM(:,day_end-day_start+1:end)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Units of radiation (newM) are Watt-hours/(m^2*day) -> [M*T^-3]
% transform it in Watt-hours/(m^2*min) and change name
phi=newM/(60*24); %[Watt-hours/(m^2*min)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This spreads the matrix phi over the time of simulation using the step of
%simulation dt_T.
 phi=reshape(repmat(phi,24*60/dt_T,1),length(phi(:,1)),length(phi(1,:))*24*60/dt_T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build a vector of zeros and ones on the dt_T time scale to consider phi
% only for the hours of daylight and then multiply it by phi to obtain
% solar radiation only during the day
tmp=zeros(length(phi(1,:)),1);
tmp1=ones(length(phi(1,:)),1);
daydlhour*60/dt_T;
tmp2=round(daydlhour*60/dt_T/2);
tmp3=0;
for j=1:day_end-day_start
    for i=tmp3+12*60/dt_T-tmp2(j):tmp3+12*60/dt_T+tmp2(j)
    tmp(i)=1;
    end
    tmp3=tmp3+24*60/dt_T;
end

for i=1:size(phi,1)
    phi(i,:)=phi(i,:).*tmp';
end

