%Energy Balance equation:

%dT/dt=(QinTin+QgwTgw-QoutTout)/V+(1-alphaw)phi/(d*rho*cp) - Lout + Lin
%      - H - lambda*rho*EV

%alphaw = albedo of water, assumed 0.15 in Van Beek et al. 2012

%Lout = outgoing longwave radiation
%Lout=eps*sigma*T^4 where T is the stream temperature, eps = emissivity
%which is chosen equal to 1 from Van Beek et al. 2012 and sigma 
%is the boltzmann constant
eps=1;
sigma=4.903*10^(-9); %[MJ/K^4/m^2/day]
%trasform sigma in Wh/K^4/m^2/min
sigma=sigma*0.000278*1000*1000/1440; %[Wh/K^4/m^2/min]
%Lin= incoming longwave radiation, is calculated here using Konzelmann et al. (1994)
% or Equations 1 and 4 in Sedlar Hock 2008 (Table 2): Lin=epscs*F*sigma*Ta_T^4. 
%epscs is the clear sky emissivity and depends on atmospheric pressure and
%water vapor: awv_T. F is the cloud cover, we will assume it equal to 1 for now
%Compute epscs:
for i=1:size(phi,2)
    epscs(i)=0.23+0.4393*(awv_T(i)/(273.16+Ta_T(i)))^(1/7);
end
F=1;
%Compute Lin
for i=1:size(phi,2) %[Wh/m^2/min] 
   Lin(i)=F*epscs(i)*sigma*(273.16+Ta_T(i))^4; %[Wh/m^2/min] 
end
%H = sensible heat flux: H=Kh*(T-Ta), T is the stream temperature. (Van Beek et al. 2012)
Kh= 20*0.000278*60; % [Wh/min/m^2/K] Turbolent heat exchange coefficient
%lambda*rho*EV is the latent heat flux, with lambda being the latent heat
%of vaporization:
lambda=2260*0.000278*1000; %[Wh/Kg] 



%cp= 4.1855 [kJ/(Kg·K)]specific heat capacity of water = 4.1855*0.000278 [kWh/(Kg*K)]
cp=0.001163569*1000; %[Wh/(Kg*K)]
rho=1000; %[Kg/m^3] water density
Tgw= 284; % [K] soil water temperature


%Intitial condition on temperature
T=zeros(size(phi));%Temperature
T(:,1)=Tgw;

flux_inT=zeros(size(Qout_netT));
Tmean=zeros(size(Qout_netT));
TmeanRK=zeros(size(Qout_netT));
%TmeanRK=zeros(length(Qout_netT(:,1)),round(length(Qout_netT(1,:))/2));
Tmean(:,:)=284;
TmeanRK(:,:)=284;

        f1 = @(tem,fluxin1,Qout1,Hill1,V1) (fluxin1-Qout1*tem+Tgw*Hill1)/(dt_T*V1);
        f2 = @(tem,d1)                      eps*sigma*(tem^4)/(d1*rho*cp);
        f3 = @(L1,d1)                       L1/(d1*rho*cp);
        f4 = @(phi1,d1)                     mean(phi1)/(d1*rho*cp);
        f5 = @(tem,Ta1,d1)                  Kh*(tem-(273.16+Ta1))/(d1*rho*cp);
        f6 = @(EV1,d1)                      lambda*rho*EV1*0.001/(dt_T*d1*rho*cp);

%         fT = @(tem,fluxin1,Qout1,Hill1,V1,d1,phi1,Ta1,EV1) ...
%             f1(tem,fluxin1,Qout1,Hill1,V1)...
%             -f2(tem,d1)...
%             +f3(L1,d1)...
%             +f4(phi1,d1)...
%             -f5(tem,Ta1,d1)...
%             -f6(EV1,d1);
          
fT = @(tem,fluxin1,Qout1,Hill1,V1,d1,L1,phi1,Ta1,EV1) ...
    f1(tem,fluxin1,Qout1,Hill1,V1)...
                       -f2(tem,d1)...
            +f3(L1,d1)...
            +f4(phi1,d1)...
            -f5(tem,Ta1,d1)...
            -f6(EV1,d1);

for in=1:length(link_no)
        in  
        
        
        T_=Tmean(in,:);
        TRK_=TmeanRK(in,:);
        fluxin_=flux_inT(in,:);
        Qout_=Qout_netT(in,:);
        Hill_=Hill_dischargeT(in,:);
        d_=d_T(in,:);
        L_=Lin(:);
        
        Ta_=Ta_T(:);
        EV_=EV_T(:);
        V_=V_T(in,:);
        
    for t=7000:length(phi(1,:))-2
        if min(d_T(:,t))>0
        
            mask=find(chanlightind1==link_no(in));
        
        T0=T_(t);
        TRK0=TRK_(t);
        fluxin0=fluxin_(t);
        fluxin1=fluxin_(t+1);
        fluxin2=fluxin_(t+2);
        Qout0=Qout_(t);
        Qout1=Qout_(t+1);
        Qout2=Qout_(t+2);
        Hill0=Hill_(t);
        Hill1=Hill_(t+1);
        Hill2=Hill_(t+2);
        d0=d_(t);
        d1=d_(t+1);
        d2=d_(t+2);
        L0=L_(t);
        L1=L_(t+1);
        L2=L_(t+2);
        phi0=phi(mask,t);
        phi1=phi(mask,t+1);
        phi2=phi(mask,t+2);
        Ta0=Ta_(t);
        Ta1=Ta_(t+1);
        Ta2=Ta_(t+2);
        EV0=EV_(t);
        EV1=EV_(t+1);
        EV2=EV_(t+2);
        V0=V_(:,t);
        V1=V_(:,t+1);
        V2=V_(:,t+2);
        
%EULER METHOD
        Tmean(in,t+1)=T0+dt_T*fT(T0,fluxin0,Qout0,Hill0,V0,d0,L0,phi0,Ta0,EV0);
% RUNGE KUTTA METHOD        
        %Check that fluxin comes from TmeanRK and not Tmean in this case
        k1=2*dt_T*fT(TRK0,fluxin0,Qout0,Hill0,V0,d0,L0,phi0,Ta0,EV0);
        k2=2*dt_T*fT(TRK0+0.5*k1,fluxin1,Qout1,Hill1,V1,d1,L1,phi1,Ta1,EV1);
        k3=2*dt_T*fT(TRK0+0.5*k2,fluxin1,Qout1,Hill1,V1,d1,L1,phi1,Ta1,EV1);
        k4=2*dt_T*fT(TRK0+k3,fluxin2,Qout2,Hill2,V2,d2,L2,phi2,Ta2,EV2);
        
        TmeanRK(in,t+1)=TRK0+(1/6)*(k1+2*k2+2*k3+k4);

        
        








        
        mask=[];
       
       
        else
            Tmean(in,t)=284;
        end 
         
    end
    % Update the following fluxes
    if in<length(link_no)
    mask=find(link_no==downstream_link_no(in));
    
    tmp3=[0 Qout_netT(in,:).*Tmean(in,:)];
    tmp3(end)=[];
    %keyboard
    flux_inT(mask,:)=flux_inT(mask,:)+tmp3;
    
    end
end

 xx=Tmean(1,:);
 yy=d_T(1,:);
  x=reshape(xx,720/dt_T,(day_end-day_start)*2); %-day_start+1
  y=reshape(yy,720/dt_T,(day_end-day_start)*2); %-day_start+1
  for i=1:length(x(1,:))
      xxx(i)=mean(x(:,i));
      yyy(i)=mean(y(:,i));
  end
  if 0
      prova=[1 10 14 16 26 46 54 56];
      plot(Tmean(prova(1),:))
hold on
plot(Tmean(prova(3),:),'r')
plot(Tmean(prova(6),:),'g')
plot(Tmean(prova(7),:),'k')
legend('1st order stream','2nd order stream','3rd order stream','4th order stream','Location','NorthWest')
  end