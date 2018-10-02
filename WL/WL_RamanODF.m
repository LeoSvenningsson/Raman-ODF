%Version 1.0
%Licenced by GPLv3
%Free to use share and adapt
%Appropriate credits given to Leo Svenningsson and relevant cited article

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input from user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Values= importRaman('WL_RamanData.txt',1, 15);  %Experimental data text file

init=[0.4,0.5,-0.5,1];                             %initialcondition for the for the iterative solver [gamma,alpha1,alpha2,alpha3]
h=0.000001;                                              %stepsize to calculate derivative for error estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conf330 = Values(1)/Values(1);
conf3390 = Values(2)/Values(1);
conf310 = Values(3)/Values(1);
conf3190 = Values(4)/Values(1);
conf3145 = Values(5)/Values(1);
IF = Values(8);
res=401;

thetaV=linspace(0,pi,res);

P2WL=zeros(1,6);
P4WL=zeros(1,6);

for i=1:5 % these are used to calculate derivatives for error estimation
    I=[conf330,conf3390,conf310,conf3190,conf3145];
    I(i)=I(i)+h; %adds delta I for derivatives calaculation
    I3300=I(1);
    I3390=I(3)*I(2)/I(4);
    I3100=IF*I(3);
    I3145=IF*I(5);
    
    f = @(vars) WLfun(vars,I3300,I3390,I3100,I3145);
    options = optimset('PlotFcns',@optimplotfval,'MaxFunEvals',10000,'TolFun',1e-14,'TolX',1e-14);%1e-6 default
    [ga1a2a3,fval2] = fminunc(f,init,options);
    
    
    WL=1/(pi*2)*sinh(ga1a2a3(1))./(cosh(ga1a2a3(1))-cos(2*thetaV));
    funWL=sin(thetaV).*WL;
    funWL1=sin(thetaV).*WL.*(3*cos(thetaV).^2 - 1)/2;
    P2WL(i+1)=trapz(thetaV,funWL1)/trapz(thetaV,funWL);
    funWL2=sin(thetaV).*WL.*(35*cos(thetaV).^4 - 30*cos(thetaV).^2 + 3)/8;
    P4WL(i+1)=trapz(thetaV,funWL2)/trapz(thetaV,funWL);
    
    
end

I3300=conf330;
I3390=conf310*conf3390/conf3190;
I3100=IF*conf310;
I3145=IF*conf3145;

f = @(vars) WLfun(vars,I3300,I3390,I3100,I3145);
options = optimset('PlotFcns',@optimplotfval,'MaxFunEvals',10000,'TolFun',1e-14,'TolX',1e-14);%1e-6 default
[ga1a2a3,fval2] = fminunc(f,init,options);
gamma=ga1a2a3(1);
WL=1/(pi)*sinh(gamma)./(cosh(gamma)-cos(2*thetaV));

funWL=sin(thetaV).*WL;
funWL1=sin(thetaV).*WL.*(3*cos(thetaV).^2 - 1)/2;
P2WL(1)=trapz(thetaV,funWL1)/trapz(thetaV,funWL);
funWL2=sin(thetaV).*WL.*(35*cos(thetaV).^4 - 30*cos(thetaV).^2 + 3)/8;
P4WL(1)=trapz(thetaV,funWL2)/trapz(thetaV,funWL);


dP2dh=zeros(1,5);
dP4dh=zeros(1,5);
for i=1:5
dP2dh(i)=(P2WL(i+1)-P2WL(1))/h;
dP4dh(i)=(P4WL(i+1)-P4WL(1))/h;
end
STDI=[Values(9),Values(10),Values(11),Values(12),Values(13)]/Values(1);

STDP2=sqrt(sum((dP2dh.*STDI).^2)); %final error estimation
STDP4=sqrt(sum((dP4dh.*STDI).^2)); %final error estimation

IntegralWL=trapz(thetaV,WL); 

P2=P2WL(1);
P4=P4WL(1); 

%%%%%%%%%%%%%%%%%%%
%cWrapped Lorentzian distribution
%%%%%%%%%%%%%%%%%%%

figure(1)
plot(thetaV(1:(res-1)/2 +1)*180/pi,WL(1:(res-1)/2 +1),'LineWidth',4);
xlabel('Angle \theta')
ylabel('ODF f(\theta)')
title('Wrapped Lorentzian ODF')


%%%%%%%%%%%%%%%%%%%
%calculate most probable distribution
%%%%%%%%%%%%%%%%%%%
f = @(lambda) legendreMP(lambda,P2,P4);

initMP= [1.1,4];

[lambda1and2,fvalMP] = fminsearch(f,initMP);
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%plot most probable distribution
%%%%%%%%%%%%%%%%%%%
MP=exp(lambda1and2(1)*(3*cos(thetaV).^2 - 1)/2 + lambda1and2(2)*(35*cos(thetaV).^4 - 30*cos(thetaV).^2 + 3)/8);
MP=MP/(integral(@(theta) exp(lambda1and2(1)*(3*cos(theta).^2 - 1)/2 + lambda1and2(2)*(35*cos(theta).^4 - 30*cos(theta).^2 + 3)/8),0,pi));
IntegralMP=trapz(thetaV,MP);
figure(2)
plot(thetaV(1:(res-1)/2 +1)*180/pi,MP(1:(res-1)/2 +1),'LineWidth',4)
xlabel('Angle \theta')
ylabel('ODF f(\theta)')
title('Most Probable ODF')
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%calculate gauss distribution
%%%%%%%%%%%%%%%%%%%
f = @(mphi) legendreGauss(mphi,P2,P4);

initGauss= [0.5,0];
con1a=[];
con1b=[];
con2a=[];
con2b=[];
lb=[0,0];
ub=[10^3,pi/2];
nonlcon=[];
options = optimoptions('fmincon','display','none');
[mandphi,fvalGauss] = fmincon(f,initGauss,con1a,con1b,con2a,con2b,lb,ub,nonlcon,options);
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%plot gauss distribution
%%%%%%%%%%%%%%%%%%%
G=exp(-mandphi(1)*(thetaV-mandphi(2)).^2);
G=G/(integral(@(theta) exp(-mandphi(1)*(theta-mandphi(2)).^2),-10^3,10^3));
IntegralGauss=trapz(thetaV,G);
figure(3)
plot(thetaV(1:(res-1)/2 +1)*180/pi,G(1:(res-1)/2 +1),'LineWidth',4,'color',[0.8500 0.3250 0.0980])
xlabel('Angle \theta')
ylabel('ODF f(\theta)')
title('Gaussian ODF')
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%plot all
%%%%%%%%%%%%%%%%%%%
figure(4)
plot(thetaV(1:(res-1)/2 +1)*180/pi,WL(1:(res-1)/2 +1),thetaV(1:(res-1)/2 +1)*180/pi,MP(1:(res-1)/2 +1),thetaV(1:(res-1)/2 +1)*180/pi,G(1:(res-1)/2 +1),'LineWidth',4)
xlabel('Angle \theta')
ylabel('ODF f(\theta)')
legend('Wrapped Lorentzian','Most Probable','Gauss')
title('All ODFs')
%%%%%%%%%%%%%%%%%%%



disp(['P2 is calculated to ',num2str(round(P2,4,'significant')),' with STD of ',num2str(round(STDP2,1,'significant'))])
disp(['P4 is calculated to ',num2str(round(P4,4,'significant')),' with STD of ',num2str(round(STDP4,1,'significant'))])
