%Version 1.0
%Licenced by GPLv3
%Free to use share and adapt
%Appropriate credits given to Leo Svenningsson and relevant cited article

clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input from user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Values= importRaman('RamanData.txt',1, 15);  %Experimental data text file
conf330 = Values(1);
conf3390 = Values(2);
conf310 = Values(3);
conf3190 = Values(4);
conf3145 = Values(5);
conf210RAS = Values(6);
conf230RAS = Values(7);
IF = Values(8);%IF is the instrumental factor as in Yang and Michielsen
STDconf330 = Values(9);
STDconf3390 = Values(10);
STDconf310 = Values(11);
STDconf3190 = Values(12);
STDconf3145 = Values(13);
STDconf210RAS = Values(14);
STDconf230RAS = Values(15);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=10; % stepsize to calculate derivative for error estimation
Idh=[0,conf330,conf3390,conf310,conf3190,conf3145,conf210RAS,conf230RAS]+h;
STDI=[STDconf330,STDconf3390,STDconf310,STDconf3190,STDconf3145,STDconf210RAS,STDconf230RAS];
IN={'void','conf330','conf3390','conf310','conf3190','conf3145','conf210RAS','conf230RAS'};
syms a1 a2 b p2 p4
P2dh=zeros(1,length(IN)-1);
P4dh=zeros(1,length(IN)-1);
dP2dh=zeros(1,length(IN)-1);
dP4dh=zeros(1,length(IN)-1);
P2=0;
P4=0;
for i=1:length(IN)
    if i==1
        
    else
        assignin('base', IN{i}, Idh(i));
    end
    
    

    I330=conf330/conf330;
    I110=conf310*conf3390/conf3190/conf330;
    I310=IF*conf310/conf330;
    I3145=IF*conf3145/conf330;
    I21RAS=IF*conf310*conf210RAS/conf230RAS/conf330;
    
    alpha22 = b*((3*a1^2 + 3*a2^2 + 3 + 2*a1*a2 + 2*a1 + 2*a2)/15 + p2*(3*a1^2 + 3*a2^2 - 6 + 2*a1*a2 - a1 - a2)/21 + 3*p4*(3*a1^2 + 3*a2^2 + 8 + 2*a1*a2 - 8*a1 - 8*a2)/280) == I110;
    
    alpha33 = b*((3*a1^2 + 3*a2^2 + 3 + 2*a1*a2 + 2*a1 + 2*a2)/15 - 2*p2*(3*a1^2 + 3*a2^2 - 6 + 2*a1*a2 - a1 - a2)/21 + p4*(3*a1^2 + 3*a2^2 + 8 + 2*a1*a2 - 8*a1 - 8*a2)/35) == I330;
    
    alpha21 = b*((a1^2 + a2^2 + 1 + - a1*a2 - a1 - a2)/15 + p2*(a1^2 + a2^2 - 2 - 4*a1*a2 + 2*a1 + 2*a2)/21 + p4*(3*a1^2 + 3*a2^2 + 8 + 2*a1*a2 - 8*a1 - 8*a2)/280) == I21RAS;
    
    alpha23 = b*((a1^2 + a2^2 + 1 - a1*a2 - a1 - a2)/15 - p2*(a1^2 + a2^2 - 2 - 4*a1*a2 + 2*a1 + 2*a2)/42 - p4*(3*a1^2 + 3*a2^2 + 8 + 2*a1*a2 - 8*a1 - 8*a2)/70) == I310;
    
    alphaI31deg45 = b*((a1^2 + a2^2 + 1 - a1*a2 - a1 - a2)/15 - 1/2*p2*(a1^2 + a2^2 - 2 - 4*a1*a2 + 2*a1 + 2*a2)/21 + 19/32*p4*(3*a1^2 + 3*a2^2 + 8 + 2*a1*a2 - 8*a1 - 8*a2)/35) == I3145;
    
    %alpha2233 = b*((a1^2 + a2^2 + 1 + 4*a1*a2 + 4*a1 + 4*a2)/15 - p2*(a1^2 +
    %a2^2 - 2 + 10*a1*a2 - 5*a1 - 5*a2)/42 - p4*(3*a1^2 + 3*a2^2 + 8 + 2*a1*a2
    %- 8*a1 - 8*a2)/70) == -0.100; not in use, used in Citra et al.
    
    [sola1, sola2, solb, solp2, solp4,] = solve([alpha22, alpha33, alpha21, alpha23, alphaI31deg45],[a1, a2, b, p2, p4]);
    
    P2all=double(solp2);
    P4all=double(solp4);
    B=double(solb);
    a1all=double(sola1);
    a2all=double(sola2);
    
    if i==1
        
        if (35*P2all(1)^2 - 10*P2all(1) - 7)/18 <= P4all(1) && P4all(1) <= (5*P2all(1) + 7)/12
            P2=P2all(1);
            P4=P4all(1);
        end
        if (35*P2all(3)^2 - 10*P2all(3) - 7)/18 <= P4all(3) && P4all(3) <= (5*P2all(3) + 7)/12 && P2==0
            P2=P2all(3);
            P4=P4all(3);
        end
        if (35*P2all(3)^2 - 10*P2all(3) - 7)/18 <= P4all(3) && P4all(3) <= (5*P2all(3) + 7)/12 && P2~=0
            disp('more then one physical solution')
            return
        end
        if P2==0 && P4==0
            disp('no physical solution')
            return
        end
        
    end
    
    if i>1
        
        if (35*P2all(1)^2 - 10*P2all(1) - 7)/18 <= P4all(1) && P4all(1) <= (5*P2all(1) + 7)/12
            P2dh(i-1)=P2all(1);
            P4dh(i-1)=P4all(1);
        end
        if (35*P2all(3)^2 - 10*P2all(3) - 7)/18 <= P4all(3) && P4all(3) <= (5*P2all(3) + 7)/12 && P2==0
            P2dh(i-1)=P2all(3);
            P4dh(i-1)=P4all(3);
        end
        if (35*P2all(3)^2 - 10*P2all(3) - 7)/18 <= P4all(3) && P4all(3) <= (5*P2all(3) + 7)/12 && P2~=0
            disp('more then one physical solution')
            return
        end
        if P2==0 && P4==0
            disp('no physical solution')
            return
        end
        dP2dh(i-1)=(P2dh(i-1)-P2)/h;
        dP4dh(i-1)=(P4dh(i-1)-P4)/h;
    end
    
end
STDP2=sqrt(sum((dP2dh.*STDI).^2));
STDP4=sqrt(sum((dP4dh.*STDI).^2));

%%%%%%%%%%%%%%%%%%%
%You can force specific specific P2 and P4 by uncommenting
%%%%%%%%%%%%%%%%%%%
%  P2=0.46;
%  P4=0.28;
%  if (35*P2^2 - 10*P2 - 7)/18 <= P4 && P4 <= (5*P2 + 7)/12
%  elseSTDconf3390
%      disp('no physical solution')
%      return
%  end
%%%%%%%%%%%%%%%%%%%


res=401;
thetaV=linspace(0,pi,res);

%%%%%%%%%%%%%%%%%%%
%calculate Wrapped Lorentzian distribution
%%%%%%%%%%%%%%%%%%%
f = @(gamma) legendreWL(gamma,P2,P4);

initWL= 1;

[gamma,fvalWL] = fminsearch(f,initWL);
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%plot wrapped lorentzian distribution
%%%%%%%%%%%%%%%%%%%
WL=1/(pi)*sinh(gamma)./(cosh(gamma)-cos(2*thetaV))/trapz(thetaV,1/(pi)*sinh(gamma)./(cosh(gamma)-cos(2*thetaV)));
IntegralWL=trapz(thetaV,WL);
figure(1)
plot(thetaV(1:(res-1)/2 +1)*180/pi,WL(1:(res-1)/2 +1),'LineWidth',4);
xlabel('Angle \theta')
ylabel('ODF f(\theta)')
title('Wrapped Lorentzian ODF')
%%%%%%%%%%%%%%%%%%%

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
G=sqrt(mandphi(1)/pi)*exp(-mandphi(1)*(thetaV-mandphi(2)).^2);
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

%%%%%%%%%%%%%%%%%%%
%plot legendre expansion
%%%%%%%%%%%%%%%%%%%
%k1 =1/(2*pi)*(4+1)/2;
%k2 =1/(2*pi)*(8+1)/2;
%Lexp = k1*P2*(3*cos(theta).^2 - 1)/2 + k2*P4*(35*cos(theta).^4 - 30*cos(theta).^2 + 3)/8;
%figure(5)
%plot(theta*180/pi,Lexp)
%%%%%%%%%%%%%%%%%%%
disp(['P2 is calculated to ',num2str(round(P2,4,'significant')),' with STD of ',num2str(round(STDP2,1,'significant'))])
disp(['P4 is calculated to ',num2str(round(P4,4,'significant')),' with STD of ',num2str(round(STDP4,1,'significant'))])
