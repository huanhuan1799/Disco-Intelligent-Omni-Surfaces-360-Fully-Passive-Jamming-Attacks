function [Hd,G,Hr,path_g,path_r,path_d,path_AJ,noise] = J_generate_DIOSchannel_nearfield(WW,HH,N,K,eb1,eb2)
%%    
%clc;close all; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NNN = WW*HH;
% for n = 1:1000
AP = [0,0,10];%%基站的位置
%Dar = 120;
xxx = 2;
yyy = 2;
IRS = [xxx,yyy,8]; %智能反射面的位置

AJ = IRS; 

%% 用户位置
Rai = 20;
aaa = rand(1,K);
clear Tha
Tha(1:K/2) = pi*( aaa(1:K/2)  ) - pi/2;
Tha(K/2+1:K) = pi/2+pi*( aaa( (K/2 + 1):K) );
 
rad(1:K) = Rai*sqrt(rand(1,K));
Dar1  = 0;
Dar2 = 180;

for k = 1:K
    Pt(k,:) = [Dar1+rad(k)*cos(Tha(k)),Dar2+rad(k)*sin(Tha(k)),0];    
end
% for k = 1:K/2
% plot(Pt(k,1),Pt(k,2),'b.'); hold on
% end
% for k = K/2+1:K
% plot(Pt(k,1),Pt(k,2),'r.'); hold on
% end
% axis([-20,20,160,200]);
% end


noise = -170+10*log10(180*1e3);%设置噪声功率
%noise = -120;

    for k = 1:K  %对用户循环
       %% channel Hr_w
        dAJ(k) = sqrt(sum(abs(Pt(k,:)-AJ).^2));%IRS到user k的距离
        LAJ(k) = 32.6 + 36.7*log10(dAJ(k));%c0=-30dB,a = 2.8;
        path_AJ(k) = 10.^((-LAJ(k))/10);
        pAJ(k) = sqrt(path_AJ(k)*10.^(-noise/10));%IRS-user的路径损耗
    end

%% AP-DIOS links
G_LOS = zeros(NNN,N);

    d_g = sqrt(sum(abs(AP-IRS).^2));%ap和IRS之间的距离
    Lg = 35.6 + 22*log10(d_g) ;%c0=-30dB,a = 2; 
    path_g = 10.^((-Lg)/10);
    pg = sqrt(path_g*10.^(-noise/10)); %AP-IRS的路径损耗
   %% channels G(r)
    G_NLOS = sqrt(1/2).*(randn(NNN,N)+1j.*randn(NNN,N));
    Ddiff = zeros(NNN,N); 
    for  w = 1:WW
        for h = 1:HH
            for n = 1:N
                r = (w-1)*HH+h;
                Ddiff(r,n) = norm([(n-1)*0.05/2 , 0 ,2] - [xxx+(w-1)*0.05/2,yyy+(h-1)*0.05/2, 2],'fro');
                G_LOS(r,n) = exp(-1j*(2*pi/0.05)*Ddiff(r,n));
            end
        end
    end
    G = eb1.*G_LOS+eb2.*G_NLOS;
    G = pg.*G;

    
%% DIOS-User links
Hr_sig = sqrt(1/2).*(randn(K,NNN)+1j.*randn(K,NNN));
    for k = 1:K  %对用户循环
       %% channel Hr_w
        dr(k) = sqrt(sum(abs(Pt(k,:)-IRS).^2));%IRS到user k的距离
        Lr(k) = 32.6 + 36.7*log10(dr(k));%c0=-30dB,a = 2.8;
        path_r(k) = 10.^((-Lr(k))/10);
        pr(k) = sqrt(path_r(k));%IRS-user的路径损耗
        Hr(k,:) = pr(k).*Hr_sig(k,:);
    end


for k = 1:K
    %% channel Hd_w    
    dk(k) = sqrt(sum(abs(Pt(k,:)-AP).^2)); %ap到user的距离
    Ld(k) = 32.6 + 36.7*log10(dk(k)) ;% 36.7
    path_d(k) = 10.^((-Ld(k))/10);
    pd(k) = sqrt(path_d(k)*10.^(-noise/10));
    Hd(k,:)=sqrt(1/2).*(randn(1,N)+1j.*randn(1,N));
    Hd(k,:) = pd(k).*Hd(k,:);
end

