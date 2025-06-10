close all;  clc;clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Simulation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle'); %Initiate the random number generators with a random seed
K = 24; %Number of users
N = 128; %Number of transmit antennas
Niter = 2000;%Number of realizations in the Monte Carlo simulations
WW = 64;
HH = 32;
NN = WW*HH; %Number of IRS elements
%Range of SNR values
Pmin = -20;
Pmax = 20;
PdBm = [Pmin:5:Pmax]; %dB scale
P = (10.^(PdBm/10)); %Linear scale
P = P*K;
AJ_dB = (10.^(5/10)); % dB
%%%%%%%%%%%%%%%%%%%%%%For Estimation
Tc = 6; %% the length of a DT phase is C(=6) times longer than that of a PT phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
User_Sig_WO = zeros(K,length(PdBm));
InterUser_Inter_WO = zeros(K,length(PdBm));
User_Sig_AJ = zeros(K,length(PdBm));
InterUser_Inter_AJ = zeros(K,length(PdBm));
User_Sig_GFPJ_CM = zeros(K,length(PdBm));
InterUser_Inter_GFPJ_CM = zeros(K,length(PdBm));
Theo_User_Sig_GFPJ_CM = zeros(K,length(PdBm));
Theo_InterUser_Inter_GFPJ_CM = zeros(K,length(PdBm));
User_Sig_GFPJ_OR = zeros(K,length(PdBm));
InterUser_Inter_GFPJ_OR = zeros(K,length(PdBm));
User_Sig_GFPJ_DM = zeros(K,length(PdBm));
InterUser_Inter_GFPJ_DM = zeros(K,length(PdBm));
Theo_User_Sig_GFPJ_DM = zeros(K,length(PdBm));
Theo_InterUser_Inter_GFPJ_DM = zeros(K,length(PdBm));
%% Generate Random Reflecting Vector
b = 1;
T_Amp = [0.78, 0.82]; % Refractive amplitudes
R_Amp = [0.62, 0.57]; % Reflective amplitudes
T_Omg = [5*pi/3,2*pi/3]; % Refractive phase shifts
R_Omg = [pi/9,7*pi/6]; % Reflective phase shifts
OR_Omg = [pi/9,7*pi/6]; % Reflective RIS model
for pv = 1:length(PdBm)
    pv
    for m = 1:Niter
        %% Channel;Only DIRS, not random wireless channels
        eb = 10; eb22 = 1/(1+eb); eb11 = eb/(1+eb); eb1 = sqrt(eb11); eb2 = sqrt(eb22); % Rician factor
        [Hd,Gr,Hr,lg,lr,ld,lAJ,noise] = J_generate_DIOSchannel_nearfield(WW,HH,N,K,eb1,eb2); % Channel generation
        Stacha_Acc_CM = lr*lg*NN* 0.5 *(10.^(-noise/10)); % mu = 0.5
        Stacha_Acc_DM(1:K/2) = lr(1:K/2)*lg*NN* 0.66 *(10.^(-noise/10)); % 1:K/2 Transmissive; K/2+1:K Reflective, mu = 0.66 (Due to T_Amp)
        Stacha_Acc_DM( (K/2+1):K ) = lr( (K/2+1):K )*lg*NN* 0.34 *(10.^(-noise/10)); % mu = 0.34 (Due to R_Amp)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% w/o Jamming and Existing Active Jammer%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        W = ( (Hd'*inv(Hd*Hd')) ) * diag(sqrt(P(pv)).*ones(1,K))./norm(Hd'*inv(Hd*Hd'),'fro'); % Zero-forcing precoding
        % norm(W,'fro')^2
        for kk = 1:K
            User_Sig_WO(kk,pv) =  User_Sig_WO(kk,pv) + norm( Hd(kk,:)*W(:,kk),'fro')^2;
            User_Sig_AJ(kk,pv) =  User_Sig_AJ(kk,pv) + norm( Hd(kk,:)*W(:,kk),'fro' )^2;
            User_set_Hd = [1:K];
            User_set_Hd(User_set_Hd == kk) = [];
            HD_age = 0;
            for u = 1:length(User_set_Hd)
                HD_age = HD_age +  norm(  Hd(kk,:)*W(:,User_set_Hd(u)),'fro' )^2;
                InterUser_Inter_AJ(kk,pv) = InterUser_Inter_AJ(kk,pv) +  norm( Hd(User_set_Hd(u),:)*W(:,kk),'fro' )^2;
            end
            InterUser_Inter_WO(kk,pv) =  InterUser_Inter_WO(kk,pv) + HD_age + 1;
            InterUser_Inter_AJ(kk,pv) = InterUser_Inter_AJ(kk,pv) +  lAJ(kk)*10.^(-noise/10)*10^(AJ_dB/10) +1;
        end
        %% Generation Pha during the RPT and DT phases
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T_betaRPT_DM = zeros(1,NN );
        R_betaRPT_DM = zeros(1,NN );
        T_phaDIRSRPT_DM = zeros(1,NN );
        R_phaDIRSRPT_DM = zeros(1,NN );
        T_betaRPT_CM = zeros(1,NN );
        R_betaRPT_CM = zeros(1,NN );
        T_phaDIRSRPT_CM = zeros(1,NN );
        R_phaDIRSRPT_CM = zeros(1,NN );
        OR_betaRPT_CM = zeros(1,NN );
        OR_phaDIRSRPT_CM = zeros(1,NN );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ri = 1:length(T_phaDIRSRPT_DM)
            iddele = randi(2^b,1,1)*randi(2^b,1,1) ;
            if iddele == 1
                T_phaDIRSRPT_DM(ri) = T_Omg(1); %%% Taking the first element of the phase shift set
                T_betaRPT_DM(ri) = T_Amp(1); %%% Taking the first amplitude element
                R_phaDIRSRPT_DM(ri) = R_Omg(1);
                R_betaRPT_DM(ri) = R_Amp(1);
                %%%% Constant Amplitude
                T_phaDIRSRPT_CM(ri) = T_Omg(1);
                T_betaRPT_CM(ri) = sqrt(2)/2;
                R_phaDIRSRPT_CM(ri) = R_Omg(1);
                R_betaRPT_CM(ri) = sqrt(2)/2;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reflective RIS
                OR_phaDIRSRPT_CM(ri) = OR_Omg(1);
                OR_betaRPT_CM(ri) = 1;
                
            else
                T_phaDIRSRPT_DM(ri) = T_Omg(2);
                T_betaRPT_DM(ri) = T_Amp(2);
                R_phaDIRSRPT_DM(ri) = R_Omg(2);
                R_betaRPT_DM(ri) = R_Amp(2);
                
                T_phaDIRSRPT_CM(ri) = T_Omg(2);
                T_betaRPT_CM(ri) = sqrt(2)/2;
                R_phaDIRSRPT_CM(ri) = R_Omg(2);
                R_betaRPT_CM(ri) = sqrt(2)/2;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reflective RIS
                OR_phaDIRSRPT_CM(ri) = OR_Omg(2);
                OR_betaRPT_CM(ri) = 1;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Precoding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T_theta_inbarRPT_DM = T_betaRPT_DM.*exp(1j.*(T_phaDIRSRPT_DM)); % The variable-amplitude IOS model
        R_theta_inbarRPT_DM = R_betaRPT_DM.*exp(1j.*(R_phaDIRSRPT_DM));
        T_HRPT_DM = Hd;  %% Channel in the PT phase, 
        W_DM = (T_HRPT_DM'*inv(T_HRPT_DM*T_HRPT_DM')) * diag(sqrt(P(pv)).*ones(1,K))./norm(T_HRPT_DM'*inv(T_HRPT_DM*T_HRPT_DM'),'fro');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T_theta_inbarRPT_CM = T_betaRPT_CM.*exp(1j.*(T_phaDIRSRPT_CM)); % The constant-amplitude IOS model
        R_theta_inbarRPT_CM = R_betaRPT_CM.*exp(1j.*(R_phaDIRSRPT_CM));
        T_HRPT_CM = Hd; %% Channel in the PT phase
        W_CM = (T_HRPT_CM'*inv(T_HRPT_CM*T_HRPT_CM'))* diag(sqrt(P(pv)).*ones(1,K))./norm(T_HRPT_CM'*inv(T_HRPT_CM*T_HRPT_CM'),'fro');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        OR_theta_inbarRPT_CM = OR_betaRPT_CM.*exp(1j.*(OR_phaDIRSRPT_CM)); % Reflective RIS
        OR_HRPT_CM = Hd; %% Channel in the PT phase
        W_CM_OR = (OR_HRPT_CM'*inv(OR_HRPT_CM*OR_HRPT_CM'))* diag(sqrt(P(pv)).*ones(1,K))./norm(OR_HRPT_CM'*inv(OR_HRPT_CM*OR_HRPT_CM'),'fro');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Precoding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T_theta_inbarDT_CM = zeros(Tc,NN );
        R_theta_inbarDT_CM = zeros(Tc,NN );
        OR_theta_inbarDT_CM = zeros(Tc,NN );
        T_theta_inbarDT_DM = zeros(Tc,NN );
        R_theta_inbarDT_DM = zeros(Tc,NN );
        for ic = 1:Tc  %% Time-varying DRIS
            T_betaDT_DM_temp = zeros(1,NN );
            T_phaDIRSDT_DM_Temp = zeros(1,NN );
            T_betaDT_CM_temp = zeros(1,NN );
            T_phaDIRSDT_CM_Temp = zeros(1,NN );

            R_betaDT_DM_temp = zeros(1,NN );
            R_phaDIRSDT_DM_Temp = zeros(1,NN );
            R_betaDT_CM_temp = zeros(1,NN );
            R_phaDIRSDT_CM_Temp = zeros(1,NN );
            
            OR_betaDT_CM_temp = zeros(1,NN );
            OR_phaDIRSDT_CM_Temp = zeros(1,NN );
            for ri = 1:length(T_phaDIRSDT_DM_Temp)
                iddele = randi(2^b,1,1)*randi(2^b,1,1) ;
                if iddele == 1
                    T_phaDIRSDT_DM_Temp(ri) = T_Omg(1);
                    T_betaDT_DM_temp(ri) = T_Amp(1);
                    T_phaDIRSDT_CM_Temp(ri) = T_Omg(1);
                    T_betaDT_CM_temp(ri) = sqrt(2)/2;
                    
                    R_phaDIRSDT_DM_Temp(ri) = R_Omg(1);
                    R_betaDT_DM_temp(ri) = R_Amp(1);
                    R_phaDIRSDT_CM_Temp(ri) = R_Omg(1);
                    R_betaDT_CM_temp(ri) = sqrt(2)/2;
                    
                    OR_phaDIRSDT_CM_Temp(ri) = R_Omg(1);
                    OR_betaDT_CM_temp(ri) = 1;
                else
                    T_phaDIRSDT_DM_Temp(ri) = T_Omg(2);
                    T_betaDT_DM_temp(ri) = T_Amp(2);
                    T_phaDIRSDT_CM_Temp(ri) = T_Omg(2);
                    T_betaDT_CM_temp(ri) = sqrt(2)/2;
                    
                    R_phaDIRSDT_DM_Temp(ri) = R_Omg(2);
                    R_betaDT_DM_temp(ri) = R_Amp(2);
                    R_phaDIRSDT_CM_Temp(ri) = R_Omg(2);
                    R_betaDT_CM_temp(ri) = sqrt(2)/2;
                    
                    OR_phaDIRSDT_CM_Temp(ri) = R_Omg(2);
                    OR_betaDT_CM_temp(ri) = 1;
                end
            end
            T_theta_inbarDT_CM(ic,:) = T_betaDT_CM_temp.*exp(1j.*(T_phaDIRSDT_CM_Temp));
            R_theta_inbarDT_CM(ic,:) = R_betaDT_CM_temp.*exp(1j.*(R_phaDIRSDT_CM_Temp));
            OR_theta_inbarDT_CM(ic,:) = OR_betaDT_CM_temp.*exp(1j.*(OR_phaDIRSDT_CM_Temp));
            T_theta_inbarDT_DM(ic,:) = T_betaDT_DM_temp.*exp(1j.*(T_phaDIRSDT_DM_Temp));
            R_theta_inbarDT_DM(ic,:) = R_betaDT_DM_temp.*exp(1j.*(R_phaDIRSDT_DM_Temp));
            %% General Fully Passive Jammer
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constant/Diff Amplitude %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            HD_DT_EachTime_CM = [ Hr(1:K/2,:)*diag( T_theta_inbarDT_CM(ic,:) )*Gr; Hr(K/2+1:K,:)*diag( R_theta_inbarDT_CM(ic,:) )*Gr ] + Hd;
            HD_DT_EachTime_DM = [ Hr(1:K/2,:)*diag( T_theta_inbarDT_DM(ic,:) )*Gr; Hr(K/2+1:K,:)*diag( R_theta_inbarDT_DM(ic,:) )*Gr ] + Hd;
            HD_DT_EachTime_CM_OR = [ 0.*Hr(1:K/2,:)*diag( T_theta_inbarDT_CM(ic,:) )*Gr; Hr(K/2+1:K,:)*diag( OR_theta_inbarDT_CM(ic,:) )*Gr ] + Hd;
            for kk = 1:K
                User_Sig_GFPJ_CM(kk,pv) = User_Sig_GFPJ_CM(kk,pv) + norm( HD_DT_EachTime_CM( kk,:)*W_CM(:,kk) ,'fro' )^2;
                User_Sig_GFPJ_DM(kk,pv) = User_Sig_GFPJ_DM(kk,pv) + norm( HD_DT_EachTime_DM( kk,:)*W_DM(:,kk),'fro' )^2;
                User_Sig_GFPJ_OR(kk,pv) = User_Sig_GFPJ_OR(kk,pv) + norm( HD_DT_EachTime_CM_OR( kk,:)*W_CM_OR(:,kk),'fro' )^2;
                Theo_User_Sig_GFPJ_CM(kk,pv) = Theo_User_Sig_GFPJ_CM(kk,pv) + 2*P(pv)*(N-K)/sum(1./ld)*10.^(-noise/10) + P(pv)/K*lg*lr(kk)*NN*10.^(-noise/10);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Theorems 1 and 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if kk<=K/2
                    Theo_User_Sig_GFPJ_DM(kk,pv) = Theo_User_Sig_GFPJ_DM(kk,pv) + P(pv)*(N-K)/sum(1./ld)*10.^(-noise/10) + P(pv)/K*lg*lr(kk)*NN*10.^(-noise/10)* 0.66 ;
                else
                    Theo_User_Sig_GFPJ_DM(kk,pv) = Theo_User_Sig_GFPJ_DM(kk,pv) + P(pv)*(N-K)/sum(1./ld)*10.^(-noise/10) + P(pv)/K*lg*lr(kk)*NN*10.^(-noise/10)* 0.34 ;
                end
                User_set = [1:K];
                User_set(User_set == kk) = [];
                Theo_InterUser_Inter_GFPJ_CM(kk,pv)  = Theo_InterUser_Inter_GFPJ_CM(kk,pv) + P(pv)/K*lg*sum(lr(User_set))*NN*10.^(-noise/10) + 2;
                if kk<=K/2
                    Theo_InterUser_Inter_GFPJ_DM(kk,pv)  = Theo_InterUser_Inter_GFPJ_DM(kk,pv) + P(pv)/K*lg*sum(lr(User_set))*NN*10.^(-noise/10)* 0.66 + 1;
                else
                    Theo_InterUser_Inter_GFPJ_DM(kk,pv)  = Theo_InterUser_Inter_GFPJ_DM(kk,pv) + P(pv)/K*lg*sum(lr(User_set))*NN*10.^(-noise/10)* 0.34 + 1;
                end
                HD_age_CM=0;
                HD_age_DM = 0;
                Hd_age_CM_OR = 0;
                for u = 1:length(User_set)
                    HD_age_CM = HD_age_CM + norm( ...
                        HD_DT_EachTime_CM( kk,:)* W_CM(:,User_set(u))...
                        ,'fro' )^2;
                    HD_age_DM = HD_age_DM +...
                        norm( HD_DT_EachTime_DM( kk,:)*W_DM(:,User_set(u)),'fro' )^2;
                    Hd_age_CM_OR   =  Hd_age_CM_OR +  norm( ...
                        HD_DT_EachTime_CM_OR( kk,:)*W_CM_OR(:,User_set(u)) ,...
                        'fro' )^2;
                end
                InterUser_Inter_GFPJ_CM(kk,pv) =  InterUser_Inter_GFPJ_CM(kk,pv) + HD_age_CM + 1;
                InterUser_Inter_GFPJ_DM(kk,pv) =  InterUser_Inter_GFPJ_DM(kk,pv)+HD_age_DM +1;
                InterUser_Inter_GFPJ_OR(kk,pv) =  InterUser_Inter_GFPJ_OR(kk,pv) + Hd_age_CM_OR + 1;
            end
        end
    end
end
%% Averaging Over Iterations
User_Sig_WO = User_Sig_WO./Niter;
InterUser_Inter_WO = InterUser_Inter_WO./Niter;
User_Sig_AJ = User_Sig_AJ./Niter;
InterUser_Inter_AJ = InterUser_Inter_AJ./Niter;
User_Sig_GFPJ_CM = User_Sig_GFPJ_CM./Niter./Tc;
InterUser_Inter_GFPJ_CM = InterUser_Inter_GFPJ_CM./Niter./Tc;
Theo_User_Sig_GFPJ_CM = Theo_User_Sig_GFPJ_CM./Niter./Tc;
Theo_InterUser_Inter_GFPJ_CM = Theo_InterUser_Inter_GFPJ_CM./Niter./Tc;
User_Sig_GFPJ_OR = User_Sig_GFPJ_OR./Niter./Tc;
InterUser_Inter_GFPJ_OR = InterUser_Inter_GFPJ_OR./Niter./Tc;
User_Sig_GFPJ_DM = User_Sig_GFPJ_DM./Niter./Tc;
InterUser_Inter_GFPJ_DM = InterUser_Inter_GFPJ_DM./Niter./Tc;
Theo_User_Sig_GFPJ_DM = Theo_User_Sig_GFPJ_DM./Niter./Tc;
Theo_InterUser_Inter_GFPJ_DM = Theo_InterUser_Inter_GFPJ_DM./Niter./Tc;

SumRate_WO_T = zeros(1,length(PdBm));
SumRate_AJ_T = zeros(1,length(PdBm));
SumRate_GFPJ_CM_T = zeros(1,length(PdBm));
Theo_SumRate_GFPJ_CM_T = zeros(1,length(PdBm));
SumRate_GFPJ_OR_T = zeros(1,length(PdBm));
SumRate_GFPJ_DM_T = zeros(1,length(PdBm));
Theo_SumRate_GFPJ_DM_T = zeros(1,length(PdBm));
SumRate_WO_R = zeros(1,length(PdBm));
SumRate_AJ_R = zeros(1,length(PdBm));
SumRate_GFPJ_CM_R = zeros(1,length(PdBm));
Theo_SumRate_GFPJ_CM_R = zeros(1,length(PdBm));
SumRate_GFPJ_OR_R = zeros(1,length(PdBm));
SumRate_GFPJ_DM_R = zeros(1,length(PdBm));
Theo_SumRate_GFPJ_DM_R = zeros(1,length(PdBm));
%% Sum Rate Calculation
% Loop through power levels and users to calculate the sum rate 
for pp = 1:length(PdBm)
    for k = 1:K
        if k<=K/2
            % Calculate sum rate for transmissive users (first half)
            SumRate_WO_T(pp) = SumRate_WO_T(pp) + log2( 1+ User_Sig_WO(k,pp)/InterUser_Inter_WO(k,pp) );
            SumRate_AJ_T(pp) = SumRate_AJ_T(pp) + log2( 1+ User_Sig_AJ(k,pp)/InterUser_Inter_AJ(k,pp) );
            SumRate_GFPJ_CM_T(pp) = SumRate_GFPJ_CM_T(pp) + log2( 1+ User_Sig_GFPJ_CM(k,pp)/InterUser_Inter_GFPJ_CM(k,pp) );
            Theo_SumRate_GFPJ_CM_T(pp) = Theo_SumRate_GFPJ_CM_T(pp) + log2( 1+ Theo_User_Sig_GFPJ_CM(k,pp)/Theo_InterUser_Inter_GFPJ_CM(k,pp) );
            SumRate_GFPJ_OR_T(pp) = SumRate_GFPJ_OR_T(pp) + log2( 1+ User_Sig_GFPJ_OR(k,pp)/InterUser_Inter_GFPJ_OR(k,pp) );
            SumRate_GFPJ_DM_T(pp) = SumRate_GFPJ_DM_T(pp) + log2( 1+ User_Sig_GFPJ_DM(k,pp)/InterUser_Inter_GFPJ_DM(k,pp) );
            Theo_SumRate_GFPJ_DM_T(pp) = Theo_SumRate_GFPJ_DM_T(pp) + log2( 1+ Theo_User_Sig_GFPJ_DM(k,pp)/Theo_InterUser_Inter_GFPJ_DM(k,pp) );
        else
            % Calculate sum rate for reflective users (second half)
            SumRate_WO_R(pp) = SumRate_WO_R(pp) + log2( 1+ User_Sig_WO(k,pp)/InterUser_Inter_WO(k,pp) );
            SumRate_AJ_R(pp) = SumRate_AJ_R(pp) + log2( 1+ User_Sig_AJ(k,pp)/InterUser_Inter_AJ(k,pp) );
            SumRate_GFPJ_CM_R(pp) = SumRate_GFPJ_CM_R(pp) + log2( 1+ User_Sig_GFPJ_CM(k,pp)/InterUser_Inter_GFPJ_CM(k,pp) );
            Theo_SumRate_GFPJ_CM_R(pp) = Theo_SumRate_GFPJ_CM_R(pp) + log2( 1+ Theo_User_Sig_GFPJ_CM(k,pp)/Theo_InterUser_Inter_GFPJ_CM(k,pp) );
            SumRate_GFPJ_OR_R(pp) = SumRate_GFPJ_OR_R(pp) + log2( 1+ User_Sig_GFPJ_OR(k,pp)/InterUser_Inter_GFPJ_OR(k,pp) );
            SumRate_GFPJ_DM_R(pp) = SumRate_GFPJ_DM_R(pp) + log2( 1+ User_Sig_GFPJ_DM(k,pp)/InterUser_Inter_GFPJ_DM(k,pp) );
            Theo_SumRate_GFPJ_DM_R(pp) = Theo_SumRate_GFPJ_DM_R(pp) + log2( 1+ Theo_User_Sig_GFPJ_DM(k,pp)/Theo_InterUser_Inter_GFPJ_DM(k,pp) );
        end
    end
end
%Plot simulation results
PdBmPerLU = PdBm; 
close all
figure;

plot(PdBmPerLU,SumRate_WO_T./K,'k-.>','LineWidth',2,'MarkerSize',8);hold on;

plot(PdBmPerLU,SumRate_GFPJ_CM_T./K,'rs','LineWidth',2,'MarkerSize',8);hold on;
plot(PdBmPerLU,Theo_SumRate_GFPJ_CM_T./K,'r-','LineWidth',2,'MarkerSize',8);hold on;

plot(PdBmPerLU,SumRate_GFPJ_DM_T./K,'bd','LineWidth',2,'MarkerSize',8);hold on;
plot(PdBmPerLU,Theo_SumRate_GFPJ_DM_T./K,'b-','LineWidth',2,'MarkerSize',8);hold on;

plot(PdBmPerLU,SumRate_GFPJ_OR_T./K,'m-.x','LineWidth',2,'MarkerSize',8);hold on;

plot(PdBmPerLU,SumRate_AJ_T./K,'c-.o','LineWidth',2,'MarkerSize',8);hold on;


grid on
axis([Pmin,Pmax,0,7.5])
xlabel('Transmit Power Per LU [dBm]');
ylabel('Refractive Rate Per LU [bits/symbol/user]');
legend('W/O Jamming','Proposed W/CA','Theorem 1','Proposed W/DA',...
    'Theorem 2','R-FPJ in [20]', 'AJ W/ 5 dBm','Location','Northwest');

figure;

plot(PdBmPerLU,SumRate_WO_R./K,'k-.>','LineWidth',2,'MarkerSize',8);hold on;

plot(PdBmPerLU,SumRate_GFPJ_CM_R./K,'rs','LineWidth',2,'MarkerSize',8);hold on;
plot(PdBmPerLU,Theo_SumRate_GFPJ_CM_R./K,'r-','LineWidth',2,'MarkerSize',8);hold on;

plot(PdBmPerLU,SumRate_GFPJ_DM_R./K,'bd','LineWidth',2,'MarkerSize',8);hold on;
plot(PdBmPerLU,Theo_SumRate_GFPJ_DM_R./K,'b-','LineWidth',2,'MarkerSize',8);hold on;

plot(PdBmPerLU,SumRate_GFPJ_OR_R./K,'m-.x','LineWidth',2,'MarkerSize',8);hold on;

plot(PdBmPerLU,SumRate_AJ_R./K,'c-.o','LineWidth',2,'MarkerSize',8);hold on;

grid on
axis([Pmin,Pmax,0,7.5])
xlabel('Transmit Power Per LU [dBm]');
ylabel('Reflective Rate Per LU [bits/symbol/user]');
legend('W/O Jamming','Proposed W/CA','Theorem 1','Proposed W/DA',...
    'Theorem 2','R-FPJ in [20]', 'AJ W/ 5 dBm','Location','Northwest');