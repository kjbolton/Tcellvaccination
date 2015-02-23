function [SEIRS_dot_new,Sout,SAout] = RHS_flow_stratified(time, iorder, m, params1, params, mo, SEIRS, itot, px, Sin, S2in, cycleA, cycleI, separateS)

% RHS_flow - flow in single 'loop' in RHS_7stage for model which accounts
% for reinfections
%
% multi-dimensional stratified model allowing for varying infectiousness
% (and potentially suceptibility) in different subsets of the population
%
% Arguments: m <int>: number of dimensions/stratifications in model
%            params1 <array of 7 arrays>: parameters defining flow 
%            mo <struct>: model options structure defined in MCMC*.m
%            SEIRS <1 X (16*m) double>: current states for loop in question
%            Sin <double X m>: additional flow rate of suceptibles from loop in
%                          a previous wave
%            Itot <double X m>: total current prevalence across all loops
%            cycleA <string>: signal for whether recovered A's are fed into
%                            S-state for this loop ('cycleA') or fed into 
%                            another loop ('NOcycleA')
%            cycleI <string>: signal for whether recovered I's are fed into
%                            S-state for this loop ('cycleI') or fed into 
%                            another loop ('NOcycleI')
% Output: SERIS_dot <1X double X m>: state flow for loop in question
%         Sout <double X m>: flow into new S-state due to people with temporary
%                        acquired immunity
%         SAout <double X m>: flow into (new?) S-state due to recovered
%                         asymptomatic infections


%SEIRS_dot = zeros(size(SEIRS));

if strcmp(mo.vacc_timing,'dynamic')
    [avTcell vBcell] = vacrate(params,SEIRS,px,time,mo);
else
    avTcell = [0 0];
    vBcell = 0;
end

beta = zeros(m);
Sout = zeros(m,1);
SAout = zeros(m,1);
if mo.symbolic == 1
  beta = sym(beta);
  Sout = sym(Sout);
  SAout = sym(SAout);
end

if (strcmp(separateS,'NOseparateS') && (sum(S2in) ~= 0.0))
    error('ERROR: Should not input flow into S2in if do not want to keep track of 2 susceptible pools in a single loop');
end

beta = params1.beta; %needs to be a m*m array

gamma = params1.gamma;
alpha_array = params1.alpha_array; %needs to be a m*1 array
nu = params1.nu; %needs to be a m*1 array
if mo.symbolic == 0
  phi = params1.phi/realpow(params1.chiTw,iorder);
else
  phi = params1.phi/(params1.chiTw.^iorder);
end

phiQ = params1.phiQ;
rho = params1.rho;
phix = params1.phix;


epsilonSQ =params1.epsilonSQ;
epsilonalphaQ =params1.epsilonalphaQ;

% each state needs to be m-dimensional

cS=linspace(1,m,m);
cE1=linspace(m+1,2*m,m);
cE2=linspace(2.*m+1,3*m,m);
cI=linspace(3.*m+1,4*m,m);
cA=linspace(4.*m+1,5*m,m);
cR=linspace(5.*m+1,6*m,m);
cT=linspace(6.*m+1,7*m,m);
cP=linspace(7.*m+1,8*m,m);
cCUM=linspace(8.*m+1,9*m,m);
cQ=linspace(9.*m+1,10*m,m);
cRA=linspace(10.*m+1,11*m,m);
cTA=linspace(11.*m+1,12*m,m);
cTQ=linspace(12.*m+1,13*m,m);
cT2=linspace(13.*m+1,14*m,m);
cTA2=linspace(14.*m+1,15*m,m);
cTQ2=linspace(15.*m+1,16*m,m);
cS2=linspace(16.*m+1,17*m,m);
cE1Q=linspace(17*m+1,18*m,m);
cE2Q=linspace(18*m+1,19*m,m);
cE1TQ=linspace(19*m+1,20*m,m);
cE2TQ=linspace(20*m+1,21*m,m);
cE1TQ2=linspace(21*m+1,22*m,m);
cE2TQ2=linspace(22*m+1,23*m,m);
cCUMA = linspace(23*m+1,24*m,m);
cD = linspace(24*m+1,25*m,m);

%initialise F# vectors

% infection and recovery flows

F1 = zeros(m,1);
F1A = zeros(m,1);
F1B = zeros(m,1);
F1C = zeros(m,1);
F2 = zeros(m,1);
F2A = zeros(m,1);
F2B = zeros(m,1);
F2C = zeros(m,1);
F3 = zeros(m,1);
F3A = zeros(m,1);
F3B = zeros(m,1);
F3C = zeros(m,1);
F4 = zeros(m,1);
F4A = zeros(m,1);
F4B = zeros(m,1);
F4C = zeros(m,1);
F5 = zeros(m,1);
F5A = zeros(m,1);
F6 = zeros(m,1);
F6A = zeros(m,1);
F7 = zeros(m,1);
F7A = zeros(m,1);
F7B = zeros(m,1);
F8 = zeros(m,1);
F9 = zeros(m,1);
F9A = zeros(m,1);
F9Ain = zeros(m,1);
F9Aout = zeros(m,1);
F10 = zeros(m,1);
F10A = zeros(m,1);
F10B = zeros(m,1);


% T cell vaccination flows

VT1 = zeros(m,1); % vaccination from S state
VT12 = zeros(m,1); % vaccination from S2 state
VT1A = zeros(m,1); % vaccination from TQ2 state
VT1B = zeros(m,1); % vaccination from TQ state
VT1C = zeros(m,1); % vaccination from Q state
VT2 = zeros(m,1); % vaccination from E1 state
VT2A = zeros(m,1); % vaccination from E1TQ2 state
VT2B = zeros(m,1); % vaccination from E1TQ state
VT2C = zeros(m,1); % vaccination from E1Q state
VT3 = zeros(m,1); % vaccination from E2 state
VT3A = zeros(m,1); % vaccination from E2TQ2 state
VT3B = zeros(m,1); % vaccination from E2TQ state
VT3C = zeros(m,1); % vaccination from E2Q state
VT5 = zeros(m,1); % vaccination from I state 
VT6 = zeros(m,1); % vaccination from A state 
VT7 = zeros(m,1); % vaccination from R state 
VT7A = zeros(m,1); % vaccination from RA state 
VT7B = zeros(m,1); % vaccination from TA state
VT7C = zeros(m,1); % vaccination from T state
VT9 = zeros(m,1); % vaccination from T2
VT9A = zeros(m,1); % vaccination from TA2
sumVT = zeros(m,1); % all vaccine flows
% loss of Tcell vaccination protection flows
LVT1 = zeros(m,1); % vaccination from S state
LVT12 = zeros(m,1); % vaccination from S2 state
LVT1A = zeros(m,1); % vaccination from TQ2 state
LVT1B = zeros(m,1); % vaccination from TQ state
LVT1C = zeros(m,1); % vaccination from Q state
LVT2 = zeros(m,1); % vaccination from E1 state
LVT2A = zeros(m,1); % vaccination from E1TQ2 state
LVT2B = zeros(m,1); % vaccination from E1TQ state
LVT2C = zeros(m,1); % vaccination from E1Q state
LVT3 = zeros(m,1); % vaccination from E2 state
LVT3A = zeros(m,1); % vaccination from E2TQ2 state
LVT3B = zeros(m,1); % vaccination from E2TQ state
LVT3C = zeros(m,1); % vaccination from E2Q state
LVT5 = zeros(m,1); % vaccination from I state 
LVT6 = zeros(m,1); % vaccination from A state 
LVT7 = zeros(m,1); % vaccination from R state 
LVT7A = zeros(m,1); % vaccination from RA state 
LVT7B = zeros(m,1); % vaccination from TA state
LVT7C = zeros(m,1); % vaccination from T state
LVT9 = zeros(m,1); % vaccination from T2
LVT9A = zeros(m,1); % vaccination from TA2
% B cell vaccination flows

VB1 = zeros(m,1); % vaccination from S state
VB12 = zeros(m,1); % vaccination from S2 state
VB1A = zeros(m,1); % vaccination from TQ2 state
VB1B = zeros(m,1); % vaccination from TQ state
VB1C = zeros(m,1); % vaccination from Q state
VB2 = zeros(m,1); % vaccination from E1 state
VB2A = zeros(m,1); % vaccination from E1TQ2 state
VB2B = zeros(m,1); % vaccination from E1TQ state
VB2C = zeros(m,1); % vaccination from E1Q state
VB3 = zeros(m,1); % vaccination from E2 state
VB3A = zeros(m,1); % vaccination from E2TQ2 state
VB3B = zeros(m,1); % vaccination from E2TQ state
VB3C = zeros(m,1); % vaccination from E2Q state
VB5 = zeros(m,1); % vaccination from I state 
VB6 = zeros(m,1); % vaccination from A state 
VB7 = zeros(m,1); % vaccination from R state 
VB7A = zeros(m,1); % vaccination from RA state 
VB7B = zeros(m,1); % vaccination from TA state
VB7C = zeros(m,1); % vaccination from T state
VB9 = zeros(m,1); % vaccination from T2
VB9A = zeros(m,1); % vaccination from TA2
sumVB = zeros(m,1); % all vaccine flows

if mo.symbolic == 1
  F1 = sym(F1);
  F1A = sym(F1A);
  F1B = sym(F1B);
  F1C = sym(F1C);
  F2 = sym(F2);
  F2A = sym(F2A);
  F2B = sym(F2B);
  F2C = sym(F2C);
  F3 = sym(F3);
  F3A = sym(F3A);
  F3B = sym(F3B);
  F3C = sym(F3C);
  F4 = sym(F4);
  F4A = sym(F4A);
  F4B = sym(F4B);
  F4C = sym(F4C);
  F5 = sym(F5);
  F6 = sym(F6);
  F7 = sym(F7);
  F7A = sym(F7A);
  F7B = sym(F7B);
  F8 = sym(F8);
  F9 = sym(F9);
  F9A = sym(F9A);
  F9Ain = sym(F9Ain);
  F9Aout = sym(F9Aout);
  F10 = sym(F10);
  F10A = sym(F10A);
  F10B = sym(F10B);
  
  % add vaccine flows here if want to convert to symbolic code
  
end

for i=1:m
    for j=1:m
        %fprintf(1,'F1 = %e beta %e itot %e S %e\n',F1(i),beta(i,j), itot(j), SEIRS(cS(i)));
        F1(i) = F1(i) + beta(i,j)*itot(j);
        F1A(i) = F1A(i) + beta(i,j)*(itot(j))*(epsilonSQ*SEIRS(cTQ2(i)));
        F1B(i) = F1B(i) + beta(i,j)*(itot(j))*(epsilonSQ*epsilonSQ*SEIRS(cTQ(i)));
        F1C(i) = F1C(i) + beta(i,j)*(itot(j))*(epsilonSQ*epsilonSQ*epsilonSQ*SEIRS(cQ(i)));
    end
    
%    fprintf(1,'itot %e F1 %e F1A %e F1B %e F1C %e Q %e E1Q %e S %e E1 %e\n',...
%        sum(itot), F1(i)*(SEIRS(cS(i)) + SEIRS(cS2(i))),F1A(i),F1B(i),F1C(i), SEIRS(cQ(i)), SEIRS(cE1Q(i)), SEIRS(cS(i)), SEIRS(cE1(i)));

       %fprintf(1,'11 beta %e itot %e S %e 12 beta %e itot %e S %e\n',beta(1,1),itot(1),SEIRS(cS(1)), beta(1,2),itot(2),SEIRS(cS(1)));
       %fprintf(1,'21 beta %e itot %e S %e 22 beta %e itot %e S %e\n',beta(2,1),itot(1),SEIRS(cS(2)), beta(2,2),itot(2),SEIRS(cS(2)));
    F2(i) = gamma*SEIRS(cE1(i));
    F2A(i) = gamma*SEIRS(cE1TQ2(i));
    F2B(i) = gamma*SEIRS(cE1TQ(i));
    F2C(i) = gamma*SEIRS(cE1Q(i));
    F3(i) = alpha_array(i)*gamma*SEIRS(cE2(i));
    F3A(i) = alpha_array(i)*epsilonalphaQ*gamma*SEIRS(cE2TQ2(i));
    F3B(i) = alpha_array(i)*epsilonalphaQ*epsilonalphaQ*gamma*SEIRS(cE2TQ(i));
    F3C(i) = alpha_array(i)*epsilonalphaQ*epsilonalphaQ*epsilonalphaQ*gamma*SEIRS(cE2Q(i));
    F4(i) = (1-alpha_array(i))*gamma*SEIRS(cE2(i));
    F4A(i) = (1-alpha_array(i)*epsilonalphaQ)*gamma*SEIRS(cE2TQ2(i));
    F4B(i) = (1-alpha_array(i)*epsilonalphaQ*epsilonalphaQ)*gamma*SEIRS(cE2TQ(i));
    F4C(i) = (1-alpha_array(i)*epsilonalphaQ*epsilonalphaQ*epsilonalphaQ)*gamma*SEIRS(cE2Q(i));
    % note that F3 plus F4 just add to gamma*E2 and do not depend on alpha
    F5(i) = nu(i)*SEIRS(cI(i));
    F6(i) = nu(i)*SEIRS(cA(i));
    
    if m==1 || m==2
        frac_vac = 1-erfc(avTcell(i));
       
        
        F5(i) = (1-frac_vac)*F5(i);
        F5A(i) = frac_vac*F5(i);
        F6(i) = (1-frac_vac)*F6(i);
        F6A(i) = frac_vac*F6(i);
    end
    
    if strcmp(mo.transient_xprotection,'transient_xprotection')
        F7(i) = (1-rho(i))*phix*SEIRS(cR(i));
        F8(i) = rho(i)*phix*SEIRS(cR(i));
    else
        F7(i) = (1-rho(i))*phi*SEIRS(cR(i));
        F8(i) = rho(i)*phi*SEIRS(cR(i));
    end
    F9(i) = phi*SEIRS(cT2(i));
    F10(i) = phiQ*SEIRS(cTQ2(i));
    F10A(i) = phiQ*SEIRS(cQ(i));
    F7A(i) = (1-rho(i))*phi*SEIRS(cRA(i));
    
    if strcmp(mo.waning,'reset')
        F10B(i) = 0;
        F7B(i) = 0;
        F7C(i) = 0;
    else
        F10B(i) = phiQ*SEIRS(cTQ(i));
        F7B(i) = phi*SEIRS(cTA(i));
        F7C(i) = phi*SEIRS(cT(i));
    end
    
    F8A(i) = rho(i)*phi*SEIRS(cRA(i));
    F9A(i) = phi*SEIRS(cTA2(i));
    
    Sout(i) = F9(i);
        
    if strcmp(cycleA,'cycleA')
        SAout(i) = 0;
        F9Ain(i) = F9A(i);
        F9Aout(i) = F9A(i);
    else
        % k: cycle recovered asymptomatic infections into another loop, and
        % remove from this loop
        SAout(i) = F9A(i);
        F9Ain(i) = 0;
        F9Aout(i) = F9A(i);
    end
    
    if (strcmp(mo.divert_between_strata,'divert') && (m==2))
        Sout(1) = sum(Sout);
        Sout(2) = 0.;
        SAout(1) = sum(SAout);
        SAout(2) = 0.;
        F9Ain(1) = sum(F9Ain);
        F9Ain(2) = 0.;
    end
    
    
    
    if strcmp(cycleI,'cycleI')
        % set Sin to be Sout
        Sin(i) = Sin(i) + Sout(i);
        Sout(i) = 0.;
    end
    
    
    % T cell vaccine flows
    
    if i==1 || i==2
        
        % flow out of unvaccinated states
        
        %set Vexposed=0 if vaccination can't alter infection experience of
        %exposed or infected hosts
        Vexposed = mo.Tcellvac_exposed;
        
        VT1(i) = avTcell(i)*SEIRS(cS(i));
        VT12(i) = avTcell(i)*SEIRS(cS2(i));
        VT1A(i) = avTcell(i)*SEIRS(cTQ2(i));% assume partial protection boosted
        VT1B(i) = avTcell(i)*SEIRS(cTQ(i));
        VT1C(i) = avTcell(i)*SEIRS(cQ(i));
        VT2(i) = avTcell(i)*SEIRS(cE1(i))*Vexposed; % assumes vaccine takes effect quickly enough to affect infection experience of exposed individuals
        VT2A(i) = avTcell(i)*SEIRS(cE1TQ2(i))*Vexposed;
        VT2B(i) = avTcell(i)*SEIRS(cE1TQ(i))*Vexposed;
        VT2C(i) = avTcell(i)*SEIRS(cE1Q(i))*Vexposed;
        VT3(i) = avTcell(i)*SEIRS(cE2(i))*Vexposed;
        VT3A(i) = avTcell(i)*SEIRS(cE2TQ2(i))*Vexposed;
        VT3B(i) = avTcell(i)*SEIRS(cE2TQ2(i))*Vexposed;
        VT3C(i) = avTcell(i)*SEIRS(cE2Q(i))*Vexposed;
        VT5(i) = avTcell(i)*SEIRS(cI(i))*Vexposed; % assumes vaccination while infected doesn't change course of infection
        VT6(i) = avTcell(i)*SEIRS(cA(i))*Vexposed;
        VT7(i) = avTcell(i)*SEIRS(cR(i)); % assumes vaccination with antibody protection doesn't change susceptibility (assumes vaccine waning rate shorter than longevity of sterilising protection)
        VT7A(i) = avTcell(i)*SEIRS(cRA(i));
        VT7B(i) = avTcell(i)*SEIRS(cTA(i)); % vaccination from TA state
        VT7C(i) = avTcell(i)*SEIRS(cT(i)); % vaccination from T state
        VT9(i) = avTcell(i)*SEIRS(cT2(i)); % vaccination from T2 
        VT9A(i) = avTcell(i)*SEIRS(cTA2(i)); % vaccination from TA2 
        
        sumVT(i) = VT1(i) + VT12(i) + VT1A(i) + VT1B(i) + VT1C(i) + VT2(i) + VT2A(i) + VT2B(i) + VT2C(i) + VT3(i) + VT3A(i) + ...
            VT3B(i) + VT3C(i) + VT5(i) + VT6(i) + VT7(i) + VT7A(i) + VT7B(i) + VT7C(i) + VT9(i) + VT9A(i);
        
        
     
        
        %fprintf('i %i SUMVT = %e vTcell = %e \n',i, sumVT(i)/params.N,avTcell(i)*px(i));
        if (mo.Tcellvac_exposed)
            if (sumVT(i)/params.N<0.95*avTcell(i)*px(i))
                fprintf('i = %i SUMVT = %e vTcell = %e \n',i, sumVT(i)/params.N,avTcell(i)*px(i));
                avTcell
                px
                fprintf('WARNING: inaccuracy in conserving total specified vaccine flow rate\n',time)
            end
        end
       % flow out of vaccinated states (Loss of Vaccine T cell protection)
       % currently assumes loss of vaccine takes host to like unvaccinated
       % state
       phiT = params.phi_Tcell;
       
       LVT1(i) = phiT*SEIRS(cS(i+2));
       LVT12(i) = phiT*SEIRS(cS2(i+2));
       LVT1A(i) = phiT*SEIRS(cTQ2(i+2));
       LVT1B(i) = phiT*SEIRS(cTQ(i+2));
       LVT1C(i) = phiT*SEIRS(cQ(i+2));
       LVT2(i) = phiT*SEIRS(cE1(i+2)); 
       LVT2A(i) = phiT*SEIRS(cE1TQ2(i+2));
       LVT2B(i) = phiT*SEIRS(cE1TQ(i+2));
       LVT2C(i) = phiT*SEIRS(cE1Q(i+2));
       LVT3(i) = phiT*SEIRS(cE2(i+2));
       LVT3A(i) = phiT*SEIRS(cE2TQ2(i+2));
       LVT3B(i) = phiT*SEIRS(cE2TQ2(i+2));
       LVT3C(i) = phiT*SEIRS(cE2Q(i+2));
       LVT5(i) = phiT*SEIRS(cI(i+2)); 
       LVT6(i) = phiT*SEIRS(cA(i+2));
       LVT7(i) = phiT*SEIRS(cI(i+2)); 
       LVT7A(i) = phiT*SEIRS(cA(i+2));
       LVT7B(i) = phiT*SEIRS(cTA(i+2)); 
       LVT7C(i) = phiT*SEIRS(cT(i+2)); 
       LVT9(i) = phiT*SEIRS(cT2(i+2)); 
       LVT9A(i) = phiT*SEIRS(cTA2(i+2)); 
        
    else
        
        % flow into vaccinated states is opposite to that out of
        % unvaccinated states (except for recovering infected or
        % asymptomatically infected hosts, who need to be diverted
        % separately)
        
        VT1(i) = -VT1(i-2);
        VT12(i) = -VT12(i-2);
        VT1A(i) = -VT1A(i-2);% assume partial protection boosted
        VT1B(i) = -VT1B(i-2);
        VT1C(i) = -VT1C(i-2);
        VT2(i) = -VT2(i); % assumes vaccine takes effect quickly enough to affect infection experience of exposed individuals
        VT2A(i) = -VT2A(i);
        VT2B(i) = -VT2B(i);
        VT2C(i) = -VT2C(i);
        VT3(i) = -VT3(i-2);
        VT3A(i) = -VT3A(i-2);
        VT3B(i) = -VT3B(i-2);
        VT3C(i) = -VT3C(i-2);
        VT5(i) = -VT5(i-2); % assumes vaccination while infected doesn't change course of infection
        VT6(i) = -VT6(i-2);
        VT7(i) = -VT7(i-2); % assumes vaccination with antibody protection doesn't change susceptibility (assumes vaccine waning rate shorter than longevity of sterilising protection)
        VT7A(i) = -VT7A(i-2);
        VT7B(i) = -VT7B(i-2); % vaccination from TA state
        VT7C(i) = -VT7C(i-2); % vaccination from T state
        VT9(i) = -VT9(i-2); % vaccination from T2 
        VT9A(i) = -VT9A(i-2); % vaccination from TA2 
        
        sumVT(i) = -sumVT(i-2);
        
        
        % loss of vaccine protection from vaccinated strata
        LVT1(i) = -LVT1(i-2);
        LVT12(i) = -LVT12(i-2);
        LVT1A(i) = -LVT1A(i-2);
        LVT1B(i) = -LVT1B(i-2);
        LVT1C(i) = -LVT1C(i-2);
        LVT2(i) = -LVT2(i); 
        LVT2A(i) = -LVT2A(i);
        LVT2B(i) = -LVT2B(i);
        LVT2C(i) = -LVT2C(i);
        LVT3(i) = -LVT3(i-2);
        LVT3A(i) = -LVT3A(i-2);
        LVT3B(i) = -LVT3B(i-2);
        LVT3C(i) = -LVT3C(i-2);
        LVT5(i) = -LVT5(i-2); 
        LVT6(i) = -LVT6(i-2);
        LVT7(i) = -LVT7(i-2); 
        LVT7A(i) = -LVT7A(i-2);
        LVT7B(i) = -LVT7B(i-2); 
        LVT7C(i) = -LVT7C(i-2); 
        LVT9(i) = -LVT9(i-2); 
        LVT9A(i) = -LVT9A(i-2); 
        
        
    end
    
   
    % B cell vaccine flows
       
    VB1(i) = vBcell*SEIRS(cS(i));
    VB12(i) = vBcell*SEIRS(cS2(i));
    VB1A(i) = vBcell*SEIRS(cTQ2(i));% assume partial protection boosted 
    VB1B(i) = vBcell*SEIRS(cTQ(i)); 
    VB1C(i) = vBcell*SEIRS(cQ(i));
    VB2(i) = 0.;
    VB2A(i) = 0.;
    VB2B(i) = 0.;
    VB2C(i) = 0.;
    VB3(i) = 0.;
    VB3A(i) = 0.;
    VB3B(i) = 0.;
    VB3C(i) = 0.;
    VB5(i) = 0.; % assumes vaccination while infected doesn't change course of infection
    VB6(i) = 0.;
    VB7(i) = 0.; % assumes vaccination with antibody protection doesn't change susceptibility (assumes vaccine waning rate shorter than longevity of sterilising protection)
    VB7A(i) = 0.;
    VB7B(i) = vBcell*SEIRS(cTA(i)); 
    VB7C(i) = vBcell*SEIRS(cT(i)); 
    VB9(i) = vBcell*SEIRS(cT2(i)); 
    VB9A(i) = vBcell*SEIRS(cTA2(i));
    
    sumVB(i) = VB1(i) + VB12(i) + VB1A(i) + VB1B(i) + VB1C(i) + VB2(i) + VB2A(i) + VB2B(i) + VB2C(i) + VB3(i) + VB3A(i) + ...
        VB3B(i) + VB3C(i) + VB5(i) + VB6(i) + VB7(i) + VB7A(i) + VB7B(i) + VB7C(i) + VB9(i) + VB9A(i);
end


% for SINGLE stage code only
if strcmp(mo.reinfection,'noreinfection')
    if strcmp(mo.divert_between_strata,'divert')
        % keep recovered infected individuals cycling in same loop
        % also divert T2(1)'s to S(2)
        if (m==2)
            Sin(1) = Sin(1) + sum(Sout);
        end
        if (m==4)
            Sin(1) = Sin(1) + sum(Sout(1:2));
            Sin(3) = Sin(3) + sum(Sout(3:4));
        end
    else
         % keep recovered infected individuals cycling in same loop
        Sin = Sin + Sout;
    end
end


% divert recovered vaccinated hosts (T cell vaccine) to stratum 3 or 4
% (vaccinated)
% diversion to experienced strata for naiive hosts is done seperately below
% when the model option 'divert' is specified
if (m==4) && strcmp(mo.vaccination,'none')==0
    % really want to do this if vaccination is not 'none' or 'Bcell'...
    F5A(3) = F5A(1) + F5A(3);
    F5A(4) = F5A(2) + F5A(4);
    F5A(1:2) = zeros(2,1);
    
    F6A(3) = F6A(1) + F6A(3);
    F6A(4) = F6A(2) + F6A(4);
    F6A(1:2) = zeros(2,1);
    
end


for i = 1:m
    if strcmp(separateS,'separateS')
        % equal (relative) proportion of susceptibles drawn from S pool and S2 pool
        % pool, recovering susceptibles fed into S2, incoming susceptibles
        % still fed into S
        % feed people from Q into S
        %fprintf(1,'S = %e S2 = %e\n',SEIRS(cS(i)), SEIRS(cS2(i)));
        if ((SEIRS(cS(i)) + SEIRS(cS2(i))) == 0.)
            % no susceptibles in either pool, can't infect new people
            SEIRS_dot(cS(i)) = F10(i) + Sin(i) - VT1(i) - VB1(i);
            SEIRS_dot(cS2(i)) = S2in(i) + F9Ain(i) - VT12(i) - VB12(i);
        else
            SEIRS_dot(cS(i)) = -F1(i)*SEIRS(cS(i)) + F10(i) + Sin(i) - VT1(i) - VB1(i);
            % feed people who have already been infected in this stage into S2
            % (Sin contains symptomatic infections, F9Ain contains people who
            % have had asymptomatic infections)
            SEIRS_dot(cS2(i)) =  -F1(i)*SEIRS(cS2(i)) + S2in(i) + F9Ain(i) - VT12(i) - VB12(i);
        end
    else
        %fprintf(1,'F9Ain = %e F10 = %e F1 = %e Sin = %e\n',F9Ain(i),F10(i),F1(i),Sin(i));
        % no susceptibles fed into S2 pool
        SEIRS_dot(cS(i)) = -F1(i)*(SEIRS(cS(i)) + SEIRS(cS2(i)))  + F10(i) + Sin(i) + F9Ain(i) - VT1(i) - VB1(i);
        SEIRS_dot(cS2(i)) = -VT12(i) - VB12(i);
    end
    %fprintf(1,'i = %i F1 %e F10 %e Sin %e F9Ain %e\n',i,F1(i),F10(i),Sin(i),F9Ain(i));
    SEIRS_dot(cE1(i)) = F1(i)*(SEIRS(cS(i)) + SEIRS(cS2(i))) - F2(i) - VT2(i) - VB2(i) + LVT2(i);
    SEIRS_dot(cE2(i)) = F2(i) - F3(i) - F4(i) - VT3(i) - VB3(i) + LVT3(i);
    SEIRS_dot(cI(i)) = F3(i) + F3A(i) + F3B(i) + F3C(i)- F5(i) - VT5(i) - VB5(i) + LVT5(i);
    SEIRS_dot(cA(i)) = F4(i) + F4A(i) + F4B(i) + F4C(i)- F6(i) - VT6(i) - VB6(i) + LVT6(i);
    SEIRS_dot(cR(i)) = F5(i) + F5A(i) - F7(i) - F8(i) - VT7(i) - VB7(i) + LVT7(i);
    SEIRS_dot(cT(i)) = F7(i) - F7C(i) - VT7C(i) - VB7C(i) + LVT7C(i);
    SEIRS_dot(cP(i)) = F8(i) + F8A(i); % vaccine does nothing to those in P box
    SEIRS_dot(cCUM(i)) = F3(i) + F3A(i) + F3B(i) + F3C(i);
    SEIRS_dot(cQ(i)) = -F10A(i) - F1C(i) - VT1C(i) - VB1C(i) + LVT1C(i);
    SEIRS_dot(cRA(i)) = F6(i) + F6A(i) - F7A(i) - F8A(i) - VT7A(i) -VB7A(i) + sumVB(i) + LVT7A(i);
    SEIRS_dot(cTA(i)) = F7A(i) - F7B(i) - VT7B(i) - VB7B(i) + LVT7B(i);
    SEIRS_dot(cTQ(i)) = F10A(i) - F10B(i) - F1B(i) - VT1B(i) - VB1B(i) + LVT1B(i);
    SEIRS_dot(cT2(i)) = F7C(i) - F9(i) - VT9(i) - VB9(i) + LVT9(i);
    SEIRS_dot(cTA2(i)) = F7B(i) - F9Aout(i) - VT9A(i) - VB9A(i) + LVT9A(i);
    SEIRS_dot(cTQ2(i)) = F10B(i) - F10(i) - F1A(i) - VT1A(i) - VB1A(i) + LVT1A(i);
    SEIRS_dot(cE1TQ2(i)) = F1A(i) - F2A(i) - VT2A(i) - VB2A(i) + LVT2A(i);
    SEIRS_dot(cE1TQ(i)) = F1B(i) - F2B(i) - VT2B(i) - VB2B(i) + LVT2B(i);
    SEIRS_dot(cE1Q(i)) = F1C(i) - F2C(i) - VT2C(i) - VB2C(i) + LVT2C(i);
    SEIRS_dot(cE2TQ2(i)) = F2A(i) - F3A(i) - F4A(i) - VT3A(i) - VB3A(i) + LVT3A(i);
    SEIRS_dot(cE2TQ(i)) = F2B(i) - F3B(i) - F4B(i) - VT3B(i) - VB3B(i) + LVT3B(i);
    SEIRS_dot(cE2Q(i)) = F2C(i) - F3C(i) - F4C(i) - VT3C(i) - VB3C(i) + LVT3C(i);
    SEIRS_dot(cCUMA(i)) = F4(i) + F4A(i) + F4B(i) + F4C(i);
    SEIRS_dot(cD(i)) = 0.;
end

% note indices here not related to strata
SEIRS_dot(cD(1)) = params1.phi_drift;
SEIRS_dot(cD(2)) = params1.phi_koelle;

% if want to divert infected naives to stratum 1 need to alter flow into
% cS2 - directing it all to stratum 1

if strcmp(mo.divert_between_strata,'divert')
    if (m==2)
        SEIRS_dot(cS2(1)) = SEIRS_dot(cS2(1)) + SEIRS_dot(cS2(2));
        SEIRS_dot(cS2(2)) = 0.;
    elseif (m==4)
        SEIRS_dot(cS2(1)) = SEIRS_dot(cS2(1)) + SEIRS_dot(cS2(2));
        SEIRS_dot(cS2(2)) = 0.;
        SEIRS_dot(cS2(3)) = SEIRS_dot(cS2(3)) + SEIRS_dot(cS2(4));
        SEIRS_dot(cS2(4)) = 0.;
    else
        error('havent specified what to do with option divert for this number of strata');
    end
end




SEIRS_dot_new = SEIRS_dot';


if (isnan(sum(SEIRS_dot_new)))
    SEIRS
    SEIRS_dot_new
    error('ERROR: at least one state is NaN');
end


%fprintf(1,'size SAout RHS_flow_stratified %i %i\n',size(SAout));
%fprintf(1,'size SEIRS_dot in RHS_flow_stratified %i %i\n',size(SEIRS_dot_new));


end
