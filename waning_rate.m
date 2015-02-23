function [phi phiQ stoteff stoteffalpha E1tot phi_drift phi_koelle SRtot Stot] = waning_rate(params,SEIRS,mo,cd,time_since_seeding)

global gthreshtime

% phi: waning rate of acquired protection between consecutive model states (inverse days)
% phiQ: waning rate of prior protection between consecutive model states (inverse days)
% stoteff: vector of effective suscepibility of population in each stratum (number of
% hosts)

if strcmp(mo.transient_xprotection,'transient_xprotection');
    nstate_waning = 2; 
else
    nstate_waning = 3;% 3 without transient cross protective immunity between R and T
end


nstate_waningQ = 3;
if strcmp(mo.protected_population_state,'TQ')
    % use this if fully protected population is put into box TQ, in
    % MCMCinit
    nstate_waningQ = 2;
end
m = mo.stratifications; % number of stratifications in model (i.e. allowing for heterogeneity in popultion), have only used two.
ns = mo.nstates; % number of states per strata per loop (15 independent states plus cumulative infections)
nmultiple = mo.nmultiple; % number of sub loops per stage (require one for primary infection and two for secondary infections).
stgmax = mo.nstgs; % number of stages required to compare with reinfection data (2*2*2 - 1 (for people never infected))
nloops = stgmax*m*nmultiple; % total number of SEIR loops in model

cS = [];
for j = 1:stgmax*nmultiple; % add two to seperate people who have one versus two or more infections
    cS = [cS linspace((j-1)*ns*m+1,((j-1)*ns+1)*m,m)];
end

cE1=cS+m;
cE2=cE1+m;
cI=cE2+m;
cA=cI+m;
cR=cA+m;
cT=cR+m;
cP=cT+m;
cCUM=cP+m;
cQ=cCUM+m;
cRA=cQ+m;
cTA=cRA+m;
cTQ=cTA+m;
cT2=cTQ+m;
cTA2=cT2+m;
cTQ2=cTA2+m;
cS2 =cTQ2+m;
cE1Q = cS2 + m;
cE2Q = cE1Q+m;
cE1TQ = cE2Q +m;
cE2TQ = cE1TQ + m;
cE1TQ2 = cE2TQ + m;
cE2TQ2 = cE1TQ2 + m;
cCUMA = cE2TQ2 + m;

itot = zeros(m,1);
stot = zeros(m,1);
srtot = zeros(m,1);
stoteff = zeros(m,1);
stoteffalpha = zeros(m,1);
E1tot = zeros(m,1);
icumtot = zeros(m,1);
icumtot_wave = zeros(m,1);
alpha_array = ones(m,1);
s_array = ones(m,1);
i_array = ones(m,1);
if mo.symbolic == 1
  itot = sym(itot);
  stot = sym(stot);
  srtot = sym(srtot);
  stoteff = sym(stoteff);
  stoteffalpha = sym(stoteffalpha);
  E1tot = sym(E1tot);
  icumtot = sym(icumtot);
  icumtot_wave = sym(icumtot_wave);
  s_array = sym(s_array);
  i_array = sym(i_array);
  alpha_array = sym(alpha_array);
end


% construct arrays with scaled infectiousness, susceptibility, symptomatic
% proportion and infectious period for experienced and naiive strata

if m==2
    i_array(1) = params.epsilonI;
    s_array(1) = params.epsilonS;
    alpha_array(1) =  alpha_array(1)*params.epsilonALPHA;
elseif m==4
    
    if strcmp(mo.drift_rate,'percapita')
        % leave i_array = 1
    elseif strcmp(mo.drift_rate,'propinfectiousness')
        % scale drift rate with total infectiousness of each case
    
        if strcmp(mo.vacc_action,'boost')
            epsilonIE = params.fIE*params.epsilonI;
        elseif strcmp(mo.vacc_action,'sat')
            epsilonIE = min(params.fIE,params.epsilonI);
        end
        epsilonIN = params.fIN;
    
        i_array(1) = params.epsilonI;
        i_array(2) = 1.;
        i_array(3) = epsilonIE;
        i_array(4) = epsilonIN;
    elseif strcmp(mo.drift_rate,'CTLdependent')
        % hosts with elevated CTLs can not spit out an antigenic variant,
        % others can
        i_array(3) = 0; % vaccinated adults exert no selection on virus and don't contribute to drift
    end
    
end


for i =1:m
    for j = 1:stgmax*nmultiple
        % sum over loops and subloops for each strata
        % total number infected
        itot(i) = itot(i) + (SEIRS(cI((j-1)*m+i)) + SEIRS(cA((j-1)*m+i)))*i_array(i);
        %total number susceptible or recovered without permanent protection
        srtot(i) = srtot(i) + s_array(i)*(SEIRS(cS((j-1)*m+i)) + SEIRS(cS2((j-1)*m+i)) + SEIRS(cT((j-1)*m+i)) + SEIRS(cT2((j-1)*m+i))...
            +  SEIRS(cTA((j-1)*m+i)) + SEIRS(cTA2((j-1)*m+i))) + SEIRS(cR((j-1)*m+i)) + SEIRS(cRA((j-1)*m+i));
        stot(i) = stot(i) + s_array(i)*(SEIRS(cS((j-1)*m+i)) + SEIRS(cS2((j-1)*m+i)));
        if (1==0)
            stot(i) = stot(i) + s_array(i)*(SEIRS(cQ((j-1)*m+i))*power(params.epsilonSQ,3)  + SEIRS(cTQ((j-1)*m+i))*power(params.epsilonSQ,2) +  SEIRS(cTQ((j-1)*m+i))*params.epsilonSQ);
            srtot(i) = srtot(i) + s_array(i)*(SEIRS(cQ((j-1)*m+i))*power(params.epsilonSQ,3)  + SEIRS(cTQ((j-1)*m+i))*power(params.epsilonSQ,2) +  SEIRS(cTQ((j-1)*m+i))*params.epsilonSQ);
        end
        
        %cumulative number infected in total (symptomatic only)
        icumtot(i) = icumtot(i) + SEIRS(cCUM((j-1)*m+i)) + SEIRS(cCUMA((j-1)*m+i));;
        % total effective susceptibility per stratum 
        stoteff(i) = stoteff(i) + (SEIRS(cS((j-1)*m+i)) + SEIRS(cS2((j-1)*m+i)) + SEIRS(cQ((j-1)*m+i))*power(params.epsilonSQ,3) +...
            SEIRS(cTQ((j-1)*m+i))*power(params.epsilonSQ,2) +  SEIRS(cTQ((j-1)*m+i))*params.epsilonSQ);
        
        % effective susceptibility modulated by symptomatic proportion for
        % calculation of Reff
        stoteffalpha(i) = stoteffalpha(i) + (SEIRS(cS((j-1)*m+i)) + SEIRS(cS2((j-1)*m+i))... 
        + SEIRS(cQ((j-1)*m+i))*power(params.epsilonSQ*params.epsilonalphaQ,3) +...
            SEIRS(cTQ((j-1)*m+i))*power(params.epsilonSQ*params.epsilonalphaQ,2) +  SEIRS(cTQ((j-1)*m+i))*params.epsilonSQ*params.epsilonalphaQ);
        
        E1tot(i) = E1tot(i) + SEIRS(cE1((j-1)*m+i)) + SEIRS(cE1TQ((j-1)*m+i)) + SEIRS(cE1TQ2((j-1)*m+i)) + SEIRS(cE1Q((j-1)*m+i));
        
    end
    
  
    loops = [1];
    R0 = params.R01;
    alpha = params.alpha1;
    alpha_array = alpha_array*alpha;
    %calculate cumulative number infected in current wave only (include
    %asymptomatic infections and weight by infectiousness)
    icumtot_wave(i) = (sum(SEIRS(cCUM((loops(:)-1)*m+i))) + sum(SEIRS(cCUMA((loops(:)-1)*m+i))))*i_array(i);
end


% default values for waning = 'normal'

phi = nstate_waning/exp(params.lgTw); % 3-state waning so phi = 3/Tw
phiQ = nstate_waningQ/exp(params.lgTwQ);


if strcmp(mo.waning,'mutwan')
    % Waning rate is assumed to increase due to viral drift
    
    %calculate total prevalence
    Iprev = sum(itot);
    %calculate total number of susceptible and temporarily recovered hosts
    SRtot = sum(srtot);
    Stot = sum(stot);
    %Calculate cumulative number of infections in current wave
    Icumtot_wave = sum(icumtot_wave);
    
    if gthreshtime < 0
        if Icumtot_wave < params.threshARkoelle*params.N
            % reference time set to 0, don't start counting cluster
            % lifetime
            threshtime = inf;
            reftime = 0;
        else
            reftime = time_since_seeding;
            threshtime = reftime;
      
            gthreshtime = threshtime;
        end
    else
        reftime = gthreshtime;
        threshtime = gthreshtime;
    end
    
    phi_drift = Iprev*R0/params.Ti*SRtot/(Stot+1)*Icumtot_wave/cd.N/cd.N;
    
    if strcmp(mo.koelle_drift,'cont')
    phi_koelle = Iprev/params.N*R0/params.Ti*power(SRtot/(Stot+1),1)*power((time_since_seeding-reftime)/params.lambda_koelle,params.kw_koelle-1)*heaviside(time_since_seeding-threshtime);
    elseif strcmp(mo.koelle_drift,'punc')
        phi_koelle = Iprev/params.N*R0/params.Ti*power((time_since_seeding-reftime)/params.lambda_koelle,params.kw_koelle-1)*heaviside(time_since_seeding-threshtime);
    else
        error('This type of estimate for Koelle-like drift is not recognised')        
    end
    if strcmp(mo.mutwan_form,'incprev')
        % all drift proportional to prevalence
        if strcmp(mo.mutwan_type,'icum')
            % form motivated by Boni et al.
            phi = nstate_waning*Iprev/(cd.N/100)*(1./exp(params.lgTw) + exp(params.kwan)*R0/params.Ti*SRtot/(Stot+1)*Icumtot_wave/cd.N);
            phiQ = nstate_waningQ*Iprev/(cd.N/100)*(1./exp(params.lgTwQ) + exp(params.kwanQ)*R0/params.Ti*SRtot/(Stot+1)*Icumtot_wave/cd.N);
        elseif strcmp(mo.mutwan_form,'iprev')
            % drift proportional to current prevalence only
            phi = nstate_waning*Iprev/(cd.N/100)*(1./exp(params.lgTw) + exp(params.kwan));
            phiQ = nstate_waningQ*Iprev/(cd.N/100)*(1./exp(params.lgTwQ) + exp(params.kwanQ));
        end
    elseif strcmp(mo.mutwan_form,'neuconst')
         % constant drift term independent of prevalence
        if strcmp(mo.mutwan_type,'icum')
            % form motivated by Boni et al.
            phi = nstate_waning*(1./exp(params.lgTw) + exp(params.kwan)*heaviside(Iprev)*R0/params.Ti*SRtot/(Stot+1)*Icumtot_wave/cd.N);
            phiQ = nstate_waningQ*(1./exp(params.lgTwQ) + exp(params.kwanQ)*heaviside(Iprev)*R0/params.Ti*SRtot/(Stot+1)*Icumtot_wave/cd.N);
        elseif strcmp(mo.mutwan_type,'iprev')
            % drift rate proportional to current prevalence only
            phi = nstate_waning*(1./exp(params.lgTw) + Iprev/(cd.N/100)*exp(params.kwan));
            phiQ = nstate_waningQ*(1./exp(params.lgTwQ) + Iprev/(cd.N/100)*exp(params.kwanQ));
        end
    else
        error('ERROR: have not chosen valid mutwan_form');
    end
    
   
end
