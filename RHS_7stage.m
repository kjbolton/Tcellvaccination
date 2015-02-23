function [SEIRS_dot requiredEqns] = RHS_7stage(time, SEIRS, params, N, mo, oo, cd, df,bf_array)
% RHS_SINGLEstage - the RHS of the single stage model (ODE specification).
% This model assumes that each host can only contract one symptomatic infection
% per wave. Recovered asymptomatic infections are cycled back to the
% suceptible pool of the same SEIRS loop until the next wave, when they are
% siphoned into the appropriate succeeding loop. Recovered symptomatic infections are held in
% the suceptible pool of the succeeding loop until the current wave period
% ends.
%
% Arguments: params <1xnum_params double>: model parameters
%            N <int>: Population size
%            inc_waning <48x1 double>: Incidence data used in waning-prop-inc model
%            mo <struct>: (see MCMC_*run.m for details)
%            oo <struct>: (see MCMC_*run.m for details)
%            cd <struct>: (see MCMC_*run.m for details)
%            time <double>: Simulation time
%            current_state <1xnum_states double>: Current state
%
% Output: Model state flows

m = mo.stratifications; % number of stratifications in model (i.e. allowing for heterogeneity in popultion), have only used two.
ns = mo.nstates; % number of states per strata per loop (15 independent states plus cumulative infections)
nmultiple = mo.nmultiple; % number of sub loops per stage (require one for primary infection and two for secondary infections).
stgmax = mo.nstgs; % number of stages required to compare with reinfection data (2*2*2 - 1 (for people never infected))
nloops = stgmax*m*nmultiple; % total number of SEIR loops in model



%assign an index to each state, keeping indices for like states (i.e. all
%susceptible states) in an array labelled with state name

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
cD = cCUMA + m;


% define array to house flows between all states
SEIRS_dot = zeros(1,stgmax*m*nmultiple*ns);
if mo.symbolic == 1
  SEIRS_dot = sym(SEIRS_dot);
end



% Calculate total instantaneous incidence across all loops in 7 stage model
% this is used both to calculate infection rate and assess global dieout
% condition below
itot = zeros(m,1);
px = zeros(m,1);
qstates = zeros(m,1);
alpha_array = ones(m,1);
nu_array = ones(m,1);
if mo.symbolic == 1
  itot = sym(itot);
  px = sym(px);
  qstates = sym(qstates);
  alpha_array = sym(alpha_array);
  nu_array = sym(nu_array);
end


% ODE equation parameters
gamma = 2/params.Te; % 2-state latent so gamma = 2/Te
nu = zeros(length(params.Ti),1);
if mo.symbolic == 1
  nu = sym(nu);
end
for i=1:length(params.Ti)
    nu(i) = 1/params.Ti(i); % 1-state infectious so nu = 1/Ti
end

% only modelling single epidemic wave
d12 = mo.ndays;
d23 = mo.ndays;

if mo.symbolic == 0
  if (time < d12)
    alpha = params.alpha1;
    R0 = params.R01;
    beta_factor = bf_array(1);
  elseif (time >= d12) && (time < d23)
    alpha = params.alpha2;
    R0 = params.R02;
    beta_factor = bf_array(2);
  else % time > d23
    alpha = params.alpha3;
    R0 = params.R03;
    beta_factor = bf_array(3);
  end
else
  if mo.wave_number == 0
    alpha = params.alpha1;
    R0 = params.R01;
    beta_factor = bf_array(1);
  elseif mo.wave_number == 1
    alpha = params.alpha2;
    R0 = params.R02;
    beta_factor = bf_array(2);
  elseif mo.wave_number == 2
    alpha = params.alpha3;
    R0 = params.R03;
    beta_factor = bf_array(3);
  else
    fprintf(1,'\nERROR: mo.wave_number not set correctly.\n')
    quit;
  end
end


% construct arrays with scaled infectiousness, susceptibility, symptomatic
% proportion and infectious period for experienced and naiive strata

alpha_array = alpha_array*alpha;
nu_array = nu_array*nu;

if (m==2)
    alpha_array(1) =  alpha_array(1)*params.epsilonALPHA; 
    nu_array(1) = nu_array(1)*params.epsilonNU;
end


for i =1:m
    for j = 1:stgmax*nmultiple
        % sum over loops and subloops for each strata
        % total number infected
        index = (j-1)*m+i;
        % modify by relative infectiousness epsilonI for experienced
        % stratum (this is done in matrix epsilon used in
        % RHS_flow_stratified.m)
        itot(i) = itot(i) + SEIRS(cI((j-1)*m+i)) + SEIRS(cA((j-1)*m+i));
        %total number susceptible or recovered without permanent protection
        px(i) = px(i) +  SEIRS(cS(index)) +  SEIRS(cS2(index)) + SEIRS(cE1(index)) + SEIRS(cE1Q(index)) + SEIRS(cE1TQ(index)) + SEIRS(cE1TQ2(index))...
            + SEIRS(cE2(index)) + SEIRS(cE2Q(index)) + SEIRS(cE2TQ(index))+ SEIRS(cE2TQ2(index)) +SEIRS(cI(index)) +...
            SEIRS(cA(index)) + SEIRS(cR(index)) + SEIRS(cT(index)) + SEIRS(cP(index)) + SEIRS(cQ(index))...
            + SEIRS(cRA(index)) + SEIRS(cTA(index)) + SEIRS(cTQ(index)) + SEIRS(cT2(index))+ SEIRS(cTA2(index)) + SEIRS(cTQ2(index));
        qstates(i) = qstates(i) + SEIRS(cQ(index)) + SEIRS(cTQ(index)) + SEIRS(cTQ2(index));
    end
end

%calculate current number of people in each strata
for loop = 1:m
    px(loop) = px(loop)/cd.N;
end

% check population is being conserved
tol = 1.e-3;
if (m==2)
    if (mo.symbolic==0)
        if (sum(px) > 1+tol || sum(px) < 1-tol)
            fprintf(1,'time = %i px(1) = %e px(2) = %e\n',time, px(1),px(2));
            params
            sum(px)
            error('ERROR: population in each strata not adding to known total');
        end
    elseif (m==1)
        %fprintf(1,'time = %e px(1) = %e  itot(1) = %e \n',time,px(1),itot(1));
        if (sum(px) > 1+tol || sum(px) < 1-tol)
            fprintf(1,'time = %i px(1) = %e \n',time, px(1));
            params
            sum(px)
            error('ERROR: population in each strata not adding to known total');
        end
    end
end

    
% add noise to itot in each strata proportionally to number of hosts,
% driven by total UK incidence

if (strcmp(mo.seeding,'RSnoise') || strcmp(mo.seeding,'noise'))
    if mo.symbolic ==1
       noise = sym('noise_t_');
    else
       noise = global_noise_lookup(time,mo,params);
    end
    
    if (isnan(noise))
        error('error calculating noise term %e\n',noise);
    end
    for i = 1:m
        itot(i) = itot(i) + heaviside(noise)*noise*px(i);
    end
end

% calculate current waning rates
time_since_seeding = time-params.Tseed1;
[phi phiQ dum1 dum2 dum3 phi_drift phi_koelle] = waning_rate(params,SEIRS,mo,cd,time_since_seeding);

%calculate matrix of relative infectiousness, susceptibility and
%symptomatic proportion (for two strata: experienced and naiive).
epsilon_array = fepsilon(params,mo);

%calculate beta matrix using beta_factor and epsilon matrix
beta = zeros(m);
if mo.symbolic == 1
  beta = sym(beta);
end
for i=1:m
    for j=1:m
        %fprintf(1,'%i %i beta_factor = %e epsilon_array = %e\n',i ,j, beta_factor,epsilon_array(i,j));
        if strcmp(mo.seasonal_forcing,'sinusoidal')
            beta(i,j) =  beta_factor*epsilon_array(i,j)/(1+params.b1*cos(2*pi/52/oo.days_factor*(time+params.SFphase*oo.days_factor)));
        elseif strcmp(mo.seasonal_forcing,'shaman')
            %need to write this in terms of average R0
            beta(i,j) =  beta_factor*epsilon_array(i,j)*(1 + params.SF_ah_factor/R0*exp(-180.*params.b1)*(exp(params.b1*180*(1.-sin(2*pi/52/oo.days_factor*(time+params.SFphase*oo.days_factor))))-1.));
        else
            beta(i,j) =  beta_factor*epsilon_array(i,j);
        end
        if mo.symbolic == 0
            if (beta(i,j) < 0.)
                fprintf(1,'b1 = %e\n',params.b1);
                error('ERROR: beta is negative\n');
            end
        end
    end 
end


% Scale beta by social distancing model
switch (mo.attenuation)
    case {'none'}
        % Nothing to do
        
    case {'SSD','SSDinc'}
        Tstart = max(1,round((time-params.Tmem-params.delta3)/oo.days_factor));
        Tend = min(48,max(1,round((time-params.delta3)/oo.days_factor)));
        M = sum(SSD_time_series(Tstart:Tend));
        beta = beta*params.kappa/(params.kappa+M);
        
    case {'exp'} % per wave attenuation
        if (time >= params.delta*oo.days_factor) && (time < params.Tseed2)
            beta = beta*exp(-k*(time-params.delta*oo.days_factor));
        elseif (time >= parmas.Tseed2+params.delta*oo.days_factor) && (time < params.Tseed3)
            beta = beta*exp(-k*(time-params.Tseed2-params.delta*oo.days_factor));
        elseif (time >= params.Tseed3+params.delta*oo.days_factor)
            beta = beta*exp(-k*(time-params.Tseed3-params.delta*oo.days_factor));
        end
        
    case {'harshSSD'} % per wave attenuation
        if (time >= params.Tseed2-params.delta2) && (time < params.Tseed2)
            beta = beta*params.harshSSDfactor;
        elseif (time >= parmas.Tseed3-params.delta2) && (time < params.Tseed3)
            beta = beta*params.harshSSDfactor;
        end
end

% create pre-equations.  these equations are evaluated before each entry into state vector is calculated.
requiredEqns = struct('varName', {}, 'varData', {} );
varCount = 0;
if mo.symbolic == 1
	varCount = varCount + 1;
	requiredEqns(varCount).varName = 'SEIRSREALEQUATION_phi';
	requiredEqns(varCount).varData = phi;
	phi = sym(requiredEqns(varCount).varName);

	varCount = varCount + 1;
	requiredEqns(varCount).varName = 'SEIRSREALEQUATION_phiQ';
	requiredEqns(varCount).varData = phiQ;
	phiQ = sym(requiredEqns(varCount).varName);
	
	for ii = 1:length(itot)
		varCount = varCount + 1;
		requiredEqns(varCount).varName = strcat('SEIRSREALEQUATION_itot_', int2str(ii));
		requiredEqns(varCount).varData = itot(ii);
		itot(ii) = sym(requiredEqns(varCount).varName);
	end
	
	% now store these in the requiredEqns struct for future processing
end

params1 = struct('beta',beta,'gamma',gamma,'alpha_array',alpha_array,'nu',nu_array,'phi',phi,'phiQ',phiQ,'phix',1/params.Tx,'phi_drift',phi_drift,'phi_koelle',phi_koelle,'rho',params.rho, 'chiTw',params.chiTw,'epsilonS',params.epsilonS,'epsilonalphaQ', params.epsilonalphaQ, 'epsilonSQ', params.epsilonSQ);


% Set initial additional flow into suceptible state to 0
% This could be non-zero if there was a net birth/death rate?

Sin0 = zeros(m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate flow rates between states during first wave %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[SEIRS_dot(cS(1):cD(m)),Sout1a,SAout1a] = RHS_flow_stratified(time, 0, m,params1, params, mo, SEIRS(cS(1):cD(m)), itot, px, Sin0, zeros(m,1), 'cycleA','cycleI','NOseparateS');

SEIRS_dot = SEIRS_dot';


end
