function [model_details] = model_setup(N,UKmort,ndays)

global days_factor;
global stratifications;
global filename_TTD;
global extra_days;

% choose 1 to generate c states corresponding to each wave for construction
% of odecpp for this particular model, or 0 to run normally with the
% ode solver specified by ode_func. Choices for ode_func are
% MATLAB_ode_wrapper (for ode45) or ode_cpp for the optimised code (latter
% needs to be generated manually first).


symbolic = 0;
ode_func = @MATLAB_ode_wrapper;
%ode_func = @ode_cpp;


% OpenMP settings
% use the following to set the number of threads you would like to use. 
% this option is only relevant for the ode_cpp solver.
number_openmp_threads = 1;
setenv('OMP_NUM_THREADS', num2str(number_openmp_threads)); 


%%%%%%%%%%%%%%%%%%%%%%%
%%% Model options %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

% Cyling, seeding, noise, dieout, waning and social distancing options
protected_population_state = 'Q'; %'TQ' to have 2 state waning to get to S, 'Q' to have 3 state waning to get to S 
globaldieout ='GD';        % 'noGD' OR 'GD' (Don't use noGD)
seeding = 'RS';            % 'RSifGD' OR 'RS' OR 'noRS' OR 'noise' OR 'RSnoise
RSoption = [1 1 1];        % entry 1, 2, 3 set to 1 if want to reseed for that wave, 0 otherwise
attenuation ='none';       % 'none' OR 'exp' OR 'SSD' OR 'SSDinc' OR 'harshSSD'
waning = 'mutwan';         % 'normal' OR 'reset' or 'mutwan'
                           % Further info: normal: Tw, reset: S = z
                           % mutwan: Waning proportional to summed mutation
                           % rate which is proportional to incidence, either
                           % the city's incidence or the whole of
                           % the UK's incidence, determined by inctype_set
mutwan_form = 'neuconst';   % 'incprev' (all terms proportional to current prevalence) OR 'neuconst' (neutral drift is independent of prevalence)
mutwan_type = 'icum';  % 'icum' to have drift term scale with number of infections in current wave
inctype_set = 'city';     
seasonal_forcing = 'sinusoidal'; % choose 'sinusoidal' or 'shaman'
% If 'waning=mutwan', use 'city' OR 'UK' incidence
reinfection = 'noreinfection'; %'reinfection' to use reinfection data and 7stage model, 
% or 'snake' to use 3 consecutive loops and not count reinfections,
% 'noreinfection' for cycling in same loop
divert_between_strata = 'divert'; %'divert': to divert recovered I's and A's to different strata, 
                                  %'nodivert': conserve population in each strata
%transient_xprotection = ''; %choose 'transient_xprotection' in order to have first stage of waning (i.e. R->T) happen on time-scale 1/Tx (Tx in days)
transient_xprotection = 'transient_xprotection';
koelle_drift = 'punc'; %'cont' continuous to use same form as Boni drift, punc (i.e. punctuated) to assume that antigenic type can change at end of season or not
drift_rate = 'percapita'; %percapita=every infected hosts has same chance of generating mutant, propinfectiousness=prob of generating mutant proportional to total viral load of infection
% CTLdependent=those with elevated CTLs (i.e. vaccinated adults) do not
% generate mutant due to rapid clearance of infection by CTLs

% parameters specifying loop and model indices, these should not need to
% change for the reinfection model
nstgs = 1; % max number of states
nmultiple = 1; % number of loops per stage
nstates = 25; % number of states per strata per loop


priorBcells = 'permanent'; % 'temporary' to start fraction of hosts defined by xI in R, 'permanent' to start them in P

transmission = 'fixReff'; % choose 'fixReff' to chose effective reproduction number (given prior immunity specified via parameters) or 'fixbeta' to choose transmission rate independently of existing immunity (using R0 and assuming naiive population)
vaccination = 'Tcell'; % choose 'none', 'Tcell', 'TBcell', 'Bcell'
vacc_strat = 'unvac'; % also choose 'unvac', in future target E 'vacE' or N 'vacN' stratum, or bias vaccination toward one stratum or the other with 'biasvacE','biasvacN'
vacc_timing = 'dynamic'; % choose 'prepandemic' to begin simulation with hosts vaccinated or 'dynamic' to distribute as epidemic progresses
vacc_action = 'boost'; % boost for Tcell efficacy to compound existing cross-protection 'sat' to saturate to vaccine efficacy
Tcellvac_exposed = 0; % 1: Tcell vaccination changes susceptibility of infected or exposed hosts, 0: does not

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose model functions appropriate for specified model options %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model_rhs = @RHS_7stage;
model_SEIRS_func = @SEIRS_model_new;


UK_inc_linear = []

model_options = struct('globaldieout', globaldieout, ...
    'seeding', seeding, ...
    'RSoption',RSoption, ...
    'attenuation', attenuation, ...
    'waning', waning, ...
    'mutwan_form',mutwan_form,...
    'mutwan_type',mutwan_type,...
    'global_noise',UK_inc_linear,...
    'inctype_set', inctype_set, ...
    'model_rhs', model_rhs, ...
    'model_SEIRS_func',model_SEIRS_func, ...
    'seasonal_forcing',seasonal_forcing,...
    'reinfection',reinfection,...
    'stratifications',stratifications,...
    'nstates',nstates,...
    'nstgs',nstgs,...
    'nmultiple',nmultiple,...
    'divert_between_strata',divert_between_strata,...
    'protected_population_state',protected_population_state,...
    'symbolic',symbolic,...
    'ode_func',ode_func,...
    'wave_number',0,...
    'cycleRAs','',...
    'transient_xprotection',transient_xprotection,...
    'ndays',ndays,...
    'transmission',transmission,...
    'vaccination',vaccination,...
    'vacc_strat',vacc_strat,...
    'vacc_action',vacc_action,...
    'vacc_timing',vacc_timing,...
    'Tcellvac_exposed',Tcellvac_exposed,...
    'priorBcells',priorBcells,...
    'koelle_drift',koelle_drift,...
    'drift_rate',drift_rate);



model_details = struct('mo', model_options);
