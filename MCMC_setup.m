function [MCMC_details iparams N UKmort ndays data_set itt_cum_tot] = MCMC_setup(extra_info,data_set)

global days_factor;
global filename_TTD;
global stratifications;
global extra_days;

% extra parameters governing calculation of likelihood 
force_zero = 15;            % if > 0 add in zero mort (weekly) data points at tail of mort data file
wave4_reinfection_data = 0; % zero to use only given reinfection data, 1 to include point Nijk4 = 0
fit_type = 'mortalityTTD'; % mortalityTTD OR incidence OR waves OR reinfection
r = 10000;                 % Negative binomial r

% Load data corresponding to chosen data set
%[N, N_survey,mort_data,inc_data,wave_data,reinfection_data, ...
%    wave_ranges,mu,UKmort] = LoadData(data_set,force_zero,wave4_reinfection_data);
UKmort = 0;
N=1000000;

% calculate number of days to simulate, to mass to model_setup
ndays = 0;

% specify parameters governing ParallelTempering options
MixingMode ='ParallelTempering';      % Local OR MultipleStarts OR ParallelTempering (PT)
if strcmp(MixingMode,'ParallelTempering')
    n_replicas = 8;            % (PT only) Number of replicas
else
    n_replicas = 1;
end
T_max = 10000;             % (PT only) Maximum temperature
itt =     100;               % (PT) Iterations per swap test
proposal_type = 'group_randomly'; %'group_randomly' (quickest, see MCMCproposal_params.m), 'only_mus_grouped' 
%(change all mu's at same time and all other parameters separately),'single_site' (for old parameter by parameter updates)
no_param_groups = 6;


% Specify number of iterations in MCMC walk (per temperature chain)
itt_tot = 500000;            % Total number of iterations to be performed in this run
itt_cum_tot =  itt_tot;   % also keep track of cumulative number of iterations

% performed with this starting point for use in MCMC_followon.m
Start = 'Provided';        % Random OR Provided or Provided_all_replicas or 'testparams' to test specific parameter set as array in file 'testparams'
testfilename = 'probparams.txt';
provided_start_fname = 'fail_params.txt'; % required for option Start = 'Provided_all_replicas'
%specify step width for each temperature chain (this is also scaled 
flutters = linspace(0.01,0.5,n_replicas); % Relative (to valid range) proposal
% widths for use with multi-chain Parallel tempering
flutters_woPT = 0.0025; %scalar value for MCMC chain without Parallel Tempering
%flutters_woPT = 0.025; %scalar value for MCMC chain without Parallel Tempering
flutters_HP = 0.006; %flutter for hyperparameter proposal function

% Some constants that should not require changing
I_0 = 1;                   % Only required if not fitting seeds
g = 1;                     % No untimed asymptomatics (see RAF/TdC paper)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Fix initial parameter values, and choose which parameters to %%%
%%%%%%%%%% vary, and over what range %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% no reseeding, single epidemic wave only
Tseed2_min =  0;
Tseed3_min = 0;


m = stratifications;
initial_pop_dist = 0.8;
pBcellE = 0.0; % proportion of experienced stratum with B-cell protection to pandemic strain

 
  parameters = struct( ...
    'xI', struct('name','xI',                        'fitted', zeros(1,m), 'range',[0.0 1.0], 'ival', [pBcellE 0 pBcellE 0]),...                            % proportion of stratum 1 (experienced) who are permanently immune
    'R01', struct ('name','R01',                     'fitted',1,'range',[0.5,10.0],'ival',2), ...                                        % Basic reproduction number in wave 1
    'R02', struct ('name','R02',                     'fitted',0,'range',[1.,10.0],'ival',9), ...                                          % Basic reproduction number in wave 2, If 0, R02 = R01
    'R03', struct ('name','R03',                     'fitted',0,'range',[1.,10.0],'ival',9), ...                                          % Basic reproduction number in wave 3, If 0, R03 = R01/R02
    'alpha1', struct ('name','alpha1',                'fitted',1,'range',[0,1.0],'ival',1.0), ...                                      % symptomatic proportion in wave 1 (must be > 0)
    'alpha2', struct ('name','alpha2',                'fitted',0,'range',[0,1.0],'ival',1.), ...                                         % symptomatic proportion in wave 2, If 0, alpha2 = alpha1 (must be > 0) 
    'alpha3', struct ('name','alpha3',                'fitted',0,'range',[0,1.0],'ival',1.), ...                                        % symptomatic proportion in wave 3, If 0, alpha3 = alpha1(must be > 0)
    'epsilonALPHA', struct ('name','epsilonALPHA',    'fitted',1,'range',[0,1.0],'ival',1.0), ...                                          % scale alphas in stratum 1 by epsilon_alpha'epsilonALPHA'
    'epsilonS', struct ('name','epsilonS',            'fitted',1,'range',[0,1.0],'ival',1.0), ...                                          % scale suceptibility in stratum 1 by epsilonS
    'epsilonI', struct ('name','epsilonI',            'fitted',1,'range',[0,1.0],'ival',1.0), ...                                          % scale infectiousness in stratum 1 by epsilonI
    'epsilonNU', struct ('name','epsilonNU',          'fitted',0,'range',[0,1.0],'ival',1.0), ...                                          % scale infectious period in stratum 1 by this factor'fSE', struct('name','fSE',                       'fitted',0,'range',[0,1.0],'ival',1.0), ...                                         % scale susceptibility of experienced stratum with T cell vaccine
    'fSE', struct('name','fSE',                       'fitted',0,'range',[0,1.0],'ival',1.0), ...                                         % scale susceptibility of experienced stratum with T cell vaccine
    'fNUE', struct('name','fNUE',                     'fitted',0,'range',[0,1.0],'ival',1.0), ...                                         % scale infectious period of experienced stratum with T cell vaccine
    'fIE', struct('name','fIE',                       'fitted',0,'range',[0,1.0],'ival',0.25), ...                                         % scale infectiousness of experienced stratum with T cell vaccine
    'fALPHAE', struct('name','fALPHAE',                'fitted',0,'range',[0,1.0],'ival',1.0), ...                                         % scale symptomatic proportion of experienced stratum with T cell vaccine (not currently used)
    'fSN', struct('name','fSN',                       'fitted',0,'range',[0,1.0],'ival',1.0), ...                                         % scale susceptibility of naiive stratum with T cell vaccine
    'fNUN', struct('name','fNUN',                     'fitted',0,'range',[0,1.0],'ival',1.0), ...                                         % scale infectious period of naiive stratum with T cell vaccine
    'fIN', struct('name','fIN',                       'fitted',0,'range',[0,1.0],'ival',0.25), ...                                         % scale infectiousness of naiive stratum with T cell vaccine
    'fALPHAN', struct('name','fALPHAN',                'fitted',0,'range',[0,1.0],'ival',1.0), ...                                         % scale symptomatic proportion of naiive stratum with T cell vaccine (not currently used)
    'epsilonSQ', struct('name','epsilonSQ',           'fitted',1,'range',[0,1],'ival',0.),...                                          % scale susceptility in prior protection states by this factor ^ (number of states to S-1)
    'epsilonalphaQ', struct('name','epsilonalphaQ',   'fitted',1,'range',[0,1],'ival',1.0),...                                            % scale symptomatic proportion in prior protection states by this factor ^ (number of states to S-1) WARNING: not currently correctly implemented.i.e. beta not scaled up accordingly
    'rho', struct ('name','rho',                      'fitted',zeros(1,m),'range',[0,1.0],'ival',[0  0 0 0]), ...                                  % proportion in stratum 1 who gain permenent protection following recovery
    'z', struct ('name','z',                          'fitted',zeros(1,m),'range',[0,1.0],'ival',[1 1 1 1]), ...                                 % proportion of stratum 1 who are initially susceptible
    'Te', struct ('name','Te',                        'fitted',0,'range',[0.5 2.0],'ival',1.6), ...                                       % latent period (days)
    'Ti', struct ('name','Ti',                        'fitted',1,'range',[0.0 3.0],'ival',1.3), ...                                      % infectious period (days)                                 
    'b1', struct ('name','b1',                        'fitted',1,'range',[0.0 0.5],'ival',0.0), ...                                    % amplitude of seasonal forcing
    'SFphase', struct('name','SFphase',               'fitted',1,'range',[14,26],'ival',19),...                                           % phase of seasonal forcing
    'SF_ah_factor',struct('name','SF_ah_factor',      'fitted',0,'range',[0,10],'ival',1),...                                             % additional sf parameter for Shaman version (not fully implemented yet)
    'Tseed1', struct ('name','Tseed1',                'fitted',1,'range',[1 100],'ival',20), ...                                 % time (days) of introduction of seed to start wave 1
    'I0',struct('name','I0',                          'fitted',0,'range',[0,1],'ival',0.00002),...    % detection threshold
    'lgTw', struct ('name','lgTw',                    'fitted',1,'range',[0,100],'ival', 100), ...      %Koelle SM (from Dushoff et al. 2004) % natural log of period of acquired temporary infection in days
    'chiTw', struct('name','chiTw',                   'fitted',1,'range',[1,100],'ival',1),...                                            % increase in period of acquired temporary infection for secondary/tertiary infections
    'lgTwQ', struct ('name','lgTwQ',                  'fitted',1,'range',[0.0 100],'ival',100), ...                                     % natural log of period of prior temporary infection in days  
    'kwan', struct ('name','kwan',                    'fitted',1,'range',[-100 2.],'ival',-100.), ...                                        % factor scaling state dependent acceleration of waning of acquired protection due to time dependent drift rate
    'kwanQ', struct ('name','kwanQ',                  'fitted',1,'range',[-3.7 2.],'ival',-3.5), ...                                     % factor scaling state dependent acceleration of waning of prior protection due to time dependent drift rate
    'Tx',struct('name','Tx',                          'fitted',0,'range',[30,187],'ival',60),...                                          % period of cross-reactive immunity (days)
    'mu1', struct ('name','mu1',                      'fitted',1,'range',[0,1],'ival',.1), ...                                 % death rate per infection wave 1
    'mu2', struct ('name','mu2',                      'fitted',1,'range',[0,1],'ival',.1), ...                                 % death rate per infection wave 2
    'mu3', struct ('name','mu3',                      'fitted',1,'range',[0,1],'ival',.1), ...                                 % death rate per infection wave 3
    'epsilonMU', struct ('name','epsilonMU',           'fitted',0,'range',[0 1],'ival',1.), ...                                           % scale death rate in stratum 1 by this factor
    'pop_dist', struct('name','pop_dist',             'fitted', ones(1,max(1,m-1)), 'range',[0.0 1.0], 'ival', [initial_pop_dist 0 0]),...      % proportion of population beginning exp, prop of exp beginning with T cell vaccine, prop of naiive beginning with T cell vaccine
    'log_noise_scale',struct('name','log_noise_scale', 'fitted', 0, 'range',[-20,-1.],'ival',-20.),... % factor scaling 'noise' term due to external epidemic
    'kw_koelle', struct('name','kw_koelle',            'fitted', 0, 'range',[1, 5],'ival',1.5),...  % shape parameter for Weillbull function describing probability of drift variant / new cluster
    'lambda_koelle', struct('name','kw_koelle',            'fitted', 0, 'range',[500, 1000],'ival',500),...  % 
    'threshARkoelle',struct('name','threshARkoelle',    'fitted',0,'range',[0,1],'ival',1/N),...
    'drift_norm',struct('name','drift_norm',           'fitted', 0, 'range',[1.,1e7],'ival',1.),... % total drift over epidemic for unimpeded epidemic
    'drift_koelle_norm',struct('name','drift_koelle_norm','fitted',0,'range',[1.,1e7],'ival',1.),...   % total drift over unimpeded epidemic for koelle function
    'N',struct('name','N',                             'fitted',0,'range',[1,1e8],'ival',1.e6),... % total population size
    'fracTcellvac',struct('name','fracTcellvac',        'fitted',0,'range',[0,1],'ival',0.0),... % fraction of poulation receiving T cell vaccin
    'periodTcellvac',struct('name','periodTcellvac',    'fitted',0,'range',[0,365],'ival',10),... % period of which above vaccine is distributed (days)
    'distTcellvac',struct('name','distTcellvac',        'fitted',0,'range',[0 0.5],'ival',1.),... % bias in vaccination likelihood, used in 'biasvacE' and 'biasvacN' vaccination strategies, failes if too small...
    'T_Tcellvac',struct('name','T_Tcellvac',            'fitted',0,'range',[0,365],'ival',1),... % day of simulation at which CTL vaccine rollout begins (not implemented yet)
    'fracBcellvac',struct('name','fracBcellvac',        'fitted',0,'range',[0,1],'ival',0.0),... % fraction of poulation given strain-specific B cell vaccine
    'periodBcellvac',struct('name','periodBcellvac',    'fitted',0,'range',[0,365],'ival',30),... % period over which above vaccine is distributed
    'T_Bcellvac',struct('name','T_Bcellvac',            'fitted',0,'range',[-1,365],'ival',1),...  % day of simulation at which above distribution begins (at least 3 months after pandemic identified?)                                    
    'phi_Tcell',struct('name','phi_Tcell',              'fitted',0,'range',[0,0.1],'ival',0),... % inverse period of protection from T cell vaccine (days)
    'phi_Bcell',struct('name','phi_Bcell',              'fitted',0,'range',[10,1000],'ival',1./158),...% inverse period of protection from B cell vaccine (days)   NOT YET IMPLEMENTED
    'beta11mix',struct('name','beta11mix',              'fitted',0,'range',[0,5],'ival',0.7),...% relative force of infection between hosts in stratum 1 due to mixing
    'beta22mix',struct('name','beta22mix',              'fitted',0,'range',[0,5],'ival',1.4),...% relative force of infection between hosts in stratum 2 due to mixing
    'beta21mix',struct('name','beta21mix',              'fitted',0,'range',[0,5],'ival',1.0),...% relative force of infection between hosts in stratum 2 on 1 due to mixing
    'beta12mix',struct('name','beta12mix',              'fitted',0,'range',[0,5],'ival',1.0)) % relative force of infection between hosts in stratum 2 on 1 due to mixing
    
% Sanity check
if xor(parameters.alpha2.fitted,parameters.alpha3.fitted)
    error('Invalid setting: Either both or neither of alpha2 and alpha3 should be fitted\n')
end

% Sanity check
if xor(parameters.alpha2.fitted,parameters.alpha3.fitted)
    error('Invalid setting: Either both or neither of alpha2 and alpha3 should be fitted\n')
end

% construct information about fitted and fixed parameters

fit_params = [];
const_params = [];
pnames = fieldnames(parameters);
param_num = 0;
for i=1:length(pnames)
    [xel yel] = size(parameters.(pnames{i}).ival);
    for x = 1:xel
        for y = 1:yel
            param_num = param_num+1;
            pnames{i}
            if (parameters.(pnames{i}).fitted(x,y))
                fit_params = [fit_params param_num];
                %check ival for parameter is in specified range
                if ((parameters.(pnames{i}).ival(x,y) < parameters.(pnames{i}).range(1)) ...
                        || (parameters.(pnames{i}).ival(x,y) > parameters.(pnames{i}).range(2)))
                    fprintf(1,'%s\n',parameters.(pnames{i}).name);
                    error('ERROR: initial parameter value not in specified range\n');
                end
            else
                const_params = [const_params param_num];
            end
        end
    end
end



% pool together parameters into structure MCMC_details

city_details = struct('N', N, ...
    'N_survey', 0, ...
    'mort_data', 0, ...
    'inc_data', 0, ...
    'wave_data', 0, ...
    'reinfection_data', 0, ...
    'wave_ranges', 0, ...
    'mu', 0);


data_fit_options = struct('fit_type', fit_type, ...
    'TTD', 0, ...
    'r', r, ...
    'MixingMode', MixingMode, ...
    'n_replicas', n_replicas, ...
    'T_max', T_max, ...
    'itt', itt, ...
    'itt_tot', itt_tot, ...
    'itt_cum_tot', itt_cum_tot, ...
    'Start', Start, ...
    'flutters', flutters, ...
    'flutters_woPT',flutters_woPT,...
    'flutters_HP',flutters_HP,...
    'force_zero',force_zero,...
    'wave4_reinfection_data',wave4_reinfection_data,...
     'proposal_type',proposal_type,...
    'no_param_groups',no_param_groups);

parameter_details = struct('names', {pnames}, ...
    'num_params', length(pnames), ...
    'num_fitted_params',length(fit_params), ...
    'fitted', fit_params, ...
    'constant', const_params);



%generate cell array of initial parameters for each chain, if provided
iparams = cell(n_replicas,1);

iparams{1} = Genparams(parameter_details,parameters);
if (n_replicas>2)
    for i =2:n_replicas
        iparams{i} = [];
    end
end
    
if strcmp(data_fit_options.Start,'Random')
    for i = 1:pd.num_fitted_params
         p_string = parameter_details.names{pd.fitted(i)};
         iparams.p_string = parameters.(p_string).range(1) + ...
             rand*(parameters.(p_string).range(2) - parameters.(p_string).range(1));
    end
elseif strcmp(data_fit_options.Start,'Provided_all_replicas')
    % open provided_start_fname
    big_matrix = importdata(provided_start_fname);
    
    % read in parameter arrays and convert to parameter structure and feed to iparams
    for i = 1:min(data_fit_options.n_replicas,size(big_matrix,1))
        params_array = big_matrix(i,:);
        iparams{i} = param_struct_from_array_v2(params_array,parameter_details,parameters);
    end
elseif strcmp(data_fit_options.Start,'Provided')
    % use initial values as input into parameter structure above
elseif strcmp(data_fit_options.Start,'testparams')
    params_array = importdata(testfilename);
    iparams{1} = param_struct_from_array_v2(params_array,parameter_details,parameters);
    %overide initial parameter values for parameters that are fixed and not
    %fitted
    for i=1:length(pnames)
        [xel yel] = size(parameters.(pnames{i}).ival);
        for x = 1:xel
            for y = 1:yel
                param_num = param_num+1;
                
                if (parameters.(pnames{i}).fitted(x,y)==0)
                    iparams{1}.(pnames{i})(x,y) = parameters.(pnames{i}).ival(x,y)
                
                end
            end
        end
    end
    
    % make now make sure that R0 and alphas are the same betweeen waves in
    % values in waves 2 and 3 aren't fitted
    
    if (parameters.alpha2.fitted == 0)
        iparams{1}.alpha2 = iparams{1}.alpha1;
    end
    if (parameters.alpha3.fitted == 0)
        iparams{1}.alpha3 = iparams{1}.alpha1;
    end
    
    if ((parameters.R02.fitted == 0) &&(parameters.R03.fitted == 0))
        iparams{1}.R02 = iparams{1}.R01;
        iparams{1}.R03 = iparams{1}.R01;
    end
    
    if (parameters.R03.fitted == 0)
        iparams{1}.R03 = iparams{1}.R01;
    end
    
    
else    
    error('ERROR: Invalid "Start" (must be either "Random" or "Provided" or "Provided_all_replicas"')
end

% Establish filenames for job
filename_base = MCMCfileinit(data_set,extra_info,itt_cum_tot);
filename_for_workspace = strcat('Workspace-',filename_base);
filename_for_PTswaps = strcat('PTswaps-',filename_base);

filename_for_params = cell(n_replicas,1);
for i = 1:n_replicas
    filename_for_params{i} = strcat('Params-',filename_base,'_T_',int2str(i));
end

other_options = struct('I_0', I_0, ...
    'g', g, ...
    'days_factor', days_factor);

% no fitting performed
prior_information = 0;

MCMC_details = struct('df', data_fit_options, ...
    'cd',city_details,...
    'pd', parameter_details, ...
    'p', parameters, ...
    'pi',prior_information,...
    'fname', {filename_for_params}, ...
    'filename_for_workspace',filename_for_workspace,...
    'filename_for_PTswaps',filename_for_PTswaps,...
    'param_last',[],...
    'data_set',data_set,...
    'oo',other_options);


    
