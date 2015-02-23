function [SEIRS_data model_time seeding_flag] = SEIRS_model_new(params,SEIRS_init,model_time,mo,oo,cd,df)
% SEIRS_model - Perform the numerical integration over time
%
% Arguments: params <1xnum_params double>: Model parameters
%            SEIRS_init <1xnum_states double>: Initial state-vector
%            model_time <1xnum_days int>: times (unit of weeks) for ode45
%                                         to return the state-vector
%            df <struct>:  See MCMC_*run.m for details
%            mo <struct>: See MCMC_*run.m for details
%            oo <struct>: See MCMC_*run.m for details
%            cd <struct>: See MCMC_*run.m for details
%            p <struct>: See MCMC_*run.m for details
%
% Output: model_inc <num_daysx1 double>: Infection incidence (per week)
%         model_prev <num_daysx1 double>: Prevalence of infection (number
%         of people)
%         model_inc_STG_final<7 double>: final cumulative incidence in each loop
%         model_inc_STG<7 double>: time_dependenct incidence in each loop
%         seeding_flag <int * 3>: 0 or 1 depending on whether seeding in
%                                 each of the 3 waves was sucessful
%
% NOTE: The unit of time is weeks and we report back model outputs 7 times
% per week (ie every day) 
% K: unit of time is now days

fprintf(1,'parameters as array\n');
fprintf(1, ['%e '],struct2array(params))
fprintf(1,'\n');

stgmax = mo.nstgs;
m = mo.stratifications;
ns = mo.nstates; % number of states per strata per loop
nmultiple = mo.nmultiple; % number of sub loops per loop (required to count secondary infections etc).
loopsmax = stgmax*m*nmultiple;



cS = [];
for j = 1:stgmax*nmultiple
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
cS2=cTQ2+m;
cE1Q = cS2 + m;
cE2Q = cE1Q+m;
cE1TQ = cE2Q +m;
cE2TQ = cE1TQ + m;
cE1TQ2 = cE2TQ + m;
cE2TQ2 = cE1TQ2 + m;
cCUMA = cE2TQ2 + m;


% set up initial state arrays for reinfection options that have multiple
% loops

SEIRS_init0 = [SEIRS_init zeros(1,ns*(stgmax-1)*mo.stratifications*nmultiple + (nmultiple-1)*mo.stratifications*ns)];


options = odeset('NonNegative',1:length(SEIRS_init0));



%calculate beta factors for each wave
symbolic = mo.symbolic;
mo.symbolic = 0;
epsilon_array = fepsilon(params,mo);
mo.symbolic = symbolic;
bf1 = fbeta_factor(cd, mo, params,SEIRS_init0,epsilon_array,params.alpha1,params.R01);
bf2 = fbeta_factor(cd, mo, params,SEIRS_init0,epsilon_array,params.alpha2,params.R02);
bf3 = fbeta_factor(cd, mo, params,SEIRS_init0,epsilon_array,params.alpha3,params.R03);
bf_array = [bf1, bf2, bf3];


%fprintf(1,'peak times %i %i %i\n',peak_time);

%initial set seeding_flag entries to one, these are set to 0 if seeding is
%not successful
seeding_flag = [1 1 1];
    
% only modelling one epidemic wave
t0 = model_time(1);
t1 = params.Tseed1;
t2 = model_time(end);


mo.wave_number = 0;
seeding_flag = [];

if t1>1
    mtimes = 1:t1;
    [junk, SEIRS_data0] = mo.ode_func(mtimes, SEIRS_init0, params, ...
        mo, oo, cd, df, bf_array, ...
        options);
    SEIRS_data0 = check_data_array(SEIRS_data0,mtimes);
    SEIRS_init1 = SEIRS_data0(end,:);
    

else
    SEIRS_init1 = SEIRS_init0;
    SEIRS_data0 = [];
end

mtimes = t1:t2;
[SEIRS_init1] = introduce_seed(SEIRS_init1,mo.wave_number,mo,params);
[junk, SEIRS_data1] = mo.ode_func(mtimes, SEIRS_init1, params, ...
    mo, oo, cd, df, bf_array, ...
    options);
SEIRS_data1 = check_data_array(SEIRS_data1,mtimes);


SEIRS_data = oo.g*[SEIRS_data0(1:end-1,:) ; SEIRS_data1(1:end,:)];
    

fprintf(1,'length SEIRS_data is %i\n',length(SEIRS_data(:,1)));

if isnan(sum(SEIRS_data))
    error('problem with SEIRS_data; at least one element is NaN');
end

end




