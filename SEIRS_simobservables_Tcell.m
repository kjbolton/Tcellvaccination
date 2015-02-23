function [model_inc_total total_infections drift drift_koelle Reffective pvac gthreshtime px infections_stratum ratioReff tend] = SEIRS_simobservables_Tcell(params,SEIRS_data,model_time,cd,df,mo,oo)
% SEIRS_simobservables - Given the parameters, initial conditions
%                        state-vector and model specification call the
%                        appropriate model and return the simulation output
%
% Arguments: params <1xnum_params double>: Model parameters
%            SEIRS_init <1xnum_states double>: Initial state-vector
%            model_time <1xnum_days int>: times (unit: days) for ode45 to
%                                         return the state-vector
%            df <struct>: See MCMC_*run.m for details
%            mo <struct>: See MCMC_*run.m for details
%            oo <struct>: See MCMC_*run.m for details
%            cd <struct>: See MCMC_*run.m for details
%
% Output: modelM <num_daysx1 double>: Mortality incidence (per day)
%         model_inc <num_daysx1 double>: Infection incidence (per day)
%         modelW <7x1 double>: Wave (and lull) totals
%
% NOTE: The unit of time is days and we report back model outputs each day


global gthreshtime



fprintf(1,'in SEIRS_simobs\n');

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
cD = cCUMA + m;

length_sim = length(SEIRS_data(:,1));

model_prev_LOOP = zeros(length_sim,loopsmax);
model_inc_LOOP = zeros(length_sim,loopsmax);
model_cum_LOOP = zeros(1,loopsmax);

model_prev_STG = zeros(length_sim,stgmax);
model_inc_STG = zeros(length_sim,stgmax);

model_inc_total = zeros(length_sim,1);

total_infections = 0;
infections_stratum = zeros(m,1);

Reffective = zeros(length_sim,1);

for strata = 1:m
    for stg = 1:stgmax
        for k = 1:nmultiple
            loop = strata + m*(k-1) + (stg-1)*m*nmultiple;
            prev_cumulative = [0 ; SEIRS_data(1:end-1,cCUM(loop)) + SEIRS_data(1:end-1,cCUMA(loop))];
            %first value of this meaningless          
            model_prev_LOOP(:,loop) = SEIRS_data(:,cI(loop)) + SEIRS_data(:,cA(loop));
            model_inc_LOOP(:,loop) = SEIRS_data(:,cCUM(loop)) + SEIRS_data(:,cCUMA(loop)) - prev_cumulative;
            model_inc_total(:) = model_inc_total(:)+model_inc_LOOP(:,loop);
            model_cum_LOOP(1,loop) = SEIRS_data(end,cCUM(loop)) + SEIRS_data(end,cCUMA(loop));
        end
    end
end

% estimate time at which epidemic dies out (total incidence less than 1)
% from saved time 
% find time of maximum incidence
[incmax tmax] = max(model_inc_total);

tend = 0;
for kk=1:length_sim
    if tend==0 && kk>=tmax && (model_inc_total(kk)<1)
        tend=kk;
    end
end


% attack rate in each stratum

% sum over strata to calculate total in each stage
for stg = 1:stgmax
    for j = 1:m
        % sum over subloops
        for k = 1:nmultiple
            model_prev_STG(:,stg) = model_prev_STG(:,stg) + model_prev_LOOP(:,m*nmultiple*(stg-1) + j + (k-1)*m);
            model_inc_STG(:,stg) = model_inc_STG(:,stg) + model_inc_LOOP(:,m*nmultiple*(stg-1) + j + (k-1)*m);
            infections_stratum(j) = infections_stratum(j) + model_cum_LOOP(1,m*nmultiple*(stg-1) + j + (k-1)*m);
        end
    end
end

total_infections = sum(infections_stratum)/params.N;

% calculate effective reproduction number throughout simulation
epsilon = fepsilon(params,mo);
beta = fbeta_factor(cd, mo, params,SEIRS_data(1,:),epsilon,params.alpha1,params.R01);

for i=1:length_sim
   Reffective(i) = fReff(cd,mo,params,SEIRS_data(i,:),epsilon,beta);
end

px = zeros(length_sim,m);

for i =1:m
    for j = 1:stgmax*nmultiple
        % sum over loops and subloops for each strata
        % total number infected
        index = (j-1)*m+i;
        % modify by relative infectiousness epsilonI for experienced
        % stratum (this is done in matrix epsilon used in
        % RHS_flow_stratified.m)
        
        %total number susceptible or recovered without permanent protection
        px(:,i) = px(:,i) +  SEIRS_data(:,cS(index)) +  SEIRS_data(:,cS2(index)) + SEIRS_data(:,cE1(index)) + SEIRS_data(:,cE1Q(index)) + SEIRS_data(:,cE1TQ(index)) + SEIRS_data(:,cE1TQ2(index))...
            + SEIRS_data(:,cE2(index)) + SEIRS_data(:,cE2Q(index)) + SEIRS_data(:,cE2TQ(index))+ SEIRS_data(:,cE2TQ2(index)) +SEIRS_data(:,cI(index)) +...
            SEIRS_data(:,cA(index)) + SEIRS_data(:,cR(index)) + SEIRS_data(:,cT(index)) + SEIRS_data(:,cP(index)) + SEIRS_data(:,cQ(index))...
            + SEIRS_data(:,cRA(index)) + SEIRS_data(:,cTA(index)) + SEIRS_data(:,cTQ(index)) + SEIRS_data(:,cT2(index))+ SEIRS_data(:,cTA2(index)) + SEIRS_data(:,cTQ2(index));
       
    end
end


drift = SEIRS_data(:,cD(1));

if strcmp(mo.koelle_drift,'cont')
    drift_koelle = SEIRS_data(:,cD(2));
elseif strcmp(mo.koelle_drift,'punc')
    Reffend = Reffective(end);
    
    % 'erase' B cell protection to calc max Reff for antigenically novel
    % strain
    SEIRS_data_new = SEIRS_data(end,:);
    for strata = 1:m
    for stg = 1:stgmax
        for k = 1:nmultiple
            
            loop = strata + m*(k-1) + (stg-1)*m*nmultiple;
            
          
            
            SEIRS_data_new(cS(loop)) = SEIRS_data_new(cS(loop)) + SEIRS_data_new(cR(loop)) + SEIRS_data_new(cRA(loop)) +...
            SEIRS_data_new(cT(loop)) + SEIRS_data_new(cP(loop)) + SEIRS_data_new(cTA(loop)) + SEIRS_data_new(cT2(loop)) + ...
            SEIRS_data_new(cTA2(loop));
            SEIRS_data_new(cR(loop)) = 0;
            SEIRS_data_new(cRA(loop)) = 0;
            SEIRS_data_new(cT(loop))  = 0.;
            SEIRS_data_new(cP(loop)) = 0;
            SEIRS_data_new(cTA(loop)) = 0;
            SEIRS_data_new(cT2(loop)) = 0;
            SEIRS_data_new(cTA2(loop)) = 0;
            
        end
    end
    end

    [phi phiQ stoteff stoteffalpha E1tot phi_drift phi_koelle SRtot Stot] = waning_rate(params,SEIRS_data(end,:),mo,cd,0)
    
    % if T cell immunity from vaccination is maintained
    Reffendmax = fReff(cd,mo,params,SEIRS_data_new,epsilon,beta)
    
    ratioReff = Reffendmax/Reffend;
    
    selection = 1;
    
    drift_koelle = SEIRS_data(:,cD(2))*selection;
end



% fraction in each strata
px = px/params.N;

%calculate prepandemic vaccinated fraction

px(1,:)
fprintf('prepandemic vac %f c.f. specified %f\n',px(1,3)+px(1,4),params.fracTcellvac)

% unvac and vac fractions
pvac = zeros(length_sim,2);
pvac(:,1) = px(:,1) + px(:,2); %unvac
pvac(:,2) = px(:,3) + px(:,4); %vac

