function MCMC_ntestrun(run_type)

close all


% MCMC procedure for LOCAL, MULTIPLE STARTS or PARALLEL TEMPERING

% set global variables
global extra_days
global days_factor
global stratifications;
global filename_TTD;
global gthreshtime;

days_factor = 7;
stratifications = 4;
extra_days = 0;
gthreshtime = -1;


extra_info = 'test';
data_set = 'composite';

%input 'run_type' is 'test' to run model once at initial parameter
%estimates and 'full' to do MCMC experiment

[MCMC_details iparams N UKmort ndays data_set] = MCMC_setup(extra_info,data_set);
% overwrite ndays
%ndays = 182;
ndays = 365;

[model_details] = model_setup(N, UKmort, ndays);

filename='Uboostassortpercapita_Reff2_Oct14';

plotstratum =0;

if strcmp(run_type,'test')


if strcmp(model_details.mo.vacc_timing,'prepandemic')
    % Plot model for initial paraeters chosen in MCMC_setup
            % explictly put people in vaccine class for prepandemic option
            % set pop_dist(2) (frac exp vaccinated at start)
            iparams{1}.pop_dist(2) = iparams{1}.fracTcellvac;
            % set pop_dist(3) (frac nai vaccinated at start)
            iparams{1}.pop_dist(3) = iparams{1}.fracTcellvac;

% target adults
%iparams{1}.pop_dist(2)=iparams{1}.fracTcellvac/iparams{1}.pop_dist(1);
%iparams{1}.pop_dist(3) = 0;

%target children
%iparams{1}.pop_dist(2) = (iparams{1}.fracTcellvac-(1-iparams{1}.pop_dist(1)))/iparams{1}.pop_dist(1);
%iparams{1}.pop_dist(3) = 1;
end

        
    tic
    [SEIRS_data model_time] = single_simulation(iparams{1},MCMC_details.cd, MCMC_details.df,model_details.mo,MCMC_details.oo);
    toc
    [model_inc total_infections drift drift_koelle Reffective pvac threshtime px infections_stratum ratioReff tend] = SEIRS_simobservables_Tcell(iparams{1},SEIRS_data,model_time,MCMC_details.cd,MCMC_details.df,model_details.mo,MCMC_details.oo);
    
    Figure = figure(3);    
    subplot(2,3,1)
    %ylim([0,0.05])
    hold on
    plot(model_inc)
    
    plot(threshtime,0,'b.');
    %plot(1,total_infections,'b.')
    xlabel('Time (days)')
    ylabel('Daily incidence')
    
    subplot(2,3,2)
    %ylim([0,2])
    hold on
    plot(drift_koelle)
    xlabel('Time (days)')
    ylabel('relative prob of displacement by drift variant (Koelle)')
    
    subplot(2,3,3)
    ylim([0,1])
    hold on
    plot(drift)
    xlabel('Time (days)')
    ylabel('relative prob of displacement by drift variant (Boni)')
    
    
    subplot(2,3,4)
    ylim([0,4])
    hold on
    plot(Reffective)
    plot(ones(1,length(Reffective)))
    xlabel('Time (days)')
    ylabel('Reff')
    
    subplot(2,3,5)
    ylim([0 1])
    plot(pvac(:,1),'r');
    hold on
    plot(pvac(:,2),'b');
    xlabel('Time (days)')
    ylabel('un/vaccinated proportions');
    plot(pvac(:,1)+pvac(:,2),'Linestyle','--')
    plot(px(:,3),'g');
    plot(px(:,4),'y');
    ylim([0,1.1])
    
    
    
   
    fprintf(1,'fraction of experienced hosts vaccinated %f\n',px(end,3)/px(1,1))
    fprintf(1,'fraction of naive       hosts vaccinated %f\n',px(end,4)/px(1,2))
    fprintf(1,'specified bias %f\n',iparams{1}.distTcellvac);
    fprintf(1,'ratio of vaccinated fractions %f\n',px(end,3)/px(1,1)/(px(end,4)/px(1,2)));
    fprintf(1,'total fraction vaccinated %f compared to specified coverage %f\n',px(end,3)+px(end,4),iparams{1}.fracTcellvac);
    
    
    fprintf(1,'total number infections / N = %f\n',total_infections);
    fprintf(1,'total number of experienced hosts infected / N = %f\n',(infections_stratum(1)+infections_stratum(3))/iparams{1}.N/px(1,1));
    fprintf(1,'total number of naiive hosts infected / N = %f\n',(infections_stratum(2)+infections_stratum(4))/iparams{1}.N/px(1,2));
drift_koelle

fprintf(1,'drift opportunity koelle = %e\n',drift_koelle(end))

%tend
    
elseif strcmp(run_type,'exploretiming')
    
    % fix vaccine coverage at 50 per cent
    % fix vaccine efficacy at 50 per cent
    
    % cycle over day of vaccine rollout (1 to 30 days)
    % cycle over duration of vaccine rollout (10 to 100 days)
    
    nsamples_beg = 20;
    nsamples_dur = 20;
   
    adrift = zeros(nsamples_beg,nsamples_dur);
    adrift_koelle = zeros(nsamples_beg,nsamples_dur);
    adrift_plot = zeros(nsamples_beg,nsamples_dur);
    adrift_koelle_plot = zeros(nsamples_beg,nsamples_dur);
    aAR = zeros(nsamples_beg,nsamples_dur);
    aARexp = zeros(nsamples_beg,nsamples_dur);
    aARnai = zeros(nsamples_beg,nsamples_dur);
    aReff = zeros(nsamples_beg,nsamples_dur);
    aratioReff = zeros(nsamples_beg,nsamples_dur);
    
    
    model_details.mo.vacc_timing = 'dynamic';
    model_details.mo.vacc_strat = 'none';
    
    
    % vary timing of first vaccination
    abegin = linspace(1,30,nsamples_beg);
    % vary duration of vaccination campaign
    aduration = linspace(10,100,nsamples_dur);
    
    for i=nsamples_beg
        for j = 1:nsamples_dur
        iparams{1}.T_Bcellvac = abegin(i);
        iparams{1}.periodTcellvac = aduration(j);
        
        
        
        
        end
    end
    
elseif strcmp(run_type,'explorecoverage')
    
    % vary vaccine coverage for fixed prior protection and epidemic
    % characteristics
    
    % choose prepandemc immunity in MCMC_ntestrun.m
    
    
    % choose dynamic vaccination in model_setup.m
    
  model_details.mo.vacc_timing='dynamic';
    
    nsamples_cov = 20;
    aeff = [0.25,0.5,0.75];
    nsamples_eff = length(aeff);
    
    adrift = zeros(nsamples_cov,nsamples_eff);
    adrift_koelle = zeros(nsamples_cov,nsamples_eff);
    adrift_plot = zeros(nsamples_cov,nsamples_eff);
    adrift_koelle_plot = zeros(nsamples_cov,nsamples_eff);
    aAR = zeros(nsamples_cov,nsamples_eff);
    aARexp = zeros(nsamples_cov,nsamples_eff);
    aARnai = zeros(nsamples_cov,nsamples_eff);
    aReff = zeros(nsamples_cov,nsamples_eff);
    aratioReff = zeros(nsamples_cov,nsamples_eff);
    atend = zeros(nsamples_cov,nsamples_eff);
    
    afrac_vac = linspace(0,1,nsamples_cov);
    
    iparams{1}.epsilonI = 0.5;
    

    
    for j=1:nsamples_eff
        
        
        for i = 1:nsamples_cov
            
            
            iparams{1}.fracTcellvac = afrac_vac(i);
            iparams{1}.fIE = aeff(j);
            iparams{1}.fIN = aeff(j);
            
            
            [SEIRS_data model_time] = single_simulation(iparams{1},MCMC_details.cd, MCMC_details.df,model_details.mo,MCMC_details.oo);
            [model_inc total_infections drift drift_koelle Reffective pvac gt px infections_stratum ratioReff tend] = SEIRS_simobservables_Tcell(iparams{1},SEIRS_data,model_time,MCMC_details.cd,MCMC_details.df,model_details.mo,MCMC_details.oo);
            
            adrift(i,j) = drift(end);
            adrift_koelle(i,j) = drift_koelle(end);
            aAR(i,j) = total_infections;
            aARexp(i,j) = (infections_stratum(1)+infections_stratum(3))/iparams{1}.N/px(1,1);
            aARnai(i,j) = (infections_stratum(2)+infections_stratum(4))/iparams{1}.N/px(1,2);
            aReff(i,j) = Reffective(end);
            aratioReff(i,j) = ratioReff;
            atend(i,j) = tend;
            
        end
        
        
        adrift_plot(:,j) = adrift(:,j)/adrift(1,j);
        adrift_koelle_plot(:,j) = adrift_koelle(:,j)/adrift_koelle(1,j);
        
        
    end
    
    figure(5)
    
    subplot(2,2,1)
    %ylim([0 1])
    xlabel('fractional coverage with CTL vaccine')
    ylabel('attack rate')
    
    
    set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|--|:')
  
    hold all
    plot(afrac_vac,aAR(:,1),afrac_vac,aAR(:,2),afrac_vac,aAR(:,3),'Linewidth',2,'color',[0 0 0])
    if (plotstratum==1)
        plot(afrac_vac,aARexp(:,1),afrac_vac,aARexp(:,2),afrac_vac,aARexp(:,3),'Linewidth',2,'color','b')
        plot(afrac_vac,aARnai(:,1),afrac_vac,aARnai(:,2),afrac_vac,aARnai(:,3),'Linewidth',2,'color','r')
    end
    
    
    legend('\epsilon_{vI} = 0.25','\epsilon_{vI} = 0.5','\epsilon_{vI} = 0.75','Location','SouthWest');
    
    
    subplot(2,2,2)
    set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|--|:')
  
    hold all
    for j=1:nsamples_eff
        plot(afrac_vac,aReff(:,j),'LineWidth',2)
    end
    legend('\epsilon_{vI} = 0.25','\epsilon_{vI} = 0.5','\epsilon_{vI} = 0.75','Location','NorthWest');
    xlabel('fractional coverage with CTL vaccine')
    ylabel('R_{eff} at end of season')
    
    subplot(2,2,3)
    hold all
    set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|--|:')
    for j = 1:nsamples_eff
        plot(afrac_vac,adrift_koelle_plot(:,j),'Linewidth',2)
    end
    xlabel('fractional coverage with CTL vaccine')
    ylabel('Relative prob of displacement by drift variant (\Delta \sim t^{\lambda})')
    legend('\epsilon_{vI} = 0.25','\epsilon_{vI} = 0.5','\epsilon_{vI} = 0.75','Location','SouthWest');
    
    subplot(2,2,4)
    set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|--|:')
    hold all
    for j = 1:nsamples_eff
        plot(afrac_vac,adrift_plot(:,j),'Linewidth',2)
    end
    xlabel('fractional coverage with CTL vaccine')
    ylabel('Relative prob of displacement by drift variant (\Delta \sim I_{cum})')
    legend('\epsilon_{vI} = 0.25','\epsilon_{vI} = 0.5','\epsilon_{vI} = 0.75','Location','SouthWest');    
    
    print('-r600','-depsc',filename) 
    
    % save to file for re-plotting
    
    txtfilename = strcat('keep',filename);
    
    dlmwrite(txtfilename,[nsamples_eff nsamples_cov])
    dlmwrite(txtfilename,[aAR], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[aARexp], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[aARnai], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[aReff], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[adrift_koelle], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[adrift], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[aratioReff], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[atend],'-append','precision','%10.10f')
    
    fprintf(1,'done\n');
    
elseif strcmp(run_type,'explorevacdist')
    
    % vary vaccine coverage for fixed prior protection and epidemic
    % characteristics
    
    acoverage = [0,0.2,0.2,0.2,0.2,0.2];
    vacc_strat = {'unvac','unvac','vacE','biasvacE','biasvacN','vacN'}; 
    adistTcellvac = [1,1,1,0.5,0.5,1];
    aeff = [0.25,0.5,0.75];
    nsamples_dist = length(adistTcellvac);
    nsamples_eff = length(aeff);
    
    adrift = zeros(nsamples_dist,nsamples_eff);
    adrift_koelle = zeros(nsamples_dist,nsamples_eff);
    adrift_plot = zeros(nsamples_dist,nsamples_eff);
    adrift_koelle_plot = zeros(nsamples_dist,nsamples_eff);
    aAR = zeros(nsamples_dist,nsamples_eff); 
    aARnai = zeros(nsamples_dist,nsamples_eff);
    aARexp = zeros(nsamples_dist,nsamples_eff);
    
    aReff = zeros(nsamples_dist,nsamples_eff);
    
    xplot = linspace(0,length(vacc_strat)-1,length(vacc_strat));
    
   

    % loop over vaccine efficacy
    for j=1:nsamples_eff
        
        % loop over vaccine strategy
        for i = 1:nsamples_dist
            
            model_details.mo.vacc_strat = (vacc_strat{i});
            
            
            iparams{1}.fracTcellvac = acoverage(i);
            iparams{1}.fIE = aeff(j);
            iparams{1}.fIN = aeff(j);
            %iparams{1}.fIN = min(1,aeff(j)*1.5);
            iparams{1}.distTcellvac = adistTcellvac(i);
            
            [SEIRS_data model_time] = single_simulation(iparams{1},MCMC_details.cd, MCMC_details.df,model_details.mo,MCMC_details.oo);
            [model_inc total_infections drift drift_koelle Reffective pvac gt px infections_stratum ratioReff] = SEIRS_simobservables_Tcell(iparams{1},SEIRS_data,model_time,MCMC_details.cd,MCMC_details.df,model_details.mo,MCMC_details.oo);
            
            adrift(i,j) = drift(end);
            adrift_koelle(i,j) = drift_koelle(end);
            aAR(i,j) = total_infections;
            aARexp(i,j) = (infections_stratum(1)+infections_stratum(3))/iparams{1}.N/px(1,1);
            aARnai(i,j) = (infections_stratum(2)+infections_stratum(4))/iparams{1}.N/px(1,2);
            aReff(i,j) = Reffective(end);
            aReff(i,j) = ratioReff; % if want ratio of drifted strain compared to resident strain
            
        end
        
        % calculate drift relative to epidemic without intervention
        adrift_plot(:,j) = adrift(:,j)/adrift(1,j);
        adrift_koelle_plot(:,j) = adrift_koelle(:,j)/adrift_koelle(1,j);
    end
    
    figure(5)
    set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','+|o|d')
    
    
  
    subplot(2,2,1)
    ylim([0.55 0.85])
    xlabel('vaccine distribution strategy')
    ylabel('attack rate')
    
    
    set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','+|o|d')
  
    hold all
    plot(xplot,aAR(:,1),xplot,aAR(:,2),xplot,aAR(:,3))
    plot(xplot,aARexp(:,1),xplot,aARexp(:,2),xplot,aARexp(:,3),'Linewidth',2,'color','b')
    plot(xplot,aARnai(:,1),xplot,aARnai(:,2),xplot,aARnai(:,3),'Linewidth',2,'color','r')
    
    set(gca,'XTickLabel',{'novac','uni','vacE','biasE','biasN','vacN'})
    
    legend('\epsilon_{vI} = 0.25','\epsilon_{vI} = 0.5','\epsilon_{vI} = 0.75','Location','SouthWest');
    
    
    subplot(2,2,2)
    ylim([0.3 0.7])
     set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','+|o|d')
  
  
    hold all
    for j=1:nsamples_eff
        plot(xplot,aReff(:,j))
    end
    legend('\epsilon_{vI} = 0.25','\epsilon_{vI} = 0.5','\epsilon_{vI} = 0.75','Location','NorthWest');
    xlabel('vaccine distribution strategy')
    ylabel('R_{eff} at end of season')
    set(gca,'XTickLabel',{'novac','uni','vacE','biasE','biasN','vacN'})
    
    subplot(2,2,3)
    ylim([0.5 1.1])
    hold all
     set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','+|o|d')
 
    for j = 1:nsamples_eff
        plot(xplot,adrift_koelle_plot(:,j))
    end
    xlabel('vaccine distribution strategy')
    ylabel('Relative prob of displacement by drift variant (\Delta \sim t^{\lambda})')
    legend('\epsilon_{vI} = 0.25','\epsilon_{vI} = 0.5','\epsilon_{vI} = 0.75','Location','SouthWest');
    set(gca,'XTickLabel',{'novac','uni','vacE','biasE','biasN','vacN'})
    
    subplot(2,2,4)
    ylim([0.4 1.1])
     set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','+|o|d')
    
     
    hold all
    for j = 1:nsamples_eff
        plot(xplot,adrift_plot(:,j))
    end
    xlabel('vaccine distribution strategy')
    ylabel('Relative prob of displacement by drift variant (\Delta \sim I_{cum})')
    legend('\epsilon_{vI} = 0.25','\epsilon_{vI} = 0.5','\epsilon_{vI} = 0.75','Location','SouthWest');    
    set(gca,'XTickLabel',{'novac','uni','vacE','biasE','biasN','vacN'});
    
    
    print('-r600','-depsc',filename) 
    
    % save to file for re-plotting
    
    txtfilename = strcat('keep',filename);
    
    dlmwrite(txtfilename,[nsamples_eff nsamples_dist])
    dlmwrite(txtfilename,[aAR], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[aARexp], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[aARnai], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[aReff], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[adrift_koelle], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[adrift], '-append','precision','%10.10f')
    
elseif strcmp(run_type,'exploreboth')
    
    % FIRST chose vaccine coverage and prepandemic immunity in MCMC_setup.m
    
    % ALSO choose drift rate and vaccine action model_setup.m
    
    
    % vary vaccine coverage for fixed prior protection and epidemic
    % characteristics
    ntrials=10;
    model_details.mo.vacc_timing = 'prepandemic';
    model_details.mo.vacc_strat = 'none';
    
    iparams{1}.fracTcellvac = 0.5;
    
    % fix coverage
    acoverage = iparams{1}.fracTcellvac*ones(ntrials);
    % explore fraction of experienced hosts vaccinated
    acoverage_exp = linspace(max(0,(iparams{1}.fracTcellvac-(1-iparams{1}.pop_dist(1)))/iparams{1}.pop_dist(1)),min(1,iparams{1}.fracTcellvac/iparams{1}.pop_dist(1)),ntrials);
    % explore vaccine efficacy in naiive hosts 
    aeff = linspace(0,1,ntrials);
    nsamples_dist = ntrials;
    nsamples_eff = ntrials;
    
    adrift = zeros(nsamples_dist,nsamples_eff);
    adrift_koelle = zeros(nsamples_dist,nsamples_eff);
    adrift_plot = zeros(nsamples_dist,nsamples_eff);
    adrift_koelle_plot = zeros(nsamples_dist,nsamples_eff);
    aAR = zeros(nsamples_dist,nsamples_eff); 
    aARnai = zeros(nsamples_dist,nsamples_eff);
    aARexp = zeros(nsamples_dist,nsamples_eff);
    aReff = zeros(nsamples_dist,nsamples_eff);
      
    % loop over vaccine efficacy
    for j=1:nsamples_eff
        % set vaccine efficacies
            iparams{1}.fIE = aeff(j);
            iparams{1}.fIN = aeff(j);
       
            % loop over vaccine strategy/coverage
        for i = 1:nsamples_dist
            
            % explictly put people in vaccine class for prepandemic option
            % set pop_dist(2) (frac exp vaccinated at start)
            iparams{1}.pop_dist(2) = acoverage_exp(i);
            % set pop_dist(3) (frac nai vaccinated at start)
            iparams{1}.pop_dist(3) = max(0,(acoverage(i,j)-iparams{1}.pop_dist(1)*acoverage_exp(i))/(1-iparams{1}.pop_dist(1)));
            
          
            
            % calculate attack rates, Reff and drift 
            [SEIRS_data model_time] = single_simulation(iparams{1},MCMC_details.cd, MCMC_details.df,model_details.mo,MCMC_details.oo);
            [model_inc total_infections drift drift_koelle Reffective pvac gt px infections_stratum ratioReff] = SEIRS_simobservables_Tcell(iparams{1},SEIRS_data,model_time,MCMC_details.cd,MCMC_details.df,model_details.mo,MCMC_details.oo);
            
            adrift(i,j) = drift(end);
            adrift_koelle(i,j) = drift_koelle(end);
            aAR(i,j) = total_infections;
            aARexp(i,j) = (infections_stratum(1)+infections_stratum(3))/iparams{1}.N/(px(1,1)+px(1,3));
            aARnai(i,j) = (infections_stratum(2)+infections_stratum(4))/iparams{1}.N/(px(1,2)+px(1,4));
            aReff(i,j) = ratioReff;
            
        end
        
        
    end
    
    
    iparams{1}.fracTcellvac = 0;
iparams{1}.pop_dist(2)=0;
iparams{1}.pop_dist(3)=0;
    
%simulate unimpeded epidemic (in order to normalise drift)

    [SEIRS_data model_time] = single_simulation(iparams{1},MCMC_details.cd, MCMC_details.df,model_details.mo,MCMC_details.oo);
            [model_inc total_infections drift drift_koelle Reffective pvac gt px infections_stratum] = SEIRS_simobservables_Tcell(iparams{1},SEIRS_data,model_time,MCMC_details.cd,MCMC_details.df,model_details.mo,MCMC_details.oo);
    
     % calculate drift relative to epidemic without intervention
     adrift_plot = adrift/drift(end);
     adrift_koelle_plot = adrift_koelle/drift_koelle(end);
     
adrift_koelle_plot
drift_koelle(end)

    % save to file for re-plotting
    
    txtfilename = strcat('keep',filename);
    
    dlmwrite(txtfilename,[nsamples_eff nsamples_dist])
   
    dlmwrite(txtfilename,[aAR], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[aARexp], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[aARnai], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[aReff], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[adrift_koelle_plot], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[adrift_plot], '-append','precision','%10.10f')
    
    dlmwrite(txtfilename,[acoverage_exp], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[aeff], '-append','precision','%10.10f')
    
elseif strcmp(run_type,'explorePI')
    
    varyfE = 0;
    varypBcellE = 1;
    
    %FIRST choose vaccine coverage 
    iparams{1}.fracTcellvac = 0.75;
    
    % also choose vaccine efficacy
    
    j=1;
    aeff = [0.25 0.5 0.75];
    %choose vaccine efficacy
    iparams{1}.fIE = aeff(j);
    iparams{1}.fIN = aeff(j);
    
    
    % vary population distribution, fixed vaccine distribution strategy
    % (general population)
    ntrials=10;
    model_details.mo.vacc_timing = 'dynamic';
    model_details.mo.vacc_strat = 'unvac';
    model_details.mo.vacc_action = 'boost';
    
    % choose humoral immunity in MCMC_setup.m (limited by popdist)
   
    
    % range of popdist to trial (this is fraction of the hosts who are
    % influenza experienced)
    apopdist = linspace(0,0.5,ntrials);
    % range of epsilonI to trial (this is the suppression of infectiousness
    % of experienced hosts)
    aepsilonI = linspace(0.01,1.0,ntrials);
    
  
    nsamples_popdist = ntrials;
    nsamples_epsilonI = ntrials;
    
    adrift = zeros(nsamples_popdist,nsamples_epsilonI);
    adrift_koelle = zeros(nsamples_popdist,nsamples_epsilonI);
    adrift_plot = zeros(nsamples_popdist,nsamples_epsilonI);
    adrift_koelle_plot = zeros(nsamples_popdist,nsamples_epsilonI);
    aAR = zeros(nsamples_popdist,nsamples_epsilonI); 
    aARnai = zeros(nsamples_popdist,nsamples_epsilonI);
    aARexp = zeros(nsamples_popdist,nsamples_epsilonI);
    aReff = zeros(nsamples_popdist,nsamples_epsilonI);
      
    % loop over population distribution
    for j=1:nsamples_popdist
        % set vaccine efficacies
        
        if varyfE
            iparams{1}.pop_dist(1) = apopdist(j);
            iparams{1}.pop_dist(2)=0;
            iparams{1}.pop_dist(3)=0;
        elseif varypBcellE
            iparams{1}.pop_dist(1) = 0.8;
            iparams{1}.xI = [apopdist(j) 0 apopdist(j) 0];            
        end
            
       
            % loop over strength of prior immunity
        for i = 1:nsamples_epsilonI
            
            iparams{1}.epsilonI = aepsilonI(i);
            
            
           
            % calculate attack rates, Reff and drift 
            [SEIRS_data model_time] = single_simulation(iparams{1},MCMC_details.cd, MCMC_details.df,model_details.mo,MCMC_details.oo);
            [model_inc total_infections drift drift_koelle Reffective pvac gt px infections_stratum ratioReff] = SEIRS_simobservables_Tcell(iparams{1},SEIRS_data,model_time,MCMC_details.cd,MCMC_details.df,model_details.mo,MCMC_details.oo);
            
            adrift(j,i) = drift(end);
            adrift_koelle(j,i) = drift_koelle(end);
            aAR(j,i) = total_infections;
            aARexp(j,i) = (infections_stratum(1)+infections_stratum(3))/iparams{1}.N/(px(1,1)+px(1,3));
            aARnai(j,i) = (infections_stratum(2)+infections_stratum(4))/iparams{1}.N/(px(1,2)+px(1,4));
            aReff(j,i) = ratioReff;
            
        end
        
        
    end
    
    
    iparams{1}.fracTcellvac = 0;

    
%simulate unimpeded epidemic (in order to normalise drift)

    [SEIRS_data model_time] = single_simulation(iparams{1},MCMC_details.cd, MCMC_details.df,model_details.mo,MCMC_details.oo);
            [model_inc total_infections drift drift_koelle Reffective pvac gt px infections_stratum] = SEIRS_simobservables_Tcell(iparams{1},SEIRS_data,model_time,MCMC_details.cd,MCMC_details.df,model_details.mo,MCMC_details.oo);
    
     % calculate drift relative to epidemic without intervention
     adrift_plot = adrift/drift(end);
     adrift_koelle_plot = adrift_koelle/drift_koelle(end);
     
adrift_koelle_plot;
drift_koelle(end);

    % save to file for re-plotting
    
    txtfilename = strcat('keep',filename);
    
    dlmwrite(txtfilename,[nsamples_epsilonI nsamples_popdist])
   
    dlmwrite(txtfilename,[aAR], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[aARexp], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[aARnai], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[aReff], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[adrift_koelle_plot], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[adrift_plot], '-append','precision','%10.10f')
    
    dlmwrite(txtfilename,[apopdist], '-append','precision','%10.10f')
    dlmwrite(txtfilename,[aepsilonI], '-append','precision','%10.10f')
    
    
else
    
    % Initialise matlab parallel processing pool if using PT
    if strcmp(MCMC_details.df.MixingMode,'ParallelTempering')
        matlabpool open;
    end
    
    % Save the workspace
    save(MCMC_details.filename_for_workspace) 
    
    fprintf(1,'about to run model...\n');
    % Run the model
    
    
    [param_last] = MCMCmixing(MCMC_details.fname, MCMC_details.filename_for_PTswaps, MCMC_details.df, model_details.mo, MCMC_details.oo, ...
        MCMC_details.cd, MCMC_details.pd, MCMC_details.p, iparams, MCMC_details.pi);
    
    % Write key parameters to the 'details' structure for use in data
    % visualisation
    MCMC_details.param_last = param_last;  
    
    % Save the workspace
    save(MCMC_details.filename_for_workspace)
    
    % Shutdown matlab parallel processing pool if using PT
    if strcmp(MCMC_details.df.MixingMode,'ParallelTempering')
        matlabpool close;
    end
    
    fprintf('finished MCMC experiment.\n')
    
    return
end

%quit;
