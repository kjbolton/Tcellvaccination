function [SEIRS model_time seeding_flag] = single_simulation(params,cd,df,mo,oo)

if strcmp(mo.vacc_strat,'vacE') || strcmp(mo.vacc_strat,'vacN')
    if ((params.pop_dist(1)==1 || params.pop_dist(1)==0) && params.distTcellvac < 1)
        error('does not make sense to bias distribution of vaccines when only one type of host')
    end
end

[SEIRS_init,model_time] = MCMCinit(params,cd.N,mo);

[SEIRS model_time seeding_flag] = mo.model_SEIRS_func(params,SEIRS_init,model_time,mo,oo,cd,df);

return
