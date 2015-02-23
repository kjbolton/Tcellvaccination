function [avTcell vBcell] = vacrate(params,SEIRS,px,time,mo)

tol = 1.e-5; % don't vaccinate when fraction left in pool is less than this

m = mo.stratifications;
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

vTcell = 0;
avTcell = zeros(1,2);

fracunvac = sum(px(1:2)); % all hosts in strata 1 and 2 are unvaccinated by definition

if (time >= params.T_Tcellvac) && (time < params.T_Tcellvac + params.periodTcellvac)  
    vTcell = params.fracTcellvac/params.periodTcellvac;
    % assume T cell vaccine is only given to those people who haven't
    % already recieved one
    if strcmp(mo.vacc_strat,'unvac')
        if (fracunvac > tol) 
                avTcell(1) = vTcell/fracunvac;
                avTcell(2) = vTcell/fracunvac;
        end
    end
    fracunvacE = px(1); % this is fraction of the total population
    fracunvacEnorm = px(1)/fracunvac;
    fracunvacN = px(2);
    fracunvacNnorm = px(2)/fracunvac;
    
    if strcmp(mo.vacc_strat,'biasvacN')
        if (fracunvacE>tol)
            avTcell(1) = vTcell*params.distTcellvac/fracunvac;
        else
            avTcell(1) = 0;
        end
        if(fracunvacN>tol)
            avTcell(2) = vTcell*(1-params.distTcellvac*fracunvacEnorm)/(fracunvacN);     
        else
            avTcell(2) = 0;
        end        
    end
    if strcmp(mo.vacc_strat,'vacN')
        % vaccinated all naiive hosts first, then the rest
        if (px(3) + px(4)) < params.fracTcellvac && (px(2)>tol)
            avTcell(1) = 0;
            avTcell(2) = vTcell/fracunvacN/(1-params.pop_dist(1));
        else
            if (px(3) + px(4)) < params.fracTcellvac && px(1)>tol
                avTcell(1) = vTcell/fracunvacE/params.pop_dist(1);
            end
        end
        
    end
    if strcmp(mo.vacc_strat,'biasvacE')
        if (fracunvacE>tol)
            avTcell(1) = vTcell*(1-params.distTcellvac*fracunvacNnorm)/(fracunvacE);
        else
            avTcell(1) = 0;
        end
        if(fracunvacN>tol)
            avTcell(2) = vTcell*params.distTcellvac/fracunvac;
            
        else
            avTcell(2) = 0;
        end       
    end
    if strcmp(mo.vacc_strat,'vacE')
        % vaccinated all experienced hosts first, then the rest
        if (px(3) + px(4)) < params.fracTcellvac && px(1)>tol
            avTcell(1) = vTcell/fracunvacE/params.pop_dist(1);
            avTcell(2) = 0;
        else
            if (px(3) + px(4)) < params.fracTcellvac && px(2)>tol
                avTcell(2) = vTcell/fracunvacN/(1-params.pop_dist(1));
            end
        end
    end
end

vBcell = 0;

if (time >= params.T_Bcellvac) && (time< params.T_Bcellvac + params.periodBcellvac)  
    vBcell = params.fracBcellvac/params.periodBcellvac;
    if strcmp(mo.vacc_strat,'unvac')
    Nunvac = 0;
    for i = 1:m
        Nunvac = Nunvac + px(i)*params.N - SEIRS(cRA(i)) - SEIRS(cTA(i)) - SEIRS(cTA2(i));
    end
    if Nunvac > 0
        vBcell = vBcell*params.N/Nunvac;
    else
        vBcell = 0;
    end
    end
end

if strcmp(mo.vaccination,'none') || strcmp(mo.vaccination,'Bcell')
    vTcell = 0;
    aTcell = aTcell*0;
end
if strcmp(mo.vaccination,'none')|| strcmp(mo.vaccination,'Tcell')
   vBcell = 0; 
end


if vBcell>1
    error('not sure you can vaccinate at this rate with approximations using to divert vaccinated infected individuals')
end


