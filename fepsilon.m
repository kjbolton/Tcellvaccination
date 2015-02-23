function [epsilon] = fepsilon(params,mo)

m = mo.stratifications;


epsilon = ones(m,m);
if (mo.symbolic ==1)
    epsilon = sym(epsilon);
end

% define 2*2 array with of mixing probabilities within and between stratum
% only need 2*2 even when m=4 since assuming that mixing propensity doesn't
% change with T cell vaccinations status

mixing  = ones(2,2);

mixing(1,1) = params.beta11mix;
mixing(1,2) = params.beta12mix;
mixing(2,1) = params.beta21mix;
mixing(2,2) = params.beta22mix;


if (m==2)
    epsilon(1,1) = params.epsilonI*params.epsilonS/params.epsilonNU;
    epsilon(2,1) = params.epsilonI;
    epsilon(1,2) = params.epsilonS/params.epsilonNU;
    epsilon(2,2) = 1;
end

if (m==4)
    
    epsilonI = params.epsilonI;
    epsilonS = params.epsilonS;
    epsilonNU = params.epsilonNU;
    
    if strcmp(mo.vacc_action,'boost')
        epsilonSE = params.fSE*params.epsilonS;
    elseif strcmp(mo.vacc_action,'sat')
        epsilonSE = min(params.fSE,params.epsilonS);
    end
    epsilonSN = params.fSN;
    if strcmp(mo.vacc_action,'boost')
        epsilonIE = params.fIE*params.epsilonI;
    elseif strcmp(mo.vacc_action,'sat')
        epsilonIE = min(params.fIE,params.epsilonI);
    end
    epsilonIN = params.fIN;
    if strcmp(mo.vacc_action,'boost')
        epsilonNUE = params.fNUE*params.epsilonNU;
    elseif strcmp(mo.vacc_action,'sat')
        epsilonNUE = min(params.fNUE,params.epsilonNU);
    end
    epsilonNUN = params.fNUN;
    
    epsilon(1,1) = epsilonI*epsilonS/epsilonNU*mixing(1,1);
    epsilon(2,1) = epsilonI*mixing(2,1);
    epsilon(3,1) = epsilonI*epsilonSE/epsilonNUE*mixing(1,1);
    epsilon(4,1) = epsilonI*epsilonSN/epsilonNUN*mixing(2,1);
    
    epsilon(1,2) = epsilonS/epsilonNU*mixing(1,2);
    epsilon(2,2) = 1*mixing(2,2);
    epsilon(3,2) = epsilonSE/epsilonNUE*mixing(1,2);
    epsilon(4,2) = epsilonSN/epsilonNUN*mixing(2,2);
    
    epsilon(1,3) = epsilonIE*epsilonS/epsilonNU*mixing(1,1);
    epsilon(2,3) = epsilonIE*mixing(2,1);
    epsilon(3,3) = epsilonIE*epsilonSE/epsilonNUE*mixing(1,1);
    epsilon(4,3) = epsilonIE*epsilonSN/epsilonNUN*mixing(2,1);
    
    epsilon(1,4) = epsilonIN*epsilonS/epsilonNU*mixing(1,2);
    epsilon(2,4) = epsilonIN*mixing(2,2);
    epsilon(3,4) = epsilonIN*epsilonSE/epsilonNUE*mixing(1,2);
    epsilon(4,4) = epsilonIN*epsilonSN/epsilonNUN*mixing(2,2);
    
    
end
