function [beta_factor] = fbeta_factor(cd, mo, params,SEIRS,epsilon,alpha,R0)

%calculate scaling of force of infection matrix for 2 strata model from R0
%and scaling parameters


if strcmp(mo.reinfection,'simpleSIR') 
    
   beta_factor = R0/cd.N/params.Ti/alpha;
 
elseif strcmp(mo.reinfection,'simpleSEIR')
    
    beta_factor = R0/cd.N/params.Ti/alpha;
    
else

m = mo.stratifications;
ns = mo.nstates;
stgmax = mo.nstgs;
nmultiple = mo.nmultiple;

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
cS2 =cTQ2+m;
cE1Q = cS2 + m;
cE2Q = cE1Q+m;
cE1TQ = cE2Q +m;
cE2TQ = cE1TQ + m;
cE1TQ2 = cE2TQ + m;
cE2TQ2 = cE1TQ2 + m;
cCUMA = cE1TQ2 + m;

S = zeros(m);

if strcmp(mo.vacc_timing,'prepandemic')
%unvaccinated any already vaccinated hosts to calculate appropriate beta
params.pop_dist(2) = 0;
params.pop_dist(3) = 0;
end


SEIRS_init = MCMCinit(params,cd.N,mo);

for j = 1:m
    
    S(j,1) = SEIRS_init(cS(j)) + SEIRS_init(cS2(j)) + SEIRS_init(cE1(j)) + SEIRS_init(cE2(j))  + SEIRS_init(cI(j)) + SEIRS_init(cA(j)) + ...
        SEIRS_init(cR(j)) + SEIRS_init(cRA(j)) + SEIRS_init(cT(j)) + SEIRS_init(cT2(j)) + SEIRS_init(cTA(j)) + SEIRS_init(cTA2(j)) +...
        SEIRS_init(cQ(j)) + SEIRS_init(cTQ(j)) + SEIRS_init(cTQ2(j))  + SEIRS_init(cE1Q(j)) + SEIRS_init(cE2Q(j)) + ...
        SEIRS_init(cE1TQ(j)) + SEIRS_init(cE2TQ(j)) + SEIRS_init(cE1TQ2(j)) + SEIRS_init(cE2TQ2(j));
    
    %S(j,1) = S(j,1) + SEIRS_init(cP(j));
  
    
end



if (m>2)
    % calculate for unvaccinated population so add in any hosts that have
    % been vaccinated pre-pandemic for that mo.vacc_timing
    if strcmp(mo.vacc_timing,'prepandemic')
       S(1,1) = S(1,1) +S(3,1);
       S(2,1) = S(2,1) +S(4,1);
    end
    S(1,2) = S(1,1);
    S(2,2) = S(2,1);
end


if (m==4)
   % calculate beta_factor for population without interventions (T cell vaccine)
   % so equations drop down to 2 dimensions
   Sunvac = S(1:2,1:2);
  
end

M = max(eig(Sunvac.*epsilon(1:2,1:2)));



fprintf(1,'M norm for beta (divided by N) = %e\n',M/cd.N)

if strcmp(mo.transmission,'fixReff')
    beta_factor = R0/M/params.Ti;
elseif strcmp(mo.transmission,'fixbeta')
    beta_factor = R0/params.N/params.Ti;
else
    error('have not chosen valid transmission option for model');
end

Sunvac
epsilon(1:2,1:2)

fprintf(1,'beta is %e\n',beta_factor);

end

return
