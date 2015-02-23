function [Reffective] = fReff(cd,mo,params,SEIRS,epsilon,beta)
   
% instantaneous growth rate (if virus present)
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




epsilon = epsilon*beta;

Sstates = zeros(m,m);

for strata = 1:m
    for stg = 1:stgmax
        for k = 1:nmultiple
            loop = strata + m*(k-1) + (stg-1)*m*nmultiple;
            Sstates(strata,1) = Sstates(strata,1) + SEIRS(cS(strata)) + SEIRS(cS2(strata));
        end
    end
    Sstates(strata,:) = Sstates(strata,1);
end

Reffective = max(eig(epsilon.*Sstates))*params.Ti;
