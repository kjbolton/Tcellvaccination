function [SEIRS_init] = introduce_seed(SEIRS_init,wave_no,mo,params)

%Introduce an extra infected person into population randomly.
%Input arguments SEIRS_init<16*6 double>: array of initial states at
%                                          seeding time
%                 wave_no <int>: integer indicating wave number wrt data
%                 (and therefore flow timings in 7 stage model)
%
%Output          SEIRS_init<16*6 double>: new array of initial states at
%                                         seeding time


m = mo.stratifications; % number of strata
stgmax = mo.nstgs; %number of reinfection stages
nmultiple = mo.nmultiple; % number of subloops per stage
ns = mo.nstates; %number of stages per loop

% relevant state indices
cS = [];
for j = 1:stgmax*nmultiple
    cS = [cS linspace((j-1)*ns*m+1,((j-1)*ns+1)*m,m)];
end
cE1=cS+m;
cE2=cE1+m;
cI=cE2+m;

fprintf(1,'************ INTRODUCING SEED **********************\n');

% choose strata in which infected person belongs (i.e. naiive or
% experienced)

iseed = 1; %seed into this strata
if (iseed>m) 
    iseed = m;
end
seeding = '';
%'all' to seed each strata equally (with less than one person...), 
%'random_stage' to choose a stage randomly and seed into the iseed(th) strata, 
%'newI_stage' to seed into stage corresponding to people getting their


%2/11/11 - introduce seedsize (<1) to adjust for fact that beta is scaled by 1/alpha to account for infections
%by asympomatically infected hosts.
seedsize = zeros(m,1);
seedsize(1,1) = params.I0*params.N;

if (1==1) 
    
    

SEIRS_init(cI(1)) = SEIRS_init(cI(1)) + seedsize(1);

else

if(1==0)
    
    %adjust size of seed according to symptomatic proportion (since beta is
    %scaled up accordingly)
    % don't need this now using I+A to calculate itot and FOI
    if(m==2)
        seedsize(1,1) = seedsize(1,1)*params.epsilonALPHA;
    end
    
    if(wave_no==1)
        seedsize = seedsize*params.alpha1;
    elseif(wave_no==2)
        seedsize = seedsize*params.alpha2;
    else
        seedsize = seedsize*params.alpha3;
    end
end


%first infection (default).

if strcmp(mo.reinfection,'reinfection')
    if (wave_no == 1) 
        if strcmp(seeding,'random_stage')
            if strcmp(seeding,'all')
                for i=1:m
                    SEIRS_init(cI(i)) = SEIRS_init(cI(i)) + seedsize(i)/m;
                end
            else
                SEIRS_init(cI(iseed)) = SEIRS_init(cI(iseed)) + seedsize(i);
            end
        else
            %seed into stage 1
            SEIRS_init(cI(iseed)) = SEIRS_init(cI(iseed)) + seedsize(iseed);
        end
    elseif (wave_no == 2) 
        if strcmp(seeding,'random_stage')
            %choose randomly between stages 2 and 5 to seed
            seed_loop = randint(1,1,[0,1]);
            if (seed_loop == 0)
                if strcmp(seeding,'all')
                    for i=1:m
                        SEIRS_init(cI(3*m+i)) = SEIRS_init(cI(3*m+i)) + seedsize(i)/m;
                    end
                else
                    SEIRS_init(cI(3*m+iseed)) = SEIRS_init(cI(3*m+iseed)) + seedsize(i);
                end
            else
                if strcmp(seeding,'all')
                    for i=1:m
                        SEIRS_init(cI(12*m+i)) = SEIRS_init(cI(12*m+i)) + seedsize(i)/m;
                    end
                else
                    SEIRS_init(cI(12*m+iseed)) = SEIRS_init(cI(12*m+iseed)) + seedsize(i);
                end
            end
        else
            %seed into stage 5
            SEIRS_init(cI(12*m+iseed)) = SEIRS_init(cI(12*m+iseed)) + seedsize(iseed);
        end
    elseif (wave_no == 3) 
        if strcmp(seeding,'random_stage')
            %choose randomly between stages 4, 8, 10, 12 to seed
            seed_loop = randint(1,1,[0,3]);
            if (seed_loop == 0)
                % seed into stage 4
                if strcmp(seeding,'all')
                    for i=1:m
                        SEIRS_init(cI(9*m+i)) = SEIRS_init(cI(9*m+i)) + seedsize(i)/m;
                    end
                else
                    SEIRS_init(cI(9*m+iseed)) = SEIRS_init(cI(9*m+iseed)) + seedsize(iseed);
                end
            end
            if (seed_loop == 1)
                %seed into stage 8
                if strcmp(seeding,'all')
                    for i=1:m
                        SEIRS_init(cI(21*m+i)) = SEIRS_init(cI(21*m+i)) + seedsize(i)/m;
                    end
                else
                    SEIRS_init(cI(21*m+iseed)) = SEIRS_init(cI(21*m+iseed)) + seedsize(iseed);
                end
            end
            if (seed_loop == 2)
                %seed into stage 10
                if strcmp(seeding,'all')
                    for i=1:m
                        SEIRS_init(cI(27*m+i)) = SEIRS_init(cI(27*m+i)) + seedsize(i)/m;
                    end
                else
                    SEIRS_init(cI(27*m+iseed)) = SEIRS_init(cI(27*m+iseed)) + seedsize(iseed);
                end
            end
            if (seed_loop == 3)
                %seed into stage 12
                if strcmp(seeding,'all')
                    for i=1:m
                        SEIRS_init(cI(33*m+i)) = SEIRS_init(cI(33*m+i)) + seedsize(i)/m;
                    end
                else
                    SEIRS_init(cI(33*m+iseed)) = SEIRS_init(cI(3*m+iseed)) + seedsize(iseed);
                end
            end
        else
            %seed into stage 12
            SEIRS_init(cI(33*m+iseed)) = SEIRS_init(cI(33*m+iseed)) + seedsize(iseed);
        end
    else
        %error('ERROR: invalid wave number in introduce_seed.m');
    end
end

if strcmp(mo.reinfection,'snake')
    error('ERROR: have not re-written this part of the code new state indices with secondary infections');
    if (wave_no == 1)
        if strcmp(seeding,'all') 
            for i=1:m
                SEIRS_init(cI(i)) = SEIRS_init(cI(i)) + seedsize(i)/m;
                %SEIRS_init(cS(i)) = SEIRS_init(cS(i)) - 1./m;
            end
        else
            SEIRS_init(cI(iseed)) = SEIRS_init(cI(iseed)) + seedsize(iseed);
            %SEIRS_init(cS(iseed)) = SEIRS_init(cS(iseed)) - 1;
        end
    end
    if(wave_no == 2) 
        if strcmp(seeding,'all')
            for i=1:m
                SEIRS_init(cI(m+i)) = SEIRS_init(cI(m+i)) + seedsize(i)/m;
                %SEIRS_init(cS(i)) = SEIRS_init(cS(i)) - 1./m;
            end
        else
            SEIRS_init(cI(m+iseed)) = SEIRS_init(cI(m+iseed)) + seedsize(iseed);
            %SEIRS_init(cS(m+iseed)) = SEIRS_init(cS(m+iseed)) - 1;
        end
    end
    if (wave_no == 3)
        if strcmp(seeding,'all')
            for i = 1:m
                SEIRS_init(cI(2*m+i)) = SEIRS_init(cI(2*m+i)) + seedsize(i)/m;
                %SEIRS_init(cS(2*m+i)) = SEIRS_init(cS(2*m+1)) - 1./m;
            end
        else
            SEIRS_init(cI(2*m+iseed)) = SEIRS_init(cI(2*m+iseed)) + seedsize(iseed);
            %SEIRS_init(cS(2*m+iseed)) = SEIRS_init(cS(2*m+iseed)) - 1;
        end
    end
end
end
