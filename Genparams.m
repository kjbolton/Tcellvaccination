function params = Genparams(pd,p)
%
% Arguments: pd <struct>: Model parameter choice structure (see MCMC_*run.m for details)
%            p <struct>: Paramter structure, each field (named by the parameter) containing
%            information for that parameter.
%
% Output: parameters <struct>: Simple parameter structure of textual names
%         storing a <double> for the initial value of the parameter.

for i=1:pd.num_params
    name = pd.names{i};
    params.(name) = p.(name).ival;
end

if ((p.R02.fitted==0) && (p.R03.fitted==0))
    params.R02 = params.R01;
    params.R03 = params.R01;
elseif (p.R03.fitted==0)
    params.R03 = params.R02;
end


% Check if alpha2 and alpha3 are fitted. If not, ensure they equal alpha1
for i = 2:3
    str = ['alpha' num2str(i)]; 
    if (p.(str).fitted == 0)
        temp = params.(str);
        if (temp ~= params.alpha1)
            params.(str) = params.alpha1;
            fprintf(1,'\nWARNING: %s was not set to alpha1 but is not fitted.\n', str)
        end
    end
end
