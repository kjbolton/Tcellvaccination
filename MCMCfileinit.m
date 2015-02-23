function filename_base = MCMCfileinit(data_set,extra_info,itt_tot)
% MCMCfileinit - provide filename for the model run
%
% Arguments: data_set <string>: Name of data_set
%            extra_info <string>: Freeform text for run identification
%            itt_tot <int>: Number of iterations performed
%
% Output: filename_base <string>: Base for filename used for model run
 
% Prepare file for writing...
tmp_clock = clock;
filename_base = strcat(data_set,'-',extra_info, ...
    num2str(itt_tot),'-',date,'-', ...
    num2str(tmp_clock(4)),'-', num2str(tmp_clock(5)));
end
