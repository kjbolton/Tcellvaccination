function SEIRS_data_mtimes = check_data_array(SEIRS_data,mtimes)

%Input arguments SEIRS_data <mtimes X (16*7) data array>: 2d array of states
%                           at times mtimes
%                 times <double array>: 1d array of state output times
%Output          SEIRS_data <mtimes  X (16*7) double>: new initial states


if (length(mtimes)==2)
    fprintf(1,'modifying SEIRS_data in check_data_array\n');
    no_states = size(SEIRS_data,2);
    SEIRS_data_mtimes = zeros(2,no_states);
    SEIRS_data_mtimes(1,:) = SEIRS_data(1,:);
    SEIRS_data_mtimes(2,:) = SEIRS_data(end,:);
else
    SEIRS_data_mtimes = SEIRS_data;
end
