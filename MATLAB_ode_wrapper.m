function [junk SEIRS_data] = MATLAB_ode_wrapper(mtimes, SEIRS_init, params, mo, oo, cd, df, bf_array, options)

[junk SEIRS_data] = ode23(@(time,curr_state) mo.model_rhs(time,curr_state,params,cd.N,mo,oo,cd,df,bf_array), ...
        mtimes, SEIRS_init, options);
end
