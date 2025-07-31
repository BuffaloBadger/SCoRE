function [ind, dep, flag, message] = solve_ivodes(ind_0, dep_0...
    , f_var, f_val, derivs_fcn, odes_are_stiff)
%solveIVODEs wrapper for solving initial value ODEs using ode45/ode15s
%   revised 5/7/25
%
%   ind_0 = initial value of the independent variable
%   dep_0 = column vector of initial values of the dependent variables
%   f_var = 0 if the final value of the independent variable is known or
%         = the index of the dependent variable whose final value is known
%   f_val = known final value of either the independent variable or one of 
%       the dependent variables
%   derivs_fcn = handle to a function that returns a column vector 
%       containing the values of the derivatives at values of the 
%       independent and dependent variables that are passed to it as 
%       arguments
%   odes_are_stiff = boolean that is true if the ODEs are stiff
%   ind = column vector containing values of the independent variable
%       spanning the range from the initial value to its final value
%   dep = matrix where each column contains the values of one of the
%       dependent variables at corresponding to the independent 
%       variable values in ind
%   flag = integer flag where
%           1 indicates a solution was found
%           0 indicates that the step size may have been too large 
%               causing the solution to be inaccurate
%          -1 indicates that the final value specified for the
%                dependent variable was not reached
%   message = string that provides details on flag
%
    if f_var == 0 
        % f_val is the final value of the independent variable
        flag = 1;
        message = 'The IVODEs were solved.';
        if odes_are_stiff
            [ind,dep] = ode15s(derivs_fcn,[ind_0 f_val],dep_0);
        else
            [ind,dep] = ode45(derivs_fcn,[ind_0 f_val],dep_0);
        end
    else 
        % f_val is the final value of dependent variable dep(f_var)
        options = odeset('Events',@stop);
        count = 0;
        flag = -1;
        message = 'The final value was not reached.';
        ind_f = ind_0 + 1.0;
        while (count < 10 && flag <= 0)
            count = count + 1;
            % solve the odes
            if odes_are_stiff
               [ind,dep] = ode15s(derivs_fcn,[ind_0 ind_f],dep_0...
                   ,options);
            else
                [ind,dep] = ode45(derivs_fcn,[ind_0 ind_f],dep_0...
                    ,options);
            end
            if ind(end) == ind_f % ind_f wasn't large enough
                flag = -1;
                message = 'The final value was not reached.';
                ind_f = (ind_0 + 1.0)*10^(count);
            elseif ind(end) < 0.1*ind_f % ind_f was too large
                ind_f = (ind_0 + 1.0)/10^(count);
                flag = 0;
                message = 'The step size may be too large.';
            else % ind_f was appropriate and the equations are solved
                flag = 1;
                message = 'The IVODEs were solved.';
            end
        end
    end
    
    % Function that provides the integration stopping criterion
    function [stop_if_this_is_zero, is_terminal, direction] = stop(~...
            ,dep_cur)
        is_terminal = 1;
        direction = 0;
        stop_if_this_is_zero = dep_cur(f_var) - f_val;
    end
end

