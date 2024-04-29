function reb_I_6
    %REB_I_6 Calculations for Example I.6 of Reaction Engineering Basics
        % given and known constants
    
        % Given and known constants
        n_A_in = 500.; % mol/h
        V_fluid = 500.; % L
        n_Z_in = 0.; % mol/h
        Cp_vol = 1170.; % cal/L/K
        V_flow = 250.; % L/h
        T_in_K = 423.; % K
        dH_rxn = 18200.; % cal/mol
        k_0 = 1.14E9; % L/mol/h
        E = 16200.; % cal/mol
        Re = 1.987; % cal/mol/K
    
        % residuals function
        function resids = residuals(guess)
            % extract the individual guesses
            n_A_guess = guess(1);
            n_Z_guess = guess(2);
            T_guess_K = guess(3);
    
            % Calculate the concentration of A
            C_A = n_A_guess/V_flow;
    
            % Calculate the rate
            r = k_0*exp(-E/Re/T_guess_K)*C_A^2;
            
            % Evaluate the residuals
            residual_1 = n_A_in - n_A_guess - V_fluid*r;
            residual_2 = n_Z_in - n_Z_guess + V_fluid*r;
            residual_3 = Cp_vol*V_flow*(T_guess_K - T_in_K)...
                + V_fluid*r*dH_rxn;
    
            % Return the residuals
            resids = [residual_1; residual_2; residual_3];
        end
    
        % reactor model
        function soln = unknowns()
    
            % set the initial guess
            init_guess = [n_A_in/2; n_A_in/2; T_in_K - 1.];
            
            % solve the ATEs
            [soln, flag, message] = solve_ates(@residuals, init_guess);
        
            % check that the solution was found
            if flag <= 0
                disp(' ')
                disp(['WARNING: The ATE solver did not converge: ',message])
            end
        end
    
        % perform the analysis
    
        % call the model function
        soln = unknowns();
    
        % Extract the solution
        n_A_out = soln(1);
        n_Z_out = soln(2);
        T_out_K = soln(3);
    
        % tabulate the results
        item = ["Flow Rate of A";"Flow Rate of Z";"Temperature"];
        value = [n_A_out;n_Z_out;T_out_K-273.15];
        units = ["mol h^-1^";"mol h^-1^";"°C"];
        results_table = table(item,value,units);
    
        % display the results
        disp(' ')
        disp(results_table)
    
        % save the results
        writetable(results_table,'results.csv');
    end