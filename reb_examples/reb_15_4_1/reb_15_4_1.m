function reb_15_4_1
    % constants available to all functions
	% given
    V_cstr = 350.0; % L
    V_pfr = 350.0; % L
    CA_0 = 2.5; % mol /L
    Vdot = 100.0; % L /min
    T_0 = 38 + 273.15; % K
    dH1 = -21500; % cal /mol
    dH2 = -24000; % cal /mol
    Cp = 1.0E3; % cal /L /K
    k0_1 = 1.2E5; % /min
    E_1 = 9100; % cal/mol
    k0_2 = 2.17E7; % L /mol /min
    E_2 = 13400; % cal /mol
	% known
    R = 1.987; % cal /mol
	% calculated
    nDotA_0 = Vdot*CA_0;

    % make CSTR inlet stream available to all functions
    nDotA_in = nan;
    nDotD_in = nan;
    nDotU_in = nan;
    T_in = nan;

    % residuals function
    function resids = residuals(guess)
        % extract the individual guesses
        nDotA_out = guess(1);
        nDotD_out = guess(2);
        nDotU_out = guess(3);
        T_out = guess(4);

        % calculate the rates
        r1 = k0_1*exp(-E_1/R/T_out)*nDotA_out/Vdot;
        r2 = k0_2*exp(-E_2/R/T_out)*(nDotA_out/Vdot)^2;

        % evaluate the residuals
        epsilon_1 = nDotA_in - nDotA_out + (-r1 -r2)*V_cstr;
        epsilon_2 = nDotD_in - nDotD_out + r1*V_cstr;
        epsilon_3 = nDotU_in - nDotU_out + r2*V_cstr;
        epsilon_4 = -Vdot*Cp*(T_out - T_in) - V_cstr*(r1*dH1 + r2*dH2);

        resids = [epsilon_1; epsilon_2; epsilon_3; epsilon_4];
    end
    
    % cstr model
    function soln = unknowns(initial_guess)
        
        % solve the ATEs
        [soln, flag, message] = solve_ates(@residuals, initial_guess);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp(['WARNING: The ATE solver did not converge: ',message])
        end
    end

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract the dependent variables for this integration step
        nDotA = dep(1);
        T = dep(4);

        % calculate rates
        r1 = k0_1*exp(-E_1/R/T)*nDotA/Vdot;
        r2 = k0_2*exp(-E_2/R/T)*(nDotA/Vdot)^2;

        % evaluate the derivatives
        dnAdV = -r1 -r2;
        dnDdV = r1;
        dnUdV = r2;
        dTdV = -(r1*dH1 + r2*dH2)/Vdot/Cp;

        % return the derivatives
        derivs = [dnAdV; dnDdV; dnUdV; dTdV];
    end

    % pfr model
    function [V, nDotA, nDotD, nDotU, T] = profiles(dep_0)
        % set the initial values
        ind_0 = 0.0;

        % define the stopping criterion
        f_var = 0;
        f_val = V_pfr;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [V, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract and return the dependent variable profiles
        nDotA = dep(:,1);
        nDotD = dep(:,2);
        nDotU = dep(:,3);
        T = dep(:,4);
    end

    % function that performs the analysis
	function perform_the_analysis()

        % case a
        disp('Starting case a')

        nDotA_in = nDotA_0;
        nDotD_in = 0.0;
        nDotU_in = 0.0;
        T_in = T_0;

        % solve cstr design equations
        initial_guess = [nDotA_in*.9; nDotA_in/10; nDotA_in/10; T_in + 20];
        cstr_soln = unknowns(initial_guess);
        nDotA_1 = cstr_soln(1);
        nDotD_1 = cstr_soln(2);
        nDotU_1 = cstr_soln(3);
        T_1 = cstr_soln(4);

        % solve pfr design equations
        dep_0 = [nDotA_1; nDotD_1; nDotU_1; T_1];
        [~, nDotA, nDotD, nDotU, T] = profiles(dep_0);

        % calculate the other quantities of interest
        fA_case_a = 100.0*(nDotA_0 - nDotA(end))/nDotA_0;
        sel_case_a = nDotD(end)/nDotU(end);
        T_case_a = T(end) - 273.15;

        % case b
        disp('Starting case b')

        % solve the pfr design equations
        dep_0 = [nDotA_0; 0.0; 0.0; T_0];
        [~, nDotA, nDotD, nDotU, T] = profiles(dep_0);

        % solve the cstr design equations
        nDotA_in = nDotA(end);
        nDotD_in = nDotD(end);
        nDotU_in = nDotU(end);
        T_in = T(end);
        initial_guess = [nDotA_in*0.9; nDotA_in/10; nDotA_in/10; T_in + 20];
        cstr_soln = unknowns(initial_guess);
        nDotA_2 = cstr_soln(1);
        nDotD_2 = cstr_soln(2);
        nDotU_2 = cstr_soln(3);
        T_2 = cstr_soln(4);

        % calculate the other quantities of interest
        fA_case_b = 100.0*(nDotA_0 - nDotA_2)/nDotA_0;
        sel_case_b = nDotD_2/nDotU_2;
        T_case_b = T_2 - 273.15;

        % tabulate the results
        item = ["case a conversion";"case a selectivity";"case a T"
            "case b conversion";"case b selectivity";"case b T"];
        value = [fA_case_a; sel_case_a; T_case_a; fA_case_b; sel_case_b
            T_case_b];
        units = ["%";"mol D per mol U"; "°C"; "%";"mol D per mol U"
            "°C"];
        results_table = table(item,value,units);

        % display the results
        disp(results_table)

        % save the results
        writetable(results_table,'results.csv');
    end

    % perform the analysis
    perform_the_analysis()

end