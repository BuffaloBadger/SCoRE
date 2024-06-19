function reb_15_4_3
%REB_15_4_3 Calculations for Example 15.4.3 of Reaction Engineering Basics
    % constants available to all functions
    % given
    T_0 = 30 + 273.15; % K
    CA_0 = 1.0; % mol /l
    CB_0 = 1.2; % mol /l
    nDotY_0 = 0.0;
    nDotZ_0 = 0.0;
    Vdot = 75; % l /min
    k0 = 8.72E5; % l /mol /min
    E = 7200; % cal /mol
    dH = -10700; % cal /mol
    Cp = 1.0; % cal /g /K
    rho = 1.0E3; % g /l
    fA = 0.9;
    % known
    R = 1.987; % cal /mol /K
    % calculated
    nDotA_0 = CA_0*Vdot;
    nDotB_0 = CB_0*Vdot;
    nDotA_2 = nDotA_0*(1-fA);

    % make V_R2 available to all functions
    V_R2 = nan;

    % reactor system model
    function soln = unknowns(init_guess)
        
        % solve the ATEs
        [soln, flag, message] = solve_ates(@residuals, init_guess);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp(['WARNING: The ATE solver did not converge: ',message])
        end
    end

    % residuals function
    function resids = residuals(guess)
        % extract the individual guesses
        nDotA_1 = guess(1);
        nDotB_1 = guess(2);
        nDotY_1 = guess(3);
        nDotZ_1 = guess(4);
        T_1 = guess(5);
        V_R1 = guess(6);
        nDotB_2 = guess(7);
        nDotY_2 = guess(8);
        nDotZ_2 = guess(9);
        T_2 = guess(10);

        % calculate the rates
        r_R1 = k0*exp(-E/R/T_1)*nDotA_1*nDotB_1/Vdot^2;
        r_R2 = k0*exp(-E/R/T_2)*nDotA_2*nDotB_2/Vdot^2;

        % evaluate the residuals
        eps1 = nDotA_0 - nDotA_1 - V_R1*r_R1;
        eps2 = nDotB_0 - nDotB_1 - V_R1*r_R1;
        eps3 = nDotY_0 - nDotY_1 + V_R1*r_R1;
        eps4 = nDotZ_0 - nDotZ_1 + V_R1*r_R1;
        eps5 = rho*Vdot*Cp*(T_1 - T_0) + V_R1*r_R1*dH;
        eps6 = nDotA_1 - nDotA_2 - V_R2*r_R2;
        eps7 = nDotB_1 - nDotB_2 - V_R2*r_R2;
        eps8 = nDotY_1 - nDotY_2 + V_R2*r_R2;
        eps9 = nDotZ_1 - nDotZ_2 + V_R2*r_R2;
        eps10 = rho*Vdot*Cp*(T_2 - T_1) + V_R2*r_R2*dH;

        % return the residuals
        resids = [eps1; eps2; eps3; eps4; eps5; eps6; eps7; eps8; eps9
            eps10];
    end

    % function that performs the analysis
    function perform_the_analysis()

        % choose a range of values for V_R2
        V_R2_range = linspace(50.0,70.0)';

        % set the initial guess
        init_guess = [nDotA_0*0.11
            nDotB_0 - nDotA_0*0.89
            nDotA_0*0.89
            nDotA_0*0.89
            T_0 + 10.0
            100.
            nDotB_0 - nDotA_0*0.9
            nDotA_0*0.9
            nDotA_0*0.9
            T_0 + 10.0];

        % allocate storage for the corresponding values of V_R1
        V_R1 = nan(100,1);

        % calculate the corresponding values of V_R1
        for i = 1:100
            % set the volume of reactor R2
            V_R2 = V_R2_range(i);
    
            % solve the reactor design equations
            soln = unknowns(init_guess);

            % save the result
            V_R1(i) = soln(6);
    
            % use the result as the next initial guess
            init_guess = soln;
        end

        % plot the results
        figure; % residuals vs T
        plot(V_R2_range,V_R1 + V_R2_range,'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Reactor 2 Volume (l)','FontSize', 14)
        ylabel('Total Reactor Volume (l)','FontSize', 14)
        saveas(gcf,"volume_plot.png")

        % find the minimum volume
        [min_total_V, minIndex] = min(V_R1 + V_R2_range);
        opt_V_R2 = V_R2_range(minIndex);
        opt_V_R1 = min_total_V - opt_V_R2;

        % tabulate the results
        item = ["Minimum Total Volume";"Optimum Reactor 1 Volume"
            "Optimum Reactor 2 Volume"];
        value = [min_total_V; opt_V_R1; opt_V_R2];
        units = ["L";"L";"L"];
        results_table = table(item,value,units);

        % display the results
        disp(' ')
        disp(results_table)

        % save the results
        writetable(results_table,'results.csv');
    end

    % perform the analysis
    perform_the_analysis();
end