function reb_13_7_4
%REB_13_7_4 Calculations for Example 13.7.4 of Reaction Engineering Basics
    % constants available to all functions
	% given
    dH = -9120.0; % cal /mol
    K0 = 0.132;
    Cp_CO = 29.3/4.184; % cal /mol /K
    Cp_H2O = 34.3/4.184; % cal /mol /K
    Cp_CO2 = 41.3/4.184; % cal /mol /K
    Cp_H2 = 29.2/4.184; % cal /mol /K
    Cp_I = 40.5/4.184; % cal /mol /K
    k0 = 0.0354; % mol /cm^3 /min /atm^2
    E = 9740.0; % cal /mol
    P = 26.0; % atm
    Tin = 320 + 273.15;
    nCO_in = 1.0; % mol /h
    nCO2_in = 0.359; % mol /h
    nH2_in = 4.44; % mol /h
    nI_in = 0.18; % mol /h
    fCO = 0.55;
	% known
    Re = 1.987; % cal /mol /K

    % make the inlet molar flow rate of H2O available to all functions
    nH2O_in = nan;

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract the dependent variables for this integration step
        nCO = dep(1);
        nH2O = dep(2);
        nCO2 = dep(3);
        nH2 = dep(4);
        nI = dep(5);
        T = dep(6);

        % calculate [quantities in the derivatives]
        k = k0*exp(-E/Re/T);
        K = K0*exp(-dH/Re/T);
        n = nCO + nH2O + nCO2 + nH2 + nI;
        P_CO = nCO/n*P;
        P_H2O = nH2O/n*P;
        P_CO2 = nCO2/n*P;
        P_H2 = nH2/n*P;
        r = k*(P_CO*P_H2O - P_CO2*P_H2/K);

        % evaluate the derivatives
        dnCOdV = -r;
        dnH2OdV = -r;
        dnCO2dV = r;
        dnH2dV = r;
        dnIdV = 0.0;
        dTdV = -r*dH/(nCO*Cp_CO + nH2O*Cp_H2O + nCO2*Cp_CO2 ...
            + nH2*Cp_H2 + nI*Cp_I);

        % return the derivatives
        derivs = [dnCOdV; dnH2OdV; dnCO2dV; dnH2dV; dnIdV; dTdV];
    end

    % reactor model
    function [V, nCO, nH2O, nCO2, nH2, nI, T] = profiles()
        % set the initial values

        ind_0 = 0.0;
        dep_0 = [nCO_in; nH2O_in; nCO2_in; nH2_in; nI_in; Tin];

        % define the stopping criterion
        f_var = 1;
        f_val = nCO_in*(1-fCO);
        
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
        nCO = dep(:,1);
        nH2O = dep(:,2);
        nCO2 = dep(:,3);
        nH2 = dep(:,4);
        nI = dep(:,5);
        T = dep(:,6);
    end

    % function that performs the analysis
	function perform_the_analysis()
        % set a range of H2O inlet flow rates
        nH2O_range = linspace(3*nCO_in,6*nCO_in,100);

        % allocate storage for the corresponding volumes and temperatures
        V_range = nan(100,1);
        T_range = nan(100,1);

        % calculate the volume and outlet temperature for each H2O flow rate
        for i=1:100
            nH2O_in = nH2O_range(i);

            % solve the reactor design equations
            [V, ~, ~, ~, ~, ~, T] = profiles();

            V_range(i) = V(end);
            T_range(i) = T(end) - 273.15;
        end

        % identify the minimum volume and the corresponding temperature
        [Vmin, i_min] = min(V_range);
        Tmin = T_range(i_min);
        nH2O_opt = nH2O_range(i_min);
        
        % tabulate the results
        item = ["Optimum H20 Flow";"Volume";"°C"];
        value = [nH2O_opt; Vmin; Tmin];
        units = ["mol h^-1^";"cm^3^";"°C"];
        results_table = table(item,value,units);

        % display the results
        disp(' ')
        disp(results_table)

        % save the results
        writetable(results_table,'results.csv');

        % display and save the graphs
        figure;
        plot(nH2O_range,V_range,'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('H_2O Feed Rate (mol h_-_1)','FontSize', 14)
        ylabel('Volume (cm^3)','FontSize', 14)
        saveas(gcf,"V_vs_H2O_in.png")

        figure;
        plot(nH2O_range,T_range,'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('H_2O Feed Rate (mol h_-_1)','FontSize', 14)
        ylabel('Outlet Temperature (°C)','FontSize', 14)
        saveas(gcf,"T_vs_H2O_in.png")  
    end

    % perform the analysis
    perform_the_analysis()
end