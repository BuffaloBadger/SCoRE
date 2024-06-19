function reb_15_4_2
%REB_15_4_2 Calculations for Example 15.4.2 of Reaction Engineering Basics
    % constants available to all functions
	% given
    dH = -9120; % cal /mol
    K0 = 0.132;
    CpCO = 29.3/4.184; % cal /mol /K
    CpH2O = 34.3/4.184; % cal /mol /K
    CpCO2 = 41.3/4.184; % cal /mol /K
    CpH2 = 29.1/4.184; % cal /mol /K
    CpI = 40.5/4.184; % cal /mol /K
    nDotCO_in = 1.0; % mol /h
    nDotCO2_in = 0.359; % mol /h
    nDotH2_in = 4.44; % mol /h
    nDotI_in = 0.18; % mol/h
    nDotH2O_in = 9.32; % mol /h
    P = 26.0; % atm
    T_in = 445 + 273.15; % K
    V_1 = 685; % cm^3
    k0_1 = 3.54E-2; % mol /cm^3 /min /atm^2
    E_1 = 9740; % cal /mol
    Tex_in = 20 + 273.15; % K
    mDotEx = 1100; % g /h
    Tex_out = 50 + 273.15; % K
    CpEx = 1.0; % cal /g /K
    V_2 = 3950; % cm^3
    k0_2 = 1.77E-3; % mol /cm^3 /min /atm^2
    E_2 = 3690; % cal /mol
	% known
    R = 1.987; % cal /mol /K

    % make the current reactor available to all functions
    reactor_number = nan;

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract the dependent variables for this integration step
        nDotCO = dep(1);
        nDotH2O = dep(2);
        nDotCO2 = dep(3);
        nDotH2 = dep(4);
        nDotI = dep(5);
        T = dep(6);

        % calculate the rate
        nDotTotal = nDotCO + nDotH2O + nDotCO2 + nDotH2 + nDotI;
        PCO = P*nDotCO/nDotTotal;
        PH2O = P*nDotH2O/nDotTotal;
        PCO2 = P*nDotCO2/nDotTotal;
        PH2 = P*nDotH2/nDotTotal;
        K = K0*exp(-dH/R/T);
        if reactor_number == 1
            r = k0_1*exp(-E_1/R/T)*(PCO*PH2O - PCO2*PH2/K);
        else
            r = k0_2*exp(-E_2/R/T)*(PCO*PH2O - PCO2*PH2/K);
        end

        % evaluate the derivatives
        dCOdV = -r;
        dH2OdV = -r;
        dCO2dV = r;
        dH2dV = r;
        dIdV = 0;
        dTdV = -r*dH/(nDotCO*CpCO + nDotH2O*CpH2O + nDotCO2*CpCO2 ...
            + nDotH2*CpH2 + nDotI*CpI);

        % return the derivatives
        derivs = [dCOdV; dH2OdV; dCO2dV; dH2dV; dIdV; dTdV];
    end

    % reactor model
    function [V, nDotCO, nDotH2O, nDotCO2, nDotH2, nDotI, T] ...
            = profiles(dep_0, f_val)
        % set the initial values
        ind_0 = 0.0;

        % define the stopping criterion
        f_var = 0;
        
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
        nDotCO = dep(:,1);
        nDotH2O = dep(:,2);
        nDotCO2 = dep(:,3);
        nDotH2 = dep(:,4);
        nDotI = dep(:,5);
        T = dep(:,6);
    end

    % function that performs the analysis
	function perform_the_analysis()

        % reactor 1
        reactor_number = 1;
        dep0 = [nDotCO_in; nDotH2O_in; nDotCO2_in; nDotH2_in; nDotI_in
            T_in];
        f_val = V_1;

        % solve the reactor design equations
        [V, nDotCO, nDotH2O, nDotCO2, nDotH2, nDotI, T] ...
            = profiles(dep0, f_val);
        T_1 = T(end) - 273.15;
        fCO_1 = 100*(nDotCO_in - nDotCO(end))/nDotCO_in;

        % heat exchanger
        T_in_2 = T(end) - mDotEx*CpEx*(Tex_out - Tex_in)...
            /(nDotCO(end)*CpCO + nDotH2O(end)*CpH2O ...
            + nDotCO2(end)*CpCO2 + nDotH2(end)*CpH2 + nDotI(end)*CpI);

        % reactor 2
        reactor_number = 2;
        dep0 = [nDotCO(end); nDotH2O(end); nDotCO2(end); nDotH2(end) 
            nDotI(end); T_in_2];
        f_val = V_2;

        % solve the reactor design equations
        [V, nDotCO, ~, ~, ~, ~, T] = profiles(dep0, f_val);

        % calculate the other quantities of interest
        T_in_2 = T_in_2 - 273.15;
        T_2 = T(end) - 273.15;
        fCO_2 = 100*(nDotCO_in - nDotCO(end))/nDotCO_in;

        % tabulate the results
        item = ["Reactor 1 Outlet T"; "Reactor 1 Conversion"
            "Reactor 2 Inlet T"; "Reactor 2 Outlet T"
            "Reactor 2 Conversion"];
        value = [T_1; fCO_1; T_in_2; T_2; fCO_2];
        units = ["°C"; "%"; "°C"; "°C"; "%"];
        results_table = table(item,value,units);

        % display the results
        disp(results_table)

        % save the results
        writetable(results_table,'results.csv');
    end

    % perform the analysis
    perform_the_analysis()
end