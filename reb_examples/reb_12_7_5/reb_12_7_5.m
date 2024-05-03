function reb_12_7_5
%REB_12_7_5 Reaction Engineering Basics Example 12.7.5
    % constants available to all functions
	% given
    V = 10; % gal
    Ve = 1.25; % gal
    CB0 = 3; % mol/gal
    T0 = 20 + 273.15; % K
    VFR = 0.5; % gal/min
    CAin = 7.5; % mol/gal
    CBin = 3; % mol/gal
    Tin = 50 + 273.15; % K
    Tein = 20 + 273.15; % K
    Te0 = 20 + 273.15; % K
    mein = 250; % g/min
    U = 190; % cal/ft^2/min/K
    A = 4; % ft^2
    rhoe = 1; % g/cc
    Cpe = 1; % cal/g/K
    Cp = 1600; % cal/gal/K
    tf = 30; % min
    k01 = 4.3E8; % gal/mol/min
    E1 = 14200; % cal/mol
    dH1 = -11000; % cal/mol
    k02 = 2.7E8; % gal/mol/min
    E2 = 16100; % cal/mol
    dH2 = -11600; % cal/mol
    k03 = 3.9E8; % gal/mol/min
    E3 = 14800; % cal/mol
    dH3 = -12100; % cal/mol
	% known
    R = 1.987; % cal/mol/K
	% calculated
    nAin = CAin*VFR;
    nBin = CBin*VFR;

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract the dependent variables for this integration step
        nA = dep(1);
        nB = dep(2);
        nW = dep(3);
        nX = dep(4);
        nY = dep(5);
        nZ = dep(6);
        T = dep(7);
        Te = dep(8);

        % calculate rate and rate of heat transfer
        r1 = k01*exp(-E1/R/T)*nA/VFR*nB/VFR;
        r2 = k02*exp(-E2/R/T)*nA/VFR*nB/VFR;
        r3 = k03*exp(-E3/R/T)*nA/VFR*nB/VFR;
        Q = U*A*(Te-T);

        % evaluate the derivatives
        dnAdt = VFR/V*(nAin - nA - V*(r1 + r2 + r3));
        dnBdt = VFR/V*(nBin - nB - V*(r1 + r2 + r3));
        dnWdt = VFR/V*(-nW + r1*V);
        dnXdt = VFR/V*(-nX + r2*V);
        dnYdt = VFR/V*(-nY + r2*V);
        dnZdt = VFR/V*(-nZ + V*(r1 + r2 + r3));
        dTdt = (Q - VFR*Cp*(T-Tin) - V*(r1*dH1 + r2*dH2 + r3*dH3))/V/Cp;
        dTedt = (-Q-mein*Cpe*(Te-Tein))/(rhoe*Ve*Cpe);

        % return the derivatives
        derivs = [dnAdt; dnBdt; dnWdt; dnXdt; dnYdt; dnZdt; dTdt; 
            dTedt];
    end

    % reactor model
    function [t, nA, nB, nW, nX, nY, nZ, T, Tex] = profiles()
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [0.0; VFR*CB0; 0.0; 0.0; 0.0; 0.0; T0; Te0];

        % define the stopping criterion
        f_var = 0;
        f_val = tf;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [t, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract and return the dependent variable profiles
        nA = dep(:,1);
        nB = dep(:,2);
        nW = dep(:,3);
        nX = dep(:,4);
        nY = dep(:,5);
        nZ = dep(:,6);
        T = dep(:,7);
        Tex = dep(:,8);
    end

    % function that performs the analysis
	function perform_the_analysis()

        % solve the reactor design equations
        [t, ~, nB, ~, ~, ~, ~, T, Tex] = profiles();

        % calculate the other quantities of interest
        CB = nB/VFR;
        T = T - 273.15;
        Tex = Tex - 273.15;

        % display and save the graphs
        figure % reactor temperature vs. time
        plot(t,T,'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Time (min)', 'FontSize', 14)
        ylabel('Reactor Temperature (°C)', 'FontSize', 14)
        saveas(gcf,"T_vs_t.png")
        
        figure % shell temperature vs. time
        plot(t,Tex,'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Time (min)', 'FontSize', 14)
        ylabel('Shell Temperature (°C)', 'FontSize', 14)
        saveas(gcf,"Te_vs_t.png")
    
        figure % CB vs. time
        plot(t,CB,'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Time (min)', 'FontSize', 14)
        ylabel('Concentration of B (mol/gal)', 'FontSize', 14)
        saveas(gcf,"CB_vs_t.png")

    end

    % perform the analysis
    perform_the_analysis()
end