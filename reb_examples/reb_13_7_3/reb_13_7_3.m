function reb_13_7_3
%REB_13_7_3 Calculations for Example 13.7.3 of Reaction Engineering Basics
    % constants available to all functions
	% given
    yA_in = 0.75;
    yI_in = 0.25;
    Vdot_in = 100; % cm^3 /s
    mDot_in = 0.44; % g /s
    P_in = 3.0; % atm
    T_in = 400 + 273.15; % K
    D = 2.5; % cm
    L = 800.0; % cm
    T_ex = 375 + 273.15; % K
    U = 187E3/3600*1e-4; % J /s /cm^2 /K
    Dp = 0.25; % cm
    phi = 0.7;
    epsilon = 0.6;
    k_0f = 9E17; % mol/cm^3 /s /atm
    E_f = 285000.0; % J /mol
    k_0r = 4.09E-4; % mol/cm^3 /s /atm^4
    E_r = 85000.0; % J /mol
    dH_1 = 200000.0; % J /mol
    Cp_A = 11.7*4.184; % J /mol /K
    Cp_Y = 8.3*4.184; % J /mol /K
    Cp_Z = 4.3*4.184; % J /mol /K
    Cp_I = 5.8*4.184; % J /mol /K
    mu = 0.027E-2; % g /cm /s
	% known
    Re = 8.314; % J /mol /K
    Rw = 82.06; % cm^3 atm /mol /K
    P_conv = 9.872E-7; % atm cm^2 /dyne 
    % calculated
    nA_in = yA_in*P_in*Vdot_in/Rw/T_in;
    nI_in = yI_in*P_in*Vdot_in/Rw/T_in;
    G = mDot_in/pi()/D^2*4;

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract the dependent variables for this integration step
        nA = dep(1);
        nY = dep(2);
        nZ = dep(3);
        nI = dep(4);
        T = dep(5);
        P = dep(6);

        % calculate [quantities in the derivatives]
        kf = k_0f*exp(-E_f/Re/T);
        kr = k_0r*exp(-E_r/Re/T);
        ntot = nA + nY + nZ + nI;
        PA = nA/ntot*P;
        PY = nY/ntot*P;
        PZ = nZ/ntot*P;
        r_1 = kf*PA - kr*PY*PZ^3;

        % calculate the density
        Vdot = ntot*Rw*T/P;
        rho = mDot_in/Vdot;

        % evaluate the derivatives
        dnAdz = -pi()*D^2/4*r_1;
        dnYdz = pi()*D^2/4*r_1;
        dnZdz = 3*pi()*D^2/4*r_1;
        dnIdz = 0.0;
        dTdz = (pi()*D*U*(T_ex - T) - pi()*D^2/4*r_1*dH_1)...
            /(nA*Cp_A + nY*Cp_Y + nZ*Cp_Z + nI*Cp_I);
        dPdz = -P_conv*(1-epsilon)/epsilon^3*G^2/rho/phi/Dp ...
            *(150*(1-epsilon)*mu/phi/Dp/G + 1.75);

        % return the derivatives
        derivs = [dnAdz; dnYdz; dnZdz; dnIdz; dTdz; dPdz];
    end

    % reactor model
    function [z, nA, nY, nZ, nI, T, P] = profiles()
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [nA_in, 0.0, 0.0, nI_in, T_in, P_in];

        % define the stopping criterion
        f_var = 0;
        f_val = L;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [z, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract and return the dependent variable profiles
        nA = dep(:,1);
        nY = dep(:,2);
        nZ = dep(:,3);
        nI = dep(:,4);
        T = dep(:,5);
        P = dep(:,6);
    end

    % function that performs the analysis
	function perform_the_analysis()
        % solve the reactor design equations
        [~, nA, ~, ~, ~, T, P] = profiles();

        % calculate the quantities of interest
        T_out = T(end) - 273.15;
        P_out = P(end);
        fA = 100*(nA_in - nA(end))/nA_in;
    
        % tabulate the results
        item = ["T out"; "P out";"Conversion"];
        value = [T_out; P_out; fA];
        units = ["°C";"atm";"%"];
        results_table = table(item,value,units);

        % display the results
        disp(' ')
        disp(results_table)

        % save the results
        writetable(results_table,'results.csv');
    end

    % perform the analysis
    perform_the_analysis()
end