function reb_13_7_5
%REB_13_7_5 Calculations for Example 13.7.5 of Reaction Engineering Basics
    % constants available to all functions
	% given
    D = 10.0; % cm
    L = 500.0; % cm
    Vdot_in = 75E3; % mv^3 /min
    T_before = 25 + 273.15; % K
    CA_in = 1.0E-3; % mol /cm^3
    CB_in = 1.2E-3; % mol /cm^3
    T_in = 30 + 273.15; % K
    k0_1 = 8.72E8; % cm^3 /mol /min
    E_1 = 7200; % cal /mol
    dH_1 = -10700; % cal /mol
    Cp = 1.0; % cal /g /K
    rho = 1.0; % g \cm^3
    t = [0.131; 0.262; 0.393; 0.524]; % min
	% known
    R = 1.987; % cal /mol /K
	% calculated
    nA_in = Vdot_in*CA_in;
    nB_in = Vdot_in*CB_in;
    %V = pi()*D^2/4*L;

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract the dependent variables for this integration step
        nA = dep(1);
        nB = dep(2);
        T = dep(5);

        % calculate [quantities in the derivatives]
        k_1 = k0_1*exp(-E_1/R/T);
        CA = nA/Vdot_in;
        CB = nB/Vdot_in;
        r_1 = k_1*CA*CB;

        % evaluate the derivatives
        dnAdz = - pi()*D^2/4*r_1;
        dnBdz = - pi()*D^2/4*r_1;
        dnYdz = pi()*D^2/4*r_1;
        dnZdz = pi()*D^2/4*r_1;
        dTdz = - pi()*D^2/4*r_1*dH_1/(Vdot_in*rho*Cp);

        % return the derivatives
        derivs = [dnAdz; dnBdz; dnYdz; dnZdz; dTdz];
    end

    % reactor model
    function [z, nA, nB, nY, nZ, T] = profiles(z_front)
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [nA_in; nB_in; 0.0; 0.0; T_in];

        % define the stopping criterion
        f_var = 0;
        f_val = z_front;
        
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
        nB = dep(:,2);
        nY = dep(:,3);
        nZ = dep(:,4);
        T = dep(:,5);
    end

    % function that performs the analysis
	function perform_the_analysis()
        % analyze each of the four specified elapsed times
        for i = 1:4
            % calculate the position of the front
            V_front = t(i)*Vdot_in;
            z_front = V_front/(pi()*D^2/4);

            % calculate the profiles up to the front
            [z, nA, ~, ~, ~, T] = profiles(z_front);

            % add the profile beyond the front
            if z(end) < L
                z = [z; z_front; L];
                nA = [nA; 0.0; 0.0];
                T = [T; T_before; T_before];
            end

            % generate, show, and save the graph
            figure;
            yyaxis left
            plot(z,nA,'LineWidth',2)
            title(['Profiles after ',num2str(t(i)),' min'],'FontSize',14)
            set(gca, 'FontSize', 14);
            xlabel('Axial Distance from Inlet (cm)','FontSize', 14)
            ylabel('Flow Rate of A (mol min^-^1)','FontSize', 14)
            yyaxis right
            plot(z,T-273.15,'LineWidth',2)
            ylabel('Temperature (°C)','FontSize', 14)
            filePath = ['t_',num2str(1000*t(i)),'_min_profiles.png'];
            saveas(gcf,filePath)
        end
    end

    % perform the analysis
    perform_the_analysis()
end