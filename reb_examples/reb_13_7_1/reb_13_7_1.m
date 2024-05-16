function reb_13_7_1
    %REB_13_7_1 Calculations for Example 13.7.1 of Reaction Engineering Basics
        % constants available to all functions
        % given
        dH_1 = 44800; % J /mol
        k0_1 = 7.22E6*60; % mol /atm^2 /cm^3 /min
        E_1 = 84100; % J /mol
        L = 10.0*12*2.54; % cm
        D = 2.54; % cm
        T_ex = 200 + 273.15; % K
        U = 7.48E4/60/12.^2/2.54^2; % J /min /cm^2 /K
        yA_in = 0.6;
        yB_in = 0.4;
        Vdot_in = 282E3; % cm^3 /min
        P = 2.5; % atm
        T_in = 175 + 273.15; % K
        CpA = 18.0*4.184; %J /mol /K
        CpB = 12.25*4.184; % J /mol /K
        CpZ = 21.2*4.184; % J /mol /K
        % known
        Re = 8.314; % J /mol /K
        Rw = 82.06; % cm^3 atm /mol /K
        % calculated
        nA_in = yA_in*P*Vdot_in/Rw/T_in;
        nB_in = yB_in*P*Vdot_in/Rw/T_in;
    
        % derivatives function
        function derivs = derivatives(~, dep)
            % extract the dependent variables for this integration step
            nA = dep(1);
            nB = dep(2);
            nZ = dep(3);
            T = dep(4);

            % calculate the rate
            k_1 = k0_1*exp(-E_1/Re/T);
            PA = nA/(nA + nB + nZ)*P;
            PB = nB/(nA + nB + nZ)*P;
            r_1 = k_1*PA*PB;

            % evaluate the derivatives
            dnAdz = -pi()*D^2/4*r_1;
            dnBdz = -pi()*D^2/4*r_1;
            dnZdz = pi()*D^2/4*r_1;
            dTdz = (pi()*D*U*(T_ex - T) - pi()*D^2/4*r_1*dH_1)...
                /(nA*CpA + nB*CpB + nZ*CpZ);
            % return the derivatives
            derivs = [dnAdz; dnBdz; dnZdz; dTdz];
        end
    
        % reactor model
        function [z, nA, nB, nZ, T] = profiles()
            % set the initial values
            ind_0 = 0.0;
            dep_0 = [nA_in; nB_in; 0.0; T_in];
    
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
            nB = dep(:,2);
            nZ = dep(:,3);
            T = dep(:,4);
        end
    
        % function that performs the analysis
        function perform_the_analysis()
    
            % solve the reactor design equations
            [~, ~, nB, ~, T] = profiles();
    
            % calculate the other quantities of interest
            T_out = T(end);
            f_B = 100*(nB_in - nB(end))/nB_in;

            % tabulate the results
            item = ["T out";"conversion"];
            value = [T_out - 273.15; f_B];
            units = ["°C";"%"];
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