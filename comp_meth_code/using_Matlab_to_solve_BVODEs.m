function using_Matlab_to_solve_BVODEs()
% Calculations for Using Matlab to Solve BVODEs
    % global constants
    Dax = 4.8E-4; % dm^2 /min
    k = 0.72; % /min
    K = 1.0;
    L = 1.25; % dm
    D = 0.07; % dm
    Vdot = 0.0023; % dm^3 /min
    nA_feed = 0.0023; % mol /min
    nZ_feed = 0.0; % mol /min

    % derivatives function
    function ddz = derivatives(~, dep)
        % extract the dependent variables
        nA = dep(1);
        nZ = dep(2);
        wA = dep(3);
        wZ = dep(4);

        % evaluate the derivatives
        dnAdz = wA;
        dnZdz = wZ;
        dwAdz = 1/Dax*4*Vdot/pi/D^2*wA + k/Dax*nA - k/Dax/K*nZ;
        dwZdz = 1/Dax*4*Vdot/pi/D^2*wZ - k/Dax*nA + k/Dax/K*nZ;
   
        % combine the derivatives in a vector and return
        ddz = [dnAdz; dnZdz; dwAdz; dwZdz];
    end

    % boundary conditions residuals function
    function epsilon = bc_residuals(dep_lb, dep_ub)
        % extract the boundary values needed to evaluate the residuals
        nAlb = dep_lb(1);
        nZlb = dep_lb(2);
        wAlb = dep_lb(3);
        wZlb = dep_lb(4);
        wAub = dep_ub(3);
        wZub = dep_ub(4);

        % evaluate the residuals
        epsilon_1 = nAlb - nA_feed - Dax*pi*D^2/4/Vdot*wAlb;
        epsilon_2 = nZlb - nZ_feed - Dax*pi*D^2/4/Vdot*wZlb;
        epsilon_3 = wAub;
        epsilon_4 = wZub;

        % combine the residuals as a vector and return
        epsilon = [epsilon_1; epsilon_2; epsilon_3; epsilon_4];
    end

    % axial dispersion reactor function
    function [z, nA, nZ, wA, wZ] = reactor_variables()
        % set the initial mesh with 20 mesh points
        ind = linspace(0, L, 20);

        % set the guess
        depGuess = [0.0, 0.0, 0.0, 0.0];
        solinit = bvpinit(ind,depGuess);

        % solve the BVODEs
        soln = bvp4c(@derivatives, @bc_residuals, solinit);

        % extract and return the profiles
        z = soln.x;
        nA = soln.y(1,:);
        nZ = soln.y(2,:);
        wA = soln.y(3,:);
        wZ = soln.y(4,:);
    end

    % deliverables function
    function deliverables()
        % get the reactor profiles
        [z, nA, nZ, ~, ~] = reactor_variables();

        % report the initial and final mesh sizes
        disp(" ")
        disp(['Size of initial mesh: 20, size of final mesh: '...
            num2str(length(z),3)])
        disp(" ")

        % plot the results
        figure;
        plot(z,nA,'k',z,nZ,'b','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('z (dm)','FontSize', 14)
        ylabel('Molar Flow Rate (mol/min)','FontSize', 14)
        legend({'A','Z'},'Location','northeast','FontSize',14)
        saveas(gcf,"results.pdf")
    end

    % perform the analysis
    deliverables();
end