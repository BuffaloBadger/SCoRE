function reb_19_5_2_linear()
%reb_16_5_2_linear Reaction Engineering Basics Example 19.5.2

    % given and known
    V = 500.0; % cc
    Re = 1.987E-3; % kcal/mol/K
    Rpv = 82.06; % cc atm/mol/K

    % function that performs the calculations
    function perform_the_calculations()
        % read the experimental data
        data_file = '../reb_19_5_2_data.csv';
        data_table = readtable(data_file, 'VariableNamingRule'...
            , 'preserve');
        data = table2array(data_table(:,:));

        % get the temperatures of the data blocks
        block_temperatures = unique(data(:,1));
        n_blocks = size(block_temperatures,1);

        % create vectors to store the fitting results
        k_fit = nan(n_blocks,1);
        kll_fit = nan(n_blocks,1);
        kul_fit = nan(n_blocks,1);
        rsq_fit = nan(n_blocks,1);

        % process the data blocks
        for iBlock = 1:n_blocks
            % create the block
            idx = table2array(data_table(:,1) == block_temperatures(iBlock));
            data_block = data(idx,:);

            % extract the data as arrays
            T = data_block(:,1);
            T_K = T + 273.15;
            PA0 = data_block(:,2);
            PB0 = data_block(:,3);
            tf = data_block(:,4);
            fA = data_block(:,5);
            nData = length(tf);

            % allocate storage for x and y
            x = nan(nData,1);
            y = nan(nData,1);

            % calculate x and y
            for i = 1:nData 
                nA0 = PA0(i)*V/Rpv/T_K(i);
                nB0 = PB0(i)*V/Rpv/T_K(i);
                nA = nA0*(1-fA(i));
                nB = nB0 - nA0 + nA;
                x(i) = -tf(i)*(Rpv*T_K(i))^2/V;
                if nA0 == nB0
                    y(i) = 1/nA0 - 1/nA;
                else
                    y(i) = 1/(nA0-nB0)*log(nA0*nB/nB0/nA);
                end
            end

            % fit y = m*x to the data
            beta = inv(x.'*x)*x.'*y;

            % calculate the model-predicted y
            y_model = x*beta;

            % calculate the 95% confidence interval
            nParams = length(beta);
            eps = y - x*beta;
            varEps = 1/(nData - nParams)*(eps.'*eps);
            covBeta = varEps*inv(x.'*x);
            tCrit = tinv(0.975,(nData - nParams));
            betaCI = sqrt(covBeta)*tCrit;

            % calculate R-squared
            ybar = sum(y)/nData;
            r_squared = (sum((y-ybar).^2) - sum((y-y_model).^2))/sum((y-ybar).^2);

            % save the fitting results 
            k_fit(iBlock) = beta;
            kll_fit(iBlock) = beta - betaCI;
            kul_fit(iBlock) = beta + betaCI;
            rsq_fit(iBlock) = r_squared;

            % create, show, and save a model plot
            figure
            hold on
            plot(x,y_model,'r','LineWidth',2)
            plot(x, y,'ok','MarkerSize',10,'LineWidth',2)
            hold off
            set(gca, 'FontSize', 14);
            xlabel('x','FontSize', 14)
            ylabel('y','FontSize', 14)
            fname = strcat('reb_19_5_2_model_'...
                , int2str(block_temperatures(iBlock)),'.png');
            saveas(gcf,fname)
        end

        % tabulate, show, and save the results
        results_file ="reb_19_5_2_linear_results.csv";
        results_table = table(block_temperatures, k_fit, kll_fit...
            , kul_fit, rsq_fit);
        results_table.Properties.VariableNames = ["T", "k", "k_ll"...
            , "k_ul", "R_sq"];
        disp(results_table)
        writetable(results_table,results_file);

        % fit the Arrhenius expression to the k vs. T data
        T = block_temperatures + 273.15;
        [k0, k0_ci, E, E_ci, r_squared] = Arrhenius_parameters(k_fit...
            ,T,Re);

        % tabulate, show, and save the results
        item = ["k0"; "k0_lower_limit"; "k0_upper_limit"; "E"; 
            "E_lower_limit"; "E_upper_limit"; "R_squared"];
        value = [k0; k0_ci(1); k0_ci(2); E; E_ci(1); E_ci(2);
             r_squared];
        units = ["min^-1^"; "min^-1^"; "min^-1^"; "kJ mol^-1^";
            "kJ mol^-1^"; "kJ mol^-1^"; ""];
        Arrhenius_table = table(item, value, units);
        disp(Arrhenius_table)
        writetable(Arrhenius_table,'reb_19_5_2_Arrhenius_params.csv')

        % create, show and save an Arrhenius plot
        k_pred = k0*exp(-E/Re./T);
        x = 1./T;
        figure
        semilogy(x,k_fit,'ok',x,k_pred,'r','MarkerSize',10 ...
            ,'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('1/T (/K)','FontSize', 14)
        ylabel('k (/min)','FontSize', 14)
        saveas(gcf,'reb_19_5_2_Arrhenius_plot.png')
    end

    % perform the calculations
    perform_the_calculations()
end