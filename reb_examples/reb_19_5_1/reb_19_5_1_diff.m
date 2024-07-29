function reb_19_5_1_diff()
%reb_16_5_1 Reaction Engineering Basics Example 19.5.1

    % given and known 
    V = 1.0; % L
    R = 8.314E-3; % kJ /mol /K

    % function that performs the calculations
    function perform_the_calculations()
        % read the experimental data
        data_file = '../reb_19_5_1_data.csv';
        data_table = readtable(data_file, 'VariableNamingRule'...
            , 'preserve');
        data = table2array(data_table(:,2:end));

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
            idx = table2array(data_table(:,2) == block_temperatures(iBlock));
            data_block = data(idx,:);

            % extract the data as arrays
            CA0 = data_block(:,2);
            tf = data_block(:,3);
            CAf = data_block(:,4);
            nData = length(tf);

            % calculate x
            x = CAf*V;

            % allocate storage for y
            y = nan(length(tf),1);

            % calculate y
            for i = 1:length(tf)
                if tf(i) == 5.0
                    x_minus = CA0(i)*V;
                    t_minus = 0;
                else
                    x_minus = x(i-1);
                    t_minus = tf(i-1);
                end
                y(i) = (x_minus - x(i))/(tf(i) - t_minus);
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
            fname = strcat('reb_19_5_1_model_'...
                , int2str(block_temperatures(iBlock)),'.png');
            saveas(gcf,fname)
        end

        % tabulate, show, and save the results
        results_file ="reb_19_5_1_diff_params.csv";
        results_table = table(block_temperatures, k_fit, kll_fit...
            , kul_fit, rsq_fit);
        results_table.Properties.VariableNames = ["T", "k", "k_ll"...
            , "k_ul", "R_sq"];
        disp(results_table)
        writetable(results_table,results_file);

        % fit the Arrhenius expression to the k vs. T data
        T = block_temperatures + 273.15;
        [k0, k0_ci, E, E_ci, r_squared] = Arrhenius_parameters(k_fit...
            ,T,R);

        % tabulate, show, and save the results
        item = ["k0"; "k0_lower_limit"; "k0_upper_limit"; "E"; 
            "E_lower_limit"; "E_upper_limit"; "R_squared"];
        value = [k0; k0_ci(1); k0_ci(2); E; E_ci(1); E_ci(2);
             r_squared];
        units = ["min^-1^"; "min^-1^"; "min^-1^"; "kJ mol^-1^";
            "kJ mol^-1^"; "kJ mol^-1^"; ""];
        Arrhenius_table = table(item, value, units);
        disp(Arrhenius_table)
        writetable(Arrhenius_table,'reb_19_5_1_Arrhenius_diff.csv')

        % create, show and save an Arrhenius plot
        k_pred = k0*exp(-E/R./T);
        x = 1./T;
        figure
        semilogy(x,k_fit,'ok',x,k_pred,'r','MarkerSize',10 ...
            ,'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('1/T (/K)','FontSize', 14)
        ylabel('k (/min)','FontSize', 14)
        saveas(gcf,'reb_19_5_1_Arrhenius_diff.png')
    end

    % perform the calculations
    perform_the_calculations()
end