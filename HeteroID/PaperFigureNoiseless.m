% Result in the paper
clear all;
close all;

processnoisestd = 0.0; measurementnoisestd = 0.0;

parameter.plot_onoff = 'off'; % turn off plot to save computation
t_paperfigure_start = tic;
replicate_index = 0; % index for number of experiments
T_index = 0; % index for the length of single experiment

gene_range = 1:8;
t_replicate_range = 1;
replicate_range =  1:t_replicate_range:10;
t_T_range = 10;
T_range = 13:t_T_range:113;

for which_var = gene_range
    % Specify which gene you are interested in. The maximum should be less than the dimension of the system
    for montecarlo_index = 1:50 % montecarlo number for a single setup
        for replicate = replicate_range
            replicate_index = round(replicate/t_replicate_range);
            for T = T_range
                T_index = round(T / t_T_range);
                run('RUN_MULTIPLE.m')
                rnmse_w_iter_1( replicate_index,T_index,montecarlo_index, which_var)=...
                    result.rnmse_w_iter_1;
                rnmse_w_iter_end(replicate_index,T_index,montecarlo_index, which_var)=...
                    result.rnmse_w_iter_end;
                fprintf('montecarlo_index is %d. \n',montecarlo_index);
                fprintf('replicate_index is %d. \n',replicate_index);
                fprintf('T_index is %d. \n',T_index);
                fprintf('rnmse_w_iter_1 is %f. \n',result.rnmse_w_iter_1);
                fprintf('rnmse_w_iter_end is %f. \n \n',result.rnmse_w_iter_end);
            end
        end
        save('rnmse_w_iter_1_noiseless.mat', 'rnmse_w_iter_1')
        save('rnmse_w_iter_end_noiseless.mat', 'rnmse_w_iter_end')
    end
end

t_paperfigure_end = toc(t_paperfigure_start);

%%

%% Plot the figure in the paper
clear all;
close all;

load rnmse_w_iter_1_noiseless.mat
load rnmse_w_iter_end_noiseless.mat

replicate_index = 0; % index for number of experiments
T_index = 0; % index for the length of single experiment
processnoisestd = 0.0; measurementnoisestd = 0.0;

parameter.plot_onoff = 'off'; % turn off plot to save computation
t_paperfigure_start = tic;
replicate_index = 0; % index for number of experiments
T_index = 0; % index for the length of single experiment

gene_range = 1:8;
t_replicate_range = 1;
replicate_range =  1:t_replicate_range:10;
t_T_range = 10;
T_range = 13:t_T_range:113;

for replicate = replicate_range
    replicate_index = round(replicate/t_replicate_range);
    for T = T_range
        T_index = round(T / t_T_range);
        rnmse_w_iter_1_mean(replicate_index, T_index) = ...
            mean(mean(rnmse_w_iter_1(replicate_index, T_index, :,gene_range)));
        
        rnmse_w_iter_end_mean(replicate_index, T_index) = ...
            mean(mean(rnmse_w_iter_end(replicate_index, T_index, :,gene_range)));
    end
end



figure;
surf(T_range(:),replicate_range(:),rnmse_w_iter_1_mean)
ylabel('Number of experiemnt: C')
xlabel('Length of single time series')
zlabel('RNMSE')
figure;
surf(T_range(:),replicate_range(:),rnmse_w_iter_end_mean)
ylabel('Number of experiemnt: C')
xlabel('Length of single time series')
zlabel('RNMSE')





