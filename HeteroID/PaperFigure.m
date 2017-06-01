% Result in the paper
clear all;
close all;
clc;

parameter.plot_onoff = 'off'; % turn off plot to save computation
t_paperfigure_start = tic;
replicate_index = 0; % index for number of experiments
T_index = 0; % index for the length of single experiment

for which_gene = 1:7
    % Specify which gene you are interested in. The maximum should be less than the dimension of the system
    for montecarlo_index = 1:50 % montecarlo number for a single setup
        for replicate = 1:1:10
            replicate_index = round(replicate/1);
            for T = 10:10:100
                
                T_index = round(T/10);
                run('RUN_MULTIPLE.m')
                rnmse_w_iter_1( replicate_index,T_index,montecarlo_index, which_gene)=...
                    result.rnmse_w_iter_1;
                rnmse_w_iter_end(replicate_index,T_index,montecarlo_index, which_gene)=...
                    result.rnmse_w_iter_end;
            end
        end
        save('rnmse_w_iter_1.mat', 'rnmse_w_iter_1')
        save('rnmse_w_iter_end.mat', 'rnmse_w_iter_end')
    end
end

t_paperfigure_end = toc(t_paperfigure_start);

%%

% Plot the figure in the paper
clear all;
close all;

load rnmse_w_iter_1.mat
load rnmse_w_iter_end.mat

replicate_index = 0; % index for number of experiments
T_index = 0; % index for the length of single experiment

for which_gene = 1:7
    for replicate = 1:1:10
        replicate_index = round(replicate/1);
        for T = 10:10:100
            T_index = round(T/10);
            rnmse_w_iter_1_mean(replicate_index, T_index, which_gene) = ...
                mean(rnmse_w_iter_1(replicate_index, T_index, :, which_gene));
            
            rnmse_w_iter_end_mean(replicate_index, T_index, which_gene) = ...
                mean(rnmse_w_iter_end(replicate_index, T_index, :, which_gene));
        end
    end
end


figure;
surf(1:10,10:10:100,rnmse_w_iter_1_mean(:,:, which_gene))
xlabel('Number of experiemnt: C')
ylabel('Length of single time series')
zlabel('RNMSE')
figure;
surf(1:10,10:10:100,rnmse_w_iter_end_mean(:,:, which_gene))
xlabel('Number of experiemnt: C')
ylabel('Length of single time series')
zlabel('RNMSE')





