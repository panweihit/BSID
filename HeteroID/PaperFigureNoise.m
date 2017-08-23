% Result in the paper
clear all;
close all;
%clc;

replicate = 10; T = 103; 

parameter.plot_onoff = 'off'; % turn off plot to save computation
t_paperfigure_start = tic;
measurement_index = 0; % index for number of experiments
process_index = 0; % index for the length of single experiment

gene_range = 1:1:8;
t_measurement_range = 0.01;
measurement_range =  0.01:t_measurement_range:0.1;
t_process_range = 0.01;
process_range = 0.01:t_process_range:0.1;

for which_var = gene_range
    % Specify which gene you are interested in. The maximum should be less than the dimension of the system
    for montecarlo_index = 1:50 % montecarlo number for a single setup
        for measurementnoisestd = measurement_range
            measurement_index = round(measurementnoisestd/t_measurement_range);
            for processnoisestd = process_range
                process_index = round(processnoisestd / t_process_range);
                run('RUN_MULTIPLE.m')
                rnmse_w_iter_1_noise( measurement_index,process_index,montecarlo_index, which_var)=...
                    result.rnmse_w_iter_1;
                rnmse_w_iter_end_noise(measurement_index,process_index,montecarlo_index, which_var)=...
                    result.rnmse_w_iter_end;
                fprintf('montecarlo_index is %d. \n',montecarlo_index);
                fprintf('measure_index is %d. \n',measurement_index);
                fprintf('process_index is %d. \n',process_index);
                fprintf('rnmse_w_iter_1_noise is %f. \n',result.rnmse_w_iter_1);
                fprintf('rnmse_w_iter_end_noise is %f. \n \n',result.rnmse_w_iter_end);
            end
        end
        save('rnmse_w_iter_1_noise.mat', 'rnmse_w_iter_1_noise')
        save('rnmse_w_iter_end_noise.mat', 'rnmse_w_iter_end_noise')
    end
end

t_paperfigure_end = toc(t_paperfigure_start);

%%

% Plot the figure in the paper
close all;

load rnmse_w_iter_1_noise.mat
load rnmse_w_iter_end_noise.mat

measurement_index = 0; % index for number of experiments
process_index = 0; % index for the length of single experiment
replicate = 10; T = 103; 

parameter.plot_onoff = 'off'; % turn off plot to save computation
t_paperfigure_start = tic;
measurement_index = 0; % index for number of experiments
process_index = 0; % index for the length of single experiment

gene_range = 1:1:8;
t_measurement_range = 0.01;
measurement_range =  0.01:t_measurement_range:0.1;
t_process_range = 0.01;
process_range = 0.01:t_process_range:0.1;


for measurementnoisestd = measurement_range
    measurement_index = round(measurementnoisestd/t_measurement_range);
    for processnoisestd = process_range
        process_index = round(processnoisestd/t_process_range);
        rnmse_w_iter_1_noise_mean(measurement_index, process_index) = ...
            mean(mean(rnmse_w_iter_1_noise(measurement_index, process_index, :, gene_range)));
        
        rnmse_w_iter_end_noise_mean(measurement_index, process_index) = ...
            mean(mean(rnmse_w_iter_end_noise(measurement_index, process_index, :, gene_range)));
    end
end



figure;
surf(process_range,measurement_range,rnmse_w_iter_1_noise_mean)
xlabel('Process Noise Standard Deviation')
ylabel('Measurement Noise Standard Deviation')
zlabel('RNMSE')
figure;
surf(process_range, measurement_range,rnmse_w_iter_end_noise_mean)
xlabel('Process Noise Standard Deviation')
ylabel('Measurement Noise Standard Deviation')
zlabel('RNMSE')





