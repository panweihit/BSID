%% Attention
% If you want to run PaperFigure.m, comment the following line, otherwise uncomment to run standalone example 
% close all; clear all; which_var = 4; processnoisestd = 0.0; measurementnoisestd = 0.0; replicate = 10; T = 103; 


%% Simulated Example and Settings
run('config.m')
parameter.T = T; % length of a single time series
parameter.replicate = replicate;  %number of experiments/conditions/pertubations/replicates; 
                           %when replicate = 1, the algorithm is consistent with our TAC 2016 paper
which_var = which_var;  % Specify which gene you are interested in. 
                                 % The maximum should be less than the dimension of the system
parameter.plot_onoff = 'on';  % turn on plot when you don't run PaperFigure.m
parameter.processnoisestd = processnoisestd;    %% dynamical noise
parameter.measurementnoisestd = measurementnoisestd; % percentage of the measurement noise

%% Computation
T_START = tic;

data = example(parameter, which_var);
partition = data.partition ; grouping = data.grouping; y = data.y; A = data.A; x_true = data. x_true;

[y, A] = truncate_multiple(y, A, grouping, partition, 1);

% run('RUN_SOLVER.m')
lambda = max(0.005, 20*max(parameter.processnoisestd, parameter.measurementnoisestd));
result = HeteroID(A, y, x_true, lambda, partition, parameter);

T_TOTAL = toc(T_START);
