%% CONMENT THE FOLLOWING LINES IF YOU RUN PaperFigure.m
close all; 
clear;
T = 100; % length of a single time series
replicate = 10;  %number of experiments/conditions/pertubations/replicates; 
                           %when replicate = 1, the algorithm is consistent with our TAC 2016 paper
which_gene = 1;  % Specify which gene you are interested in. 
                                 % The maximum should be less than the dimension of the system
parameter.plot_onoff = 'on';  % turn on plot when you don't run PaperFigure.m
parameter.dynamicnoise = 0;    %% dynamical noise
parameter.measurementnoiseperc = 0; % percentage of the measurement noise

%% Parameter Initialization
% You have to tune the tradeoff parameter lambda!!
parameter.lambda = 0.001;

parameter.MAX_ITER_Reweight =6;
parameter.Gammai_threshold = [1e-3  1e-3  1e-3 1e-3 1e-3 1e-3 1e-2 1e-1 1e-1];
parameter.replicate = replicate; 
parameter.T = T; 
parameter.diff_onoff = 'off';
parameter.truncate = 'on';
parameter.tmp = 'off';
parameter.info = cpuinfo;
parameter.core = 5;
parameter.rho_update = 'off';
parameter.rho =  1;
parameter.alpha = 0.7;
parameter.epsilon = 1e-10;
parameter.scale = 1e0;
parameter.scale_w = 1e0;
parameter.dimblock = 5e3; %dimension of the each block
parameter.central_threshold = 5e3;   % threshold for using centralised algorithm even distributed_onoff is 'on'
parameter.scale_w = 1e0;
parameter.Anorm_onoff = 'off';
diaginv = {'woodbury', 'lanzcos', 'full'};
prune = {'gamma','sensitivity','rank'};
method = {'woodbury', 'cholesky', 'lcg', 'full'};
parameter.diaginv =diaginv{1};
parameter.prune = prune{1};
parameter.method = method{1};
parameter.group_type = 'SMV';
parameter.multiple = 'on';
parameter.single = 'off';
parameter.group_onoff = parameter.multiple;
parameter.derivative ='solver'; % change this accordingly
parameter.derivativekth = 1;
parameter.prop  = 1; 
prop = parameter.prop;
distributed_onoff = 'off'; % If you switch to admm solver, you can enable distributed computation by specifying multi cores
solver = {'cvx','admm'}; % Switch which solver you want to use; 
                                             % CVX is using CVX toolbox calling Sedumi;
                                             % admm is customed implemented, which is faster

%% Example
solver = solver{2};
parameter.solver = solver;
% select = 'generalrepressilator_continuous_even'; % continuous system with even number of state variables, ode45 is used for simulation
select = 'generalrepressilator_continuous_odd'; % continuous system with odd number of state variables, ode45 is used for simulation
% select = 'generalrepressilator_discrete_even'; % discrete system with even number of state variables, direct simulation
% select = 'generalrepressilator_discrete_odd'; % discrete system with odd number of state variables, direct simulation

%% Computation
T_START = tic;
t_start = tic;

data = example_multiple(select, parameter, which_gene);
partition = data.partition ; grouping = data.grouping; y = data.y; A = data.A; x_true = data. x_true;

[y, A] = truncate_multiple(y, A, grouping, partition, 1);
run('RUN_SOLVER.m')


T_TOTAL = toc(T_START)
