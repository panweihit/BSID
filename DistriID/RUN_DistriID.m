%% initialization
close all;
clear all;
t_start = tic;

%% Distirbuted Computation Setup
% Do you want to distributed computation or centralised compuataion?
distributed_onoff = 'on';
if strcmp(distributed_onoff,'on')
    delete(gcp);
    parpool(15);
end

%dimension of the each block you want to put
parameter.dimblock = 1000;
% threshold for using centralised algorithm even distributed_onoff is 'on'
parameter.central_threshold = 1000;

%%
diaginv = {'woodbury', 'lanzcos', 'full'};
prune = {'gamma'};
method = {'woodbury', 'cholesky', 'full'};
parameter.truncate = 'on';
parameter.Anorm_onoff = 'off';
parameter.tmp = 'off';
parameter.info = cpuinfo;
parameter.core = 5;
parameter.rho_update = 'off';
parameter.rho =  1;
parameter.alpha = 1.5;
parameter.epsilon = 1e-12;
parameter.scale = 1e0;
parameter.scale_w = 1e0;
parameter.diff_onoff = 'off';
parameter.MAX_ITER_Reweight = 5;

solver = {'cvx','admm'};
% example = 'Kuramoto'
% example = 'Repressilator'
example = 'random'

% proportion for traning set and test set.
%5/10 means half for training and half for testing
parameter.prop = 5/10;
prop = parameter.prop;
%%
switch example
    case 'Kuramoto'
        %%
        solver = solver{2};
        parameter.solver = solver;
        parameter.group_onoff = 'off';
        parameter.group_type = 'SMV';
        parameter.scale_w = 1e0;
        parameter.Anorm_onoff = 'on';
        parameter.diaginv =diaginv{1};
        parameter.prune = prune{1};
        parameter.method = method{1};
        parameter.MAX_ITER_Reweight = 5;
        parameter.truncate = 'on';
        parameter.Gammai_threshold = 1e0* [1e-3 1e-3 1e-4 1e-4 1e-4 1e-2 1e-2 1e-1 1e-1 1e-1 1e-1];
        sigma_it =-10;
        %         sigma =10^(sigma_it*1e-1);
        sigma = .05;
        N = 1000; % dimension of the system
        [y, A, x_true, snr]=discrete_K(N,0.1, 20,sigma, pi*(1+randn(N,1)),1,1,5);
        parameter.lambda =  0.1*sigma*norm(A'*y,inf);
        %         parameter.lambda = 1;
        parameter.diaginv = diaginv{1};
        parameter.prune = prune{1};
        parameter.method = method{1};
        
    case 'random'
        %%
        solver = solver{2};
        parameter.solver = solver;
        parameter.group_onoff = 'off';
        parameter.group_type = 'SMV';
        parameter.scale_w = 1e0;
        parameter.Anorm_onoff = 'on';
        parameter.diaginv =diaginv{1};
        parameter.prune = prune{1};
        parameter.method = method{1};
        parameter.MAX_ITER_Reweight = 5;
        parameter.truncate = 'on';
        parameter.Gammai_threshold = 1e0* [1e-3 1e-3 1e-4 1e-4 1e-4 1e-2 1e-2 1e-1 1e-1 1e-1 1e-1];
        sigma_it =-10;
        %         sigma =10^(sigma_it*1e-1);
        sigma = .05;
        N = 10000; % dimension of the system
        A = randn(5000,N);
        x_true = sprandn(N,1,0.01);
        y = A*x_true;
        
        parameter.lambda =  0.01*sigma*norm(A'*y,inf);
        %         parameter.lambda = 1;
        parameter.diaginv = diaginv{1};
        parameter.prune = prune{1};
        parameter.method = method{1};
end

%%
run('RUN_SOLVER.m')
%%


