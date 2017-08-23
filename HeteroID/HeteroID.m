function result = HeteroID(A, y, x_true, lambda, partition, parameter)
t_start = tic;
if ~exist('x_true')
    x_true = zeros(size(A, 2),1);
end
run('config.m');
parameter.lambda =  lambda;
%% =========================== DATA PREPARATION ===========================
parameter.QUIET = 1;
if ~norm(y)
    error('The output is zero vector, please check!!');
end
%% data processing
% plot the oringnal and transformed data
noise = [];
[y_d, A_d] = diffdata(y, A, noise, parameter.diff_onoff);
% fprintf('Simulation Finishes ......\n');
% fprintf('Initialization and data processing finish. \n');
% random shuffling

%%
% A_oringinal = A;
[parameter.observation_dim,parameter.basis_dim] = size(A_d);
train_length = round(parameter.prop*parameter.observation_dim);

A_train = A_d(1:train_length,:);
y_train = y_d(1:train_length,:);

% index = randperm(length(y_train));
% y_train = y_train(index);
% A_train = A_train(index, :);

A_test = A_d(1+train_length:end,:);
y_test = y_d(1+train_length:end,:);

%% ===========================   SOLVER BEGIN ===========================

t_solver = tic;

result = crwgl1(A_train, y_train(:), partition, parameter);
result.t_solver = toc(t_solver);

result.partition = partition;
result = resultreport(result, y_train(:), A_train, y_test(:), A_test, x_true, parameter);


result.t_total = toc(t_start);
err = norm( result.Compare(:,end) - x_true , 'fro') / norm(result.Compare(:,1) , 'fro');

if strcmp(parameter.plot_onoff, 'on')
    fprintf('Relative Estimation error = %g\n', err);
    fprintf('Data preparation time is %g seconds. \n', result.t_total-result.t_solver);
    fprintf('Solver computation time = %g seconds\n', result.t_solver);
    fprintf('Likelihood train normal is %f. \n',result.likelihood_train_normal );
    fprintf('M=%d, N=%d, K=%d, lambda=%d \n',...
        round(parameter.observation_dim),parameter.basis_dim,result.sparsity,parameter.lambda);
    fprintf('Likelihood train normal is %f. \n',result.likelihood_train_normal );
    fprintf('rnmse_w_iter_1 is %f. \n',result.rnmse_w_iter_1);
    fprintf('rnmse_w_iter_end is %f. \n',result.rnmse_w_iter_end);
end

c = datestr(datenum(clock), 'yyyymmddHHMMSS');
% parameter.savename=strcat('data/result/',select,'_',solver,'_Distr',distributed_onoff,...
%     '_N', num2str(result.size(2)), '_M', num2str(result.size(1)),'_K',num2str(nnz(x_true)),'_',c);
parameter.savename=strcat('data/result/',solver,'_Distr',distributed_onoff,...
    '_N', num2str(result.size(2)), '_M', num2str(result.size(1)),'_K',num2str(nnz(x_true)),'_',c);
figurename=strcat(parameter.savename,'.png');
savename = strcat(parameter.savename,'.mat');

% Plotting
if strcmp(parameter.plot_onoff, 'on')
    plotdata(result, y_d,  A_d, parameter, figurename);
end

clear interval

result.w = result.Compare(:,end);
result = orderfields(result);
result.parameter = orderfields(result.parameter);

