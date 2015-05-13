fprintf('Simulation Finishes ......\n');
%% =========================== DATA PREPARATION ===========================

parameter.QUIET = 0;

if ~norm(y)
    error('The output is zero vector, please check!!');
end
%% data processing
% plot the oringnal and transformed data
noise = [];
[y_d, A_d] = diffdata(y, A, noise, parameter.diff_onoff);

fprintf('Initialization and data processing finish. \n');
%%
% A_oringinal = A;
[parameter.observation_dim,parameter.basis_dim] = size(A_d);
train_length = round(parameter.prop*parameter.observation_dim);
A_train = A_d(1:train_length,:);
y_train = y_d(1:train_length,:);
A_test = A_d(1+train_length:end,:);
y_test = y_d(1+train_length:end,:);

% x_ls = lscov(A_d,y_d);
% [S,V,D] = svd(A_train);
% % y = 1./diag(svd(A))*S'*y;
% % A = spones(V)*D;
% y = S'*y;
% A = V*D;

%% ===========================   SOLVER BEGIN ===========================

t_solver = tic;
parameter = orderfields(parameter);

if size(A_train,1) > size(A_train,2)
    parameter.method =  'cholesky';
end

switch solver
    case 'cvx'
        result = cvxrwl1(A_train, y_train, parameter);
        % result = rwl1(A_train,y_train, parameter, distributed_onoff);
    case 'admm'
        result = rwl1(A_train,y_train, parameter, distributed_onoff);
end

result.t_solver = toc(t_solver);
result = resultreport(result, y_train, A_train, x_true, parameter);

result.rnmse_prediction_test_true = norm(y_test/parameter.scale - A_test*result.Compare(:,1),2);
result.rnmse_prediction_test_iter_1 = norm(y_test/parameter.scale - A_test*result.Compare(:,2),2);
result.rnmse_prediction_test_iter_end = norm(y_test/parameter.scale - A_test*result.Compare(:,end),2);
result.rnmse_prediction_test_true_normal = norm(y_test/parameter.scale - A_test*result.Compare(:,1),2)/...
    norm(y_test/parameter.scale, 2);
result.rnmse_prediction_test_iter_1_normal = norm(y_test/parameter.scale - A_test*result.Compare(:,2),2)/...
    norm(y_test/parameter.scale, 2);
result.rnmse_prediction_test_iter_end_normal = ...
    norm(y_test/parameter.scale - A_test*result.Compare(:,end),2)/norm(y_test/parameter.scale, 2);

%                         result.fileinfo = fileinfo;
result.t_total = toc(t_start);
nmse_w_lasso = norm(result.Compare(:,2) - x_true)/norm(x_true);
nmse_w_sbl = norm(result.Compare(:,end) - x_true)/norm(x_true);
if ~parameter.QUIET
    fprintf('Data preparation time is %g seconds. \n', result.t_total-result.t_solver);
    fprintf('Solver computation time = %g seconds\n', result.t_solver);
    fprintf('Likelihood train normal is %f. \n',result.likelihood_train_normal );
    fprintf('nmse_w_lasso  = %g\n', nmse_w_lasso);
    fprintf('nmse_w_sbl  = %g\n', nmse_w_sbl);
end
c = datestr(datenum(clock), 'yyyymmddHHMMSS');
parameter.savename=strcat('Data/result/',example,'_result_',solver,'_D',distributed_onoff,...
    '_N', num2str(result.size(2)), '_M', num2str(result.size(1)),'_K',num2str(nnz(x_true)),'_',c);
figurename=strcat(parameter.savename,'.png');

% Plotting
if ~parameter.QUIET
    parameter.plot_onoff = 'on';
    plotdata(result, y_d,  A_d, parameter, figurename);
end

savename=strcat(parameter.savename,'.mat');
result.sol = result.Compare(:,end);
result = orderfields(result);
% save(savename)
