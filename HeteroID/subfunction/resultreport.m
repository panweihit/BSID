function result = resultreport(result, y_train, A_train, y_test, A_test, x_true, parameter)

if strcmp(parameter.group_type,'MMV')   
    [p,L] = size(y_train);
    y_train = vec(y_train');
    A_train = kron(A_train,eye(L));
    x_true = vec(x_true'); 
end

QUIET = 0;
scale = 1;
y_train = y_train/scale;
result.Estimate = result.Estimate/scale;
x_true = x_true/scale;
x_Estimate = sparse(result.Estimate);
result.Estimate = [];
[m,n] = size(A_train);


%% Result


% diplay result statistics
if ~QUIET
    fprintf('Estimated sparsity is %d in R^%d.\n',nnz(x_Estimate(:,end)),n);
    fprintf('True sparsity is %d in R^%d.\n',nnz(x_true),n);
end

result.rnmse_w_iter_1 = norm(x_true - x_Estimate(:,1),2)/norm(x_true,2);
result.rnmse_w_iter_end = norm(x_true - x_Estimate(:,end),2)/norm(x_true,2);
result.size = size(A_train);
result.parameter = parameter;
result.sparsity = nnz(x_true);
result.sparsity_estimate = nnz(x_Estimate(:,end));
% result.Compare = sparse([x_true,x_Estimate]) ;
result.Compare = full([x_true,x_Estimate]) ;
compare= sparse([x_true,x_Estimate(:,end)]);
result.likelihood_train = norm(y_train - A_train*compare(:,end),2);
result.likelihood_train_normal = norm(y_train - A_train*compare(:,end),2)/norm(y_train,2);
% result.x_true = sparse(x_true);
result.sparsity_estimate = nnz(compare(:,end));


result.rnmse_prediction_test_true = norm(y_test/parameter.scale - A_test*result.Compare(:,1),2);
result.rnmse_prediction_test_iter_1 = norm(y_test/parameter.scale - A_test*result.Compare(:,2),2);
result.rnmse_prediction_test_iter_end = norm(y_test/parameter.scale - A_test*result.Compare(:,end),2);
result.rnmse_prediction_test_true_normal = norm(y_test/parameter.scale - A_test*result.Compare(:,1),2)/...
    norm(y_test/parameter.scale, 2);
result.rnmse_prediction_test_iter_1_normal = ...
    norm(y_test/parameter.scale - A_test*result.Compare(:,2),2)/...
    norm(y_test/parameter.scale, 2);
result.rnmse_prediction_test_iter_end_normal = ...
    norm(y_test/parameter.scale - A_test*result.Compare(:,end),2)/...
    norm(y_test/parameter.scale, 2);


result = orderfields(result);


