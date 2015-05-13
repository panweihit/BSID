function result = resultreport(result, y, A, x_true, parameter)


if strcmp(parameter.group_type,'MMV')
    
    [p,L] = size(y);
    y = vec(y');
    A = kron(A,eye(L));
    x_true = vec(x_true');
    
end

QUIET = 0;
scale = 1;
y = y/scale;
result.Estimate = result.Estimate/scale;
x_true = x_true/scale;
x_Estimate = sparse(result.Estimate);
result.Estimate = [];
[m,n] = size(A);

% %   normalisation
% Anorm = spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalise columns
% A0= A*Anorm; % new feature matrix after normalisation
% A1 = A0; % initialisae featrue matrix for the following
% x_true1 = x_true./diag(Anorm); % convert the new weight vector


%% Result


% diplay result statistics
if ~QUIET
    fprintf('Estimated sparsity is %d in R^%d.\n',nnz(x_Estimate(:,end)),n);
    fprintf('True sparsity is %d in R^%d.\n',nnz(x_true),n);
end

result.size = size(A);
result.parameter = parameter;
result.sparsity = nnz(x_true);
result.sparsity_estimate = nnz(x_Estimate(:,end));
% result.Compare = sparse([x_true,x_Estimate]) ;
result.Compare = full([x_true,x_Estimate]) ;
result.rnmse_w_iter_1 = norm(x_true - x_Estimate(:,1),2)/norm(x_true,2);
result.rnmse_w_iter_end = norm(x_true - x_Estimate(:,end),2)/norm(x_true,2);
compare= sparse([x_true,x_Estimate(:,end)]);
result.likelihood_train = norm(y - A*compare(:,end),2);
result.likelihood_train_normal = norm(y - A*compare(:,end),2)/norm(y,2);
% result.x_true = sparse(x_true);
result.sparsity_estimate = nnz(compare(:,end));


result = orderfields(result);


% snr2 = 20*log10(norm(AA*x_true)/norm(noise));

%     figure;
%     plot( result.compare_admm(:,1),'r','LineWidth', 1.5);
%     hold on;
%     plot(result.compare_admm(:,2),'b o-','LineWidth', 1.5)
%     legend('true','estimate');
