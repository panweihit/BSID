function [result, history] = cvxrwl1(A,y, parameter, distributed_onoff)
% rwl1  Solve reweighted lasso problem
% If the problem size(A,2) <10^2, it will be solved by centralised algorithm
% otherwise it's solved by ADMM splitting
%
% = = = = = = = = = = = = = = = = = = = = = = = Input Argument = = = = = = = = = = = = = = = = = = = = = = = =
% ni:       number of mpi cores
%
%This code is mainly on the update of the weight, for the solver, please
% refer to dwl1
%
% [x, history] = rwl1(A,y, x_true, ni, snr, lambda, rho,alpha);
%
% solves the following problem via ADMM:
%
%   minimize 1/2*|| Ax - b ||_2^2 + \lambda ||U^{k} x||_1
%
%
% The solution is returned in the vector x.
%
% history is a structure that contains the objective value, the primal and
% dual residual norms, and the tolerances for the primal and dual residual
% norms at each iteration.
%
% rho is the augmented Lagrangian parameter.
%
% alpha is the over-relaxation parameter (typical values for alpha are
% between 1.0 and 1.8).
%
%
% w.pan11@imperial.ac.uk
%

% % start of the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

QUIET2 = 1;

if ~QUIET2
    fprintf('This is rwl1 \n');
end

% global parameter
QUIET = 1; % plot or not

%% distributed computation setting

% switch distributed_onoff
%     case 'on'
%
%         if parameter.info.NumProcessors >2
%             v = version('-release');
%             if strcmp(v,'2014a') == 1
%                 parpool(parameter.core);
%             elseif strcmp(v,'2012a') == 1
%                 matlabpool close
%                 matlabpool(parameter.core);
%             end
%             [m,n] = size(A_d);
%
%             if m>n
%                 parameter.dimblock = floor( n/(2*core));
%             elseif m<=n
%                 tmp1 = floor(n/m);
%                 tmp2 = floor(tmp1/parameter.core);
%                 parameter.dimblock = m*tmp2;
%             end
%         else
%             parameter.core = 1;
%         end
% end

t_start = tic;
%% function parameters
solver  = parameter.solver;
solver = 'cvx';


Gammai_threshold = parameter.Gammai_threshold;
%sensitivity_threshold = parameter.sensitivity_threshold;
central_threshold = parameter.central_threshold;
lambda = parameter.lambda;
epsilon = parameter.epsilon;
diaginv = parameter.diaginv;
prune = parameter.prune;
scale = parameter.scale;
scale_w = parameter.scale_w;
info = parameter.info;
truncate = parameter.truncate;
MAX_ITER_Reweight = parameter.MAX_ITER_Reweight; % re-weighted iteration number
Anorm_onoff = parameter.Anorm_onoff;
group.group_onoff = parameter.group_onoff;
group.group_type = parameter.group_type ;


%% Parameter Initialisation

[m,n] = size(A);
INDEX = [1:n];
x = zeros(n,MAX_ITER_Reweight);
result.Estimate = zeros(n,MAX_ITER_Reweight);



if ~QUIET2
    fprintf('Hi,so far good \n');
end

%% Normalize
switch Anorm_onoff
    case 'on'
        Anorm = spdiags(1./sqrt(sum(A.^2))',0,size(A,2),size(A,2)); % normalize columns
        A = A*Anorm;
    case 'off'
        Anorm = 1;
        A = A*Anorm; % initializae featrue matrix for the following
        % x_true1 = x_true./diag(Anorm); % convert the new weight vector
end

switch parameter.tmp
    case 'on'
        tmp = randn(m,m);
        tmpnorm = spdiags(1./sqrt(sum(tmp.^2))',0,m,m); % normalize columns
        tmp = tmp*tmpnorm; % initializae featrue matrix for the following
        
    case 'off'
        tmp = 1;
end
% tmp = tmp*spdiags(1./sqrt(sum(tmp.^2))',0,m,m);
A = tmp*scale*A;
y= tmp*scale*y;
y = scale_w*y;
clear tmp
% result.sol_ls =  (A'*A)\(A'*y); % this is the least square solution, which we should use for comprision
%%%%%%%%%%%%%%%%%%% Specification of B%%%%%%%%%%%%%%%%%%%%%
Btype = 'vector';

switch  Btype
    case 'matrix'
        B =speye(n);
        u = ones(n,1);   % f should have the same size with size(B,1)
        F = diag(u)*B;
    case 'vector'
        B =1;
        u = ones(n,1);   % f should have the same size with size(B,1)
        F = u;
end

if ~QUIET2
    fprintf('Hi, variable initilization complete. \n As you imagin, memory ....\n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result.t_reweight = zeros(MAX_ITER_Reweight, 1);
result.t_optimization = zeros(MAX_ITER_Reweight, 1);
result.t_initialization = toc(t_start);

for iter = 1:1:MAX_ITER_Reweight
    
    
    if ~QUIET
        fprintf('Iteration %d: .\n',iter);
        fprintf('The size of F is %d.\n', size(F,1));
    end
    % = = = = = = =  = = = = = Solver = = = = = =  = = = = = = = = = = = =
    
    switch solver
        case 'admm'
            %%
            switch distributed_onoff
                case 'on'
                    if info.NumProcessors > 2
                        if size(A,2)<=central_threshold
                            % Attention: this is maninly for the updated weight,
                            % typically it will be small than 1000, then we use the centralised
                            % algoithm
                            if ~QUIET
                                fprintf('Centralised weighted lasso\n');
                            end
                            t = tic;
                            [w, history] = wl1(A, y,F, parameter);
                            result.history{iter} = history;
                            result.t_optimization(iter)  = toc(t);
                        else
                            % This is the distriburted weighted lasso
                            if ~QUIET
                                fprintf('Distriburted weighted lasso\n');
                            end
                            [T, history] = dwl1(A, y, F,  parameter);
                            result.history{iter} = history;
                            w = T.w;
                            result.T{iter} = T;
                            result.t_optimization(iter)  = T.t_dwl1-T.t_parfor; % deduct the parfor time, so just count others
                        end
                    else
                        if ~QUIET
                            fprintf('centralised weighted lasso\n');
                        end
                        t = tic;
                        [w, history] = wl1(A, y,F, parameter);
                        result.history{iter} = history;
                        result.t_optimization(iter)  = toc(t);
                    end
                    
                case 'off'
                    
                    if ~QUIET
                        fprintf('centralised weighted lasso\n');
                    end
                    t = tic;
                    [w, history] = wl1(A, y,F, parameter);
                    result.history{iter} = history;
                    result.t_optimization(iter)  = toc(t);
            end
            
            result.history{iter} = history;
            
        case 'cvx'
            %%
            w = cvxwl1(A,y,F,lambda);
    end
    
    t_reweight_start = tic;
    
    % = = = = = = = = = = = = = = = = Prune parameters = = = = = = = = = = = = = = = = = = = = =
    % Prune parameters with low sensitivities
    
    if strcmp(prune,'gamma')
        Gammai=abs(B*w)./u;
        result.Gammai{iter} = Gammai;
        if strcmp(truncate,'on')
            Gammai(Gammai<Gammai_threshold(iter)) = 0;
        end
        index=find(Gammai);
    end
    
    
    %% ====================Update the re-weights========================
    
    % Prune size and generate new feature lighter feature matrix
    
    if strcmp(truncate,'on')
        Gammai = Gammai(index);
        
        if norm(B-eye(size(B))) == 0
            A=A(:,index);
            B =1;
        else
            A = A;
            B = B(index,:);
        end
        % B = B(index,:);
        %  For Bw
    end
    
    % Calculation of the weights $\sqrt{\alpha}$
    
    % when B is diagnol or ones
    if strcmp(diaginv,'woodbury')
        %         size(A)
        %         size(Gammai)
        %     Dic0inv=inv(lambda*speye(m)+A*diag(Gammai)*A');
        VarDiag = woodbury2(A,  Gammai, lambda);
        % get the weights
        u=sqrt(diag(A'*VarDiag*A));
        
    elseif strcmp(diaginv,'lanzcos')
        
        lanzcos_iter = 200;
        P = diag(1./Gammai);
        % lanczos: Computes diag( B*inv(VarDiag)*B' ) with A = A'*R*A + B'*P*B.
        VarDiag = diaginv_lanczos(A, lambda^-1,B,P, lanzcos_iter);
        u= sqrt(- VarDiag./Gammai.^2+1./Gammai);
        
    elseif strcmp(diaginv,'full')
        
        P = diag(1./Gammai);
        % full
        VarDiag = diag(B*inv(lambda^-1*A'*A+B'*P*B)*B');
        u= sqrt(- VarDiag./Gammai.^2+1./Gammai);
        
    end
    
    %     result.VarDiag{iter} = VarDiag;
    result.u{iter} = u;
    if size(B,2)>1
        F = diag(u)*B;
    elseif size(B,2)==1
        F = u.*B;
    end
    
    % ======= Truncate ==========
    
    if strcmp(truncate,'on')
        
        Index{iter}=index;
        
        ind1=Index{iter};
        for j=1:iter-1
            ind2=Index{iter-j};
            ind1=insideback(ind2,ind1); % check the function in the folder
        end
        
        x(INDEX(ind1),iter)=w(index);
        
    elseif strcmp(truncate,'off')
        x (:,iter)= w;
    end
    
    Estimate(:,iter) = sparse(Anorm*x(:,iter))/scale_w;
    
    if ~QUIET
        reportnnz(iter) = nnz(x(:,iter));
        fprintf('rwl1: found %d nonzeros out of R^%d \n',reportnnz(iter),n);
    end
    
    result.t_reweight(iter) = toc(t_reweight_start );
end


result.Estimate = sparse(Estimate);
result.time_initialization = result.t_initialization;
result.time_reweight = sum(result.t_reweight);
result.time_without_split = sum(result.t_optimization);
result.time_rwl1 = result.t_initialization + sum(result.t_reweight)+sum(result.t_optimization);
%% Report

% if ~QUIET
%     reportnnz(iter) = length(find(x(:,iter)));
%     fprintf(1,['Found a feasible x in R^%d that has %d nonzeros ' ...
%         'using reweighted l1 at the weight space.\n'],n,reportnnz(iter));
% end


% ===================== break iteration ================================
%     result.nmse_y(iter,1)= norm(y - A*w(index))^2./norm(y)^2;

%     if result.nmse_y(iter) < 1e-2
%         break;
%     end

%% Result

fprintf('The rwl1 solver time is %f seconds. \n',result.time_rwl1);
result.time_look_at_this = result.time_rwl1;
end



function W = cvxwl1(A,y,F, lambda)
fprintf('You are using CVX wl1. \n');
% = = = = = = =  = = = = = Solver CVX = = = = = =  = = = = = = = = = = = =
n = size(A,2);
F = diag(F);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cvx_begin quiet
cvx_solver sedumi
%     cvx_solver SDPT3
variable w(n)
minimize    (lambda*norm( F*w, 1 )+ 0.5*sum((A* w-y).^2) )
% minimize    (lambda*norm( F*w, 1 )+ 0.5*quad_form(A* w-y, Covinv) )
%                      subject to

cvx_end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% = = = = = = =  = = = = = CVX end = = = = = =  = = = = = = = = = = = =
W = w;

end
% snr2 = 20*log10(norm(AA*x_true)/norm(noise));
