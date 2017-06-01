function [x, history] = wgl1(A,  y, F, partition, parameter)
% function [x, history] = wl1(A, y,F, parameter)
% Weighted Group Lasso

% lasso  Solve lasso problem via ADMM
%
% [z, history] = lasso(A, b, lambda, rho, alpha);
%
% Solves the following problem via ADMM:
%
%   minimize 1/2*|| Ax - b ||_2^2 + \lambda ||F* x ||_2
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
%Contact: w.pan11@imperial.ac.uk
%%
% f is a VECTOR NOT MATRIX!!!!!!!!!!!!  you have to use diag(VECTOR)
%  if f is a vector or diagonal matrix


% check if F is diagonal or not, if it is, we abstract the vector f
if size(F,2) >1 && nnz(F) == nnz(diag(F))
    F = diag(F);
end
%%
lambda = parameter.lambda;
% F = lambda*F;
rho = parameter.rho;
alpha = parameter.alpha;
method =  parameter.method;
rho_update = parameter.rho_update;

t_start = tic;

%% Global constants and defaults

QUIET    = 1;
MAX_ITER =1000;
ABSTOL   = 1e-7;
RELTOL   = 1e-7;

%% Data preprocessing
[m_A, n_A] = size(A);

% % check that sum(p) = total number of elements in x
% if (sum(partition) ~= n)
%     error('invalid partition');
% end



f= []; 
for i =1:length(partition)
    fi = F(i)*ones(partition(i),1);
    f = [f;fi];
end
nF = size(f,1);

% cumulative partition
cum_part = cumsum(partition);


%% save a matrix-vector multiply

Atb = A'*y;
AtA = A'*A;
%% ADMM solver

x = zeros(n_A,1);
z = zeros(nF,1);
u = zeros(n_A,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(F,2) ==1 % VECTOR or DIAGONAL MATRIX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pre-factor

    F=sparse(diag(f));
    btb = f.*f;
    
    
    if  strcmp(method,'woodbury') 
        Abarinv = woodbury(A, f, rho);
    elseif  strcmp(method,'cholesky') 
        [L, R] = factor(A, f, rho);
    end
    
    
    if ~QUIET
        fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
            'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
    end
    
    Rho.rho = rho;
    Rho.rho_early = rho;
    Rho.rho_update = rho_update;
    Rho.method = method;
    
    for iter = 1:MAX_ITER
        
        % x-update
        q = Atb + rho*f.*(z-u);      % temporary value
        
        if  strcmp(method,'woodbury') 
            x = x_update_woodbury(A,  Abarinv, q, f, Rho);
        elseif  strcmp(method,'cholesky') 
            x = x_update_cholesky(A, L,R, q, f, btb, Rho);
        elseif  strcmp(method,'full') 
            Abar=AtA + rho*diag(btb);
            x = Abar\q;
        elseif  strcmp(method,'lcg') 
            x = lcg(A,1,F,rho,q);
        end
        
        % z-update with relaxation
        
        zold = z;
        start_ind = 1;
        Fx_hat = alpha*f.*x +(1-alpha)*zold;
        for i = 1:length(partition),
            sel = start_ind:cum_part(i);
            z(sel) = shrinkage(Fx_hat(sel) + u(sel), lambda/rho);
            start_ind = cum_part(i) + 1;
        end
        
        
        
        % y-update
        u = u + Fx_hat - z;
        
        % diagnostics, reporting, termination checks
        history.objval(iter)  = objective(A, y, F, lambda, cum_part, x, z);
        
        history.r_norm(iter)  = norm(f.*x - z);
        history.s_norm(iter)  = norm(-rho*f.*(z - zold));
        
        %     rho updata
        if strcmp(rho_update,'on')
            
            history.rho(iter)=rho;
            r_norm(iter) = history.r_norm(iter);
            s_norm(iter) = history.s_norm(iter);
            
            
            miu = 10;
            tau_incr = 2;
            tau_decr = 2;
            
            if r_norm(iter)>miu*s_norm(iter)
                rho=tau_incr*rho;
            elseif s_norm(iter)>miu*r_norm(iter)
                rho=rho/tau_decr;
            else
                rho=rho;
            end
        end
        
        history.rho(iter) = rho;
        if iter > 1
            rho_early  = history.rho(iter-1);
        end
        
        
        history.eps_pri(iter) = sqrt(n_A)*ABSTOL + RELTOL*max(norm(f.*x), norm(-z));
        history.eps_dual(iter)= sqrt(n_A)*ABSTOL + RELTOL*norm(rho*f.*u);
        
        if ~QUIET
            fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', iter, ...
                history.r_norm(iter), history.eps_pri(iter), ...
                history.s_norm(iter), history.eps_dual(iter), history.objval(iter));
        end
        
        if (history.r_norm(iter) < history.eps_pri(iter) && ...
                history.s_norm(iter) < history.eps_dual(iter))
            break;
        end
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif size(F,2) >1  % GENERAL MATRIX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % pre-factor
    if  strcmp(method,'woodbury') 
        F = full(F);
        if size(F,1)~=size(F,2) || norm(diag(diag(F)) -F)~=0
            error('You can not use Woodbury since the Weight matrix is not diagnol.')
        end
        Abarinv = woodbury(A, F, rho);
    elseif  strcmp(method,'cholesky') 
        [L, R] = factor(A, F, rho);
    end
    
    
    if ~QUIET
        fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
            'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
    end
    
    Rho.rho = rho;
    Rho.rho_early = rho;
    Rho.rho_update = rho_update;
    Rho.method = method;
    
    for iter = 1:MAX_ITER
        
        % x-update
        fprintf('check size\n');
        size(Atb)
        size(F)
        size(z)
        size(u)
        q = Atb + rho*F'*(z-u);      % temporary value
        
        if  strcmp(method,'woodbury') 
            x = x_update_woodbury(A,  Abarinv, q, F, Rho);
        elseif  strcmp(method,'cholesky') 
            x = x_update_cholesky(A, L,R, q, F, 0, Rho);
        elseif  strcmp(method,'full') 
            Abar=AtA + rho*F'*F;
            x = Abar\q;
        elseif  strcmp(method,'lcg') 
            x = lcg(A,1,F,rho,q);
        end
        
        % z-update with relaxation
        zold = z;
        start_ind = 1;
        Fx_hat = alpha*F*x +(1-alpha)*zold;
        for i = 1:length(partition),
            sel = start_ind:cum_part(i);
            z(sel) = shrinkage(Fx_hat(sel) + u(sel), lambda/rho);
            start_ind = cum_part(i) + 1;
        end
        
        
        
        
        % z-update
        zold = z;
        start_ind = 1;
        Fx_hat = alpha*F*x +(1-alpha)*zold;
        for i = 1:length(partition),
            sel = start_ind:cum_part(i);
            z(sel) = shrinkage(Fx_hat(sel) + u(sel), lambda/rho);
            start_ind = cum_part(i) + 1;
        end
        
        
        
        
        
        % y-update
        u = u + Fx_hat - z;
        
        % diagnostics, reporting, termination checks
        %         history.objval(iter)  = objective(A, y, lambda, x, z);
        history.objval(iter)  =  objective(A, y, F, lambda, cum_part, x, z);
        
        history.r_norm(iter)  = norm(F*x - z);
        history.s_norm(iter)  = norm(-rho*F'*(z - zold));
        
        %     rho updata
        if strcmp(rho_update,'on')
            
            history.rho(iter)=rho;
            r_norm(iter) = history.r_norm(iter);
            s_norm(iter) = history.s_norm(iter);
            
            
            miu = 10;
            tau_incr = 2;
            tau_decr = 2;
            
            if r_norm(iter)>miu*s_norm(iter)
                rho=tau_incr*rho;
            elseif s_norm(iter)>miu*r_norm(iter)
                rho=rho/tau_decr;
            else
                rho=rho;
            end
        end
        
        history.rho(iter) = rho;
        if iter > 1
            rho_early  = history.rho(iter-1);
        end
        
        
        history.eps_pri(iter) = sqrt(n_A)*ABSTOL + RELTOL*max(norm(F*x), norm(-z));
        history.eps_dual(iter)= sqrt(n_A)*ABSTOL + RELTOL*norm(rho*F'*u);
        
        if ~QUIET
            fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', iter, ...
                history.r_norm(iter), history.eps_pri(iter), ...
                history.s_norm(iter), history.eps_dual(iter), history.objval(iter));
        end
        
        if (history.r_norm(iter) < history.eps_pri(iter) && ...
                history.s_norm(iter) < history.eps_dual(iter))
            break;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


t_iter =toc(t_start);


Rho.rho = rho;
% Rho.rho_early = rho_early;
Rho.rho_update = rho_update;

%[x,z]
end

% function p = objective(A, y, lambda, x, z)
% p = ( 1/2*norm(A*x - y)^2 + lambda*norm(z,1) );
% end
function p = objective(A, y, F, lambda, cum_part, x, z)
obj = 0;
start_ind = 1;
for i = 1:length(cum_part),
    sel = start_ind:cum_part(i);
    obj = obj + F(i)*norm(z(sel),2);
    start_ind = cum_part(i) + 1;
end
p = ( 1/2*sum((A*x - y).^2) + lambda*obj );
end



% function z = shrinkage(x, kappa)
% z = max( 0, x - kappa ) - max( 0, -x - kappa );
% end
function z = shrinkage(x, kappa)
z = pos(1 - kappa/norm(x))*x;
end


function [L, R] = factor(A, f, rho)

%  if f is a vector or diagonal matrix
if size(f,2) ==1
    ftf = f.*f;
    [m, n] = size(A);
    if ( m >= n )    % if skinny
        L = chol( A'*A + rho*sparse(diag(ftf)), 'lower' );
    else            % if fat
        L = chol( speye(m) + 1/rho*(A*sparse(diag(1./ftf))*A'), 'lower' );
    end
    % force matlab to recognize the upper / lower triangular structure
    L = sparse(L);
    R = sparse(L');
    
    % if f is a general matrix
elseif size(f,2) >1
    
    ftf = sparse(f'*f);
    
    [m, n] = size(A);
    if ( m >= n )    % if skinny
        L = chol( A'*A + rho*ftf, 'lower' );
    else            % if fat
        L = chol( A'*A + rho*ftf, 'lower' );
        % the same as above! You can't use the expression below since inv(f'*f) does not
        % exist
        %        L = chol( speye(m) + 1/rho*(A*inv(f'*f)*A'), 'lower' );
    end
    % force matlab to recognize the upper / lower triangular structure
    L = sparse(L);
    R = sparse(L');
    
end
end




function x = x_update_woodbury(A,  Abarinv, q, f, Rho)

% x-update
%     Abarinv=inv(AtA + rho*FtF);
%     x = Abarinv*(Atb + rho*F'*(z-u));

% F must be diagonal!!

% if there is a new rho updated

rho = Rho.rho;
rho_early = Rho.rho_early;
rho_update = Rho.rho_update;
method = Rho.method;

if abs(rho - rho_early) >1e-2
    if strcmp(rho_update,'on')
        % pre-factor
        if  strcmp(method,'woodbury') 
            F=sparse(diag(f));
            Abarinv = woodbury(A, f, rho);
        end
    end
end

% if there is NO new rho updated

if strcmp(method,'woodbury') 
    %      q = Atb + rho*f.*(z-u);    % temporary value
    %      Abarinv = woodbury(A, f, rho);
    x = Abarinv*q;
end

end



function x = x_update_cholesky(A, L,R, q, f, ftf, Rho)
% x-update

% if there is a new rho updated
rho = Rho.rho;
rho_early = Rho.rho_early;
rho_update = Rho.rho_update;
method = Rho.method;
[m,n] = size(A);

if abs(rho - rho_early) >1e-2
    if strcmp(rho_update,'on')
        % pre-factor
        if  strcmp(method,'cholesky') 
            [L, R] = factor(A, f, rho);
        end
    end
end


% if there is NO new rho updated

%  if f is a vector or diagonal matrix
if size(f,2) ==1
    if  strcmp(method,'cholesky') 
        if( m >= n )    % if skinny
            x = R \ (L \ q);
        else            % if fat
            ff = 1./(rho*ftf);
            ffq = ff.*q;
            tmp = A'*(R \ ( L \( A*ffq) ));
            x =ffq - ff.*tmp;
        end
    end
    % if f is a general matrix
elseif size(f,2) >1
    if  strcmp(method,'cholesky') 
        if( m >= n )    % if skinny
            x = R \ (L \ q);
        else            % if fat
            x = R \ (L \ q);
        end
    end
    
end


end







