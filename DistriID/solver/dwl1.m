function [T,history] = dwl1(A, b, F, parameter)
% dwl1  Solve group distributed weighted lasso problem via ADMM feature splitting
% It should be noticed here there is no weight updating here
% Weight update is included in rwl1
%
% The block sizes n_i, so that x_i is in R^{n_i}.
%
% solves the following problem via ADMM:
%
%   minimize 1/2*|| Ax - b ||_2^2 + \lambda ||Fs x||_1
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
% ===================Author
% w.pan11@imperial.ac.uk
fprintf('you are using dwl1\n');

if size(F,2) >1 && nnz(F) == nnz(diag(F))
    F = diag(F);
end

%% Parameter Initialisation
lambda = parameter.lambda;
rho = parameter.rho;
alpha = parameter.alpha;
rho_update = parameter.rho_update;
dimblock = parameter.dimblock;
% Global constants and defaults

QUIET    = 1;
MAX_ITER = 200;
RELTOL  = 1e-2;
ABSTOL   = 1e-2;

%% Data preprocessing



[m, n] = size(A); % size of the whole system
% dimblock = 10; % dimension of subsystem

N = floor(n/dimblock); % number of subsystems
if N == 0
    ni = n;
    N = 1;
else
    for i=1:N
        ni(i) = dimblock;
    end
    res = rem(n, dimblock);
    ni(N) = ni(N)+res
end

% check that ni divides in to n

if (sum(ni) ~= n)
    error('invalid block size');
end

N
%% ADMM solver
t_start = tic;

% initialization
% x = zeros(sum(ni),N);
z = zeros(m,1);
uu = zeros(m,1);
Axbar = zeros(m,1);

zs = zeros(m,N);
Aixi = zeros(m,N);

if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
        'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end

% pre-factor

t0=tic;
start_ind = 1;
cum_part = cumsum(ni);
for i = 1:N
    % in Matlab, transposing costs space and flops
    % so we save a transpose operation everytime
    
    %     fi{i}=f((i-1)*ni(i) + 1:i*ni(i),1);
    fi{i} = F(start_ind:cum_part(i),1);
    %     Ai = A(:,(i-1)*ni(i) + 1:i*ni(i));
    Ai = A(:,start_ind:cum_part(i));
    start_ind = cum_part(i)+1;
    
    At{i}=Ai';
    AA{i} = Ai;
    
    m_tmp = size(Ai,1);
    
    if m_tmp ~=m
        error('m is wrong');
    end
    
    % check dimension use woodbury identity
    %     fi=ui{i};
    %     AAibarinv{i}=inv(Ai'*Ai + rho*sparse(diag(fi.*fi)));
    %     AAibarinv{i} = woodbury(Ai, fi{i}, rho);
    
end
t_cache=toc(t0);

t_parfor = zeros(MAX_ITER,1);
t_gather = zeros(MAX_ITER,1);


fprintf('You have cached everything. \n');

for iter = 1:MAX_ITER
    %% =============Parallel Computation===============
    if ~QUIET
        fprintf('You run parellell of iteration %d. \n',iter);
    end
    % x-update (to be done in parallel)
    t_parfor_start=tic; % par_update
    
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%     for i = 1:N
    parfor i = 1:N
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        % x-update
        
        Ai=AA{i};
        %         Aibarinv=AAibarinv{i};
        %         fprintf(1,'\n This is round %d for x-update in the %d-th iteration.\n',i,iter);
        xx = x_update(Ai,  Aixi(:,i) + z - Axbar - uu, fi{i}, parameter);
        x{i} = xx;
        Aixi(:,i) = Ai*xx;
        
    end
    
    t_parfor(iter)=toc(t_parfor_start);
    
    if ~QUIET
        fprintf('You run garther of iteration %d. \n',iter);
    end
    
    t_gather_start=tic;  % communication
    % z-update
    zold = z;
    Axbar = 1/N*A*vectorize(x);
    
    Axbar_hat = alpha*Axbar + (1-alpha)*zold;
    z = (b + rho*(Axbar_hat + uu))/(N+rho);
    
    % u-update
    uu = uu + Axbar_hat - z;
    
    % compute the dual residual norm square
    s = 0; q = 0;
    zsold = zs;
    zs = z*ones(1,N) + Aixi - Axbar*ones(1,N);
    for i = 1:N,
        % dual residual norm square
        s = s + norm(-rho*At{i}*(zs(:,i) - zsold(:,i)))^2;
        % dual residual epsilon
        q = q + norm(rho*At{i}*uu)^2;
    end
    
    % diagnostics, reporting, termination checks
    history.objval(iter)  = objective(A, b, lambda, cum_part, N, x, z);
    history.r_norm(iter)  = sqrt(N)*norm(z - Axbar);
    history.s_norm(iter)  = sqrt(s);
    
    
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
    
    
    
    history.eps_pri(iter) = sqrt(n)*ABSTOL + RELTOL*max(norm(Aixi,'fro'), norm(-zs, 'fro'));
    history.eps_dual(iter)= sqrt(n)*ABSTOL + RELTOL*sqrt(q);
    
    
    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', iter, ...
            history.r_norm(iter), history.eps_pri(iter), ...
            history.s_norm(iter), history.eps_dual(iter), history.objval(iter));
    end
    
    if history.r_norm(iter) < history.eps_pri(iter) && ...
            history.s_norm(iter) < history.eps_dual(iter);
        break
    end
    t_gather(iter)=toc(t_gather_start);
end


t_total=toc(t_start);
T.t_total = t_total;
T.t_cache = t_cache;
T.t_gather = sum(t_gather);
T.t_parfor = sum(t_parfor);

T.t_dwl1 = T.t_cache+T.t_gather + T.t_parfor;
% T1=[t_cache, t_end]
% t_cache-> record the the time for splitting to A_i and the inversion time; t_end -> total time of dwl1
% T2=[t_xupdate, t_center]
% It records the x-update and centralised computation in each iteration

T.w=vectorize(x);
end




% function p = objective(A, b, lambda, N, x, z)
%     p = ( 1/2*sum_square(N*z - b) + lambda*sum(norms(x)) );
% end

function p = objective(A, b, lambda, cum_part, N, x, z)
obj = 0;
xx = vectorize(x);
start_ind = 1;
for i = 1:length(cum_part),
    sel = start_ind:cum_part(i);
    obj = obj + norm(xx(sel),1);
    start_ind = cum_part(i) + 1;
end
p = ( 1/2*sum_square(N*z - b) + lambda*obj );

end


function x = x_update(A,  v, f, parameter)

[m,n] = size(A);

% A = multidiaginv(A,1./f);

kappa=parameter.lambda/parameter.rho;

q = A'*v;

if (norm(q) <= kappa)
    x = zeros(n,1);
else
    %     x=l1_ls(A,v,2*kappa);
    x =wl1(A, v,f, parameter); % you can replace this by any other $\ell_1$ algorithm
end

end


