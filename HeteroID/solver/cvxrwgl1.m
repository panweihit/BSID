function result=cvxrwgl1(A, y, partition, parameter)

% global parameter
QUIET = 0; % plot or not
%% function parameters

group.group_onoff = parameter.group_onoff;
group.group_type = parameter.group_type ;
L = size(y,2);
[m_a,n_a] = size(A);
Aa= A; 
if strcmp(group.group_type,'MMV') == 1 && L>1
    tmpa = 1;%randn(m_a);
    Aa = tmpa*Aa;
    y = tmpa*y;
    group.Aa = Aa;
    group.L = L;
    [y, A,~, ~, partition]= expandmatrix( y, Aa, 0, 0, group);
end

lambda = parameter.lambda;
diaginv = parameter.diaginv;
prune = parameter.prune;
scale = parameter.scale;
MAX_ITER_Reweight = parameter.MAX_ITER_Reweight; % re-weighted iteration number
Gammai_threshold = parameter.Gammai_threshold;
truncate = parameter.truncate;
lambda = lambda*scale;
%% Parameter Initialisation
[m,n]=size(A);
M = m;
p = length(partition);
partition_tmp = partition;
INDEX=[1:n];
x=zeros(n,MAX_ITER_Reweight);

%%   Normalisation
switch Anorm_onoff
    case 'on'
        Anorm = spdiags(1./sqrt(sum(A.^2))',0,size(A,2),size(A,2)); % normalize columns
        A0 = A*Anorm; 
        A = A0;
    case 'off'
        Anorm = 1;
        A0 = A*Anorm; % initializae featrue matrix for the following
        A = A0;
        % x_true1 = x_true./diag(Anorm); % convert the new weight vector
end


switch parameter.tmp
    case 'on'
        tmp = randn(m,m);
    case 'off'
        tmp = 1;
end


A = tmp*A;
y= tmp*scale*y;
        Anorm2 = diag(1./sqrt(norms( A, 2,2 )));
        A = Anorm2*A;
        y = Anorm2*y;
%%%%%%%%%%%%%
% u=abs(randn(p,1)); %
u=ones(p,1); %
%%%%%%%%%%%%

% This is intentionally left with correspondence to the paper
B =1;
% or
% B = eye(n);
%%%%%%

% Pre - Partition sequence

%%  algorithm
t1=tic;
for iter=1:1:MAX_ITER_Reweight
    
    cum_part = cumsum(partition_tmp);
    f = sqrt(u)/M;
    result.f{iter} = f;
    result.u{iter} = u;
    
    if ~QUIET
        fprintf('Iteration %d: \n ',iter);
    end
    [m_tmp,n_tmp] = size(A);
    Noblock_tmp = length(partition_tmp); % number of blocks
    Doblcok_tmp = round(n_tmp/Noblock_tmp); % dimention of each block
    v_size = Doblcok_tmp*Noblock_tmp;
    u_size = size(u);
    Doblcok_tmp = n_tmp/Noblock_tmp;
    
    % = = = = = = =  = = = = = Solver CVX = = = = = =  = = = = = = = = = = = =
    t2 = tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cvx_begin quiet
    cvx_solver sedumi
    %     cvx_solver SDPT3
    variable v(Doblcok_tmp,Noblock_tmp)
    %     minimize    (lambda*sum(f.*norms( v, 2 ,1 )')+ 0.5*sum_square(A* vec(v)-y) )
    minimize    (2*sum(f.*norms( v, 2 ,1 )')/M+ 1/lambda*sum_square(A* vec(v)-y) )
    %                      subject to
    
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    w = vec(v);
    result.t_cvxwgl1{iter} = toc(t2);
    % = = = = = = =  = = = = = CVX end = = = = = =  = = = = = = = = = = = =
    
    % ===========================Partition===========================
    start_ind = 1;
    for i = 1:length(partition_tmp)
        sel = start_ind:cum_part(i);
        gammai(i,1) = norm(w(sel),2)/(f(i)*M);  %%gammai update !!!!!!!!!!!!
        Gammai(sel,1) = gammai(i,1)*ones(length(sel),1);
        start_ind = cum_part(i) + 1;
    end
    
    % = = = = = = = = = = = = = = = = Prune parameters = = = = = = = = = = = = = = = = = = = = =
    % Prune parameters with low sensitivities
    switch prune
        case 'sensitivity'
            [w, ~]=prunedic(A1, y_d, w, epsilon); % prune the parameters
            Gammai=gammai./u;
            index=find(w);
        case 'gamma'
            result.gammai{iter} = full(gammai);
            result.Gammai{iter} = full(Gammai);
            if strcmp(truncate,'on') == 1
                gammai(gammai<Gammai_threshold) = 0;
                Gammai(Gammai<Gammai_threshold) = 0;
                index_g = find(gammai);
                index_G=find(Gammai);
                partition_tmp = partition(index_g);
                
                gammai = gammai(index_g);
                Gammai = Gammai(index_G);
                result.sparsity_block(iter,1) = length(index_g);
            end
            
    end
    %% ====================Update the re-weights========================
    
    % Prune size and generate new feature lighter feature matrix
    if strcmp(truncate,'on') == 1
        A=A(:,index_G);
        Aa = Aa(:,index_g);
    end
    B = 1/M;
    
    %  For Bw
  
    switch diaginv
        case 'woodbury'
            switch parameter.tmp
                case 'on'
                    VarDiag = woodbury2(A,  Gammai, lambda);
                case 'off'
                    vardiag = woodbury2(Aa,  gammai, lambda);
                    VarDiag = kron(vardiag,eye(L));
            end
            % get the weights
            udiag=diag(A'*VarDiag*A);
        case 'lanzcos'
            % lanczos: Computes diag( B*inv(A)*B' ) with A = X'*R*X + B'*P*B.
            P = diag(1./Gammai);
            VarDiag = diaginv_lanczos(A, lambda^-1,B,P, 50);
            %             VarDiag = diaginv_lanczos(A, lambda^-1,1,P, 50);
            udiag= - VarDiag./Gammai.^2+1./Gammai;
        case 'full'
            % full
            P = diag(1./Gammai);
            VarDiag = diag(B*inv(lambda^-1*A'*A+B'*P*B)*B');
            %             VarDiag = diag(inv(lambda^-1*A'*A+P));
            udiag= - VarDiag./Gammai.^2+1./Gammai;
            
    end
    %% sum of diagonal
    cum_part = cumsum(partition_tmp);
    start_ind = 1;
    Noblock_tmp = length(partition_tmp);
    u = zeros(Noblock_tmp,1);
    for i = 1:Noblock_tmp
        sel = start_ind:cum_part(i);
        u(i,1) = sum(udiag(sel));
        start_ind = cum_part(i) + 1;
    end
    result.uiag{iter} = udiag;
    result.VarDiag{iter} = VarDiag;
    
    % ======= Truncate ==========
    switch truncate
        case 'on'
            Index{iter}=index_G;
            ind1=Index{iter};
            for j=1:iter-1
                ind2=Index{iter-j};
                ind1=insideback(ind2,ind1); % check the function in the folder
            end
            x(INDEX(ind1),iter)=w(index_G);
        case 'off'
            x (:,iter)= w;
    end
    
    result.Estimate(:,iter) = full(Anorm*x(:,iter)/scale);
    result.nmse_y(iter,1)= norm(y - A0*x(:,iter))^2./norm(y)^2;
    %% ====================Report====================
    if ~QUIET
        reportnnz(iter) = length(find(x(:,iter)));
        fprintf('cvxrwgl1: %d nonzeros out of %d.\n',reportnnz(iter),n);
    end
end
%% Result
result.iter = iter;
x = x/scale;
result.x = x;
result.estimate = result.Estimate(:,iter);
result.sparsity = find(result.estimate);
result.elapsedTime_cvxrwgl1=toc(t1);
fprintf('\ncvxrwgl1: elapsedTime is %f seconds.\n',result.elapsedTime_cvxrwgl1);

