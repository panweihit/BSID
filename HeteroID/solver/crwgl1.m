function result = crwgl1(A, y,  partition, parameter)
fprintf('This is crwgl1 \n');
% Centralised Reweighted Group Lasso

% global parameter
QUIET = 0; % plot or not

%% function parameters

group.group_onoff = parameter.group_onoff;
group.group_type = parameter.group_type ;
L = size(y,2);
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
scale = parameter.scale;
scale_w = parameter.scale_w;
MAX_ITER_Reweight = parameter.MAX_ITER_Reweight; % re-weighted iteration number
Gammai_threshold = parameter.Gammai_threshold;
truncate = parameter.truncate;
epsilon = parameter.epsilon;
Anorm_onoff = parameter.Anorm_onoff;

%% Parameter Initialisation

[m,n]=size(A);
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
        y = y/norm(y);
    case 'off'
        Anorm = 1;
        A0 = A*Anorm; % initializae featrue matrix for the following
        A = A0;
end

switch parameter.tmp
    case 'on'
        tmp = 100;%randn(m,m);
    case 'off'
        tmp = 1;
end

A = tmp*scale*A;
y= tmp*scale*y;

y = scale_w*y;


%%%%%%%%%%%
% u=abs(randn(p,1)); 
u=ones(p,1); %
B =1; 
%%%%%%%%%%%

% Pre - Partition sequence

%%  algorithm

t1=tic;

cum_part = cumsum(partition_tmp);

for iter=1:1:MAX_ITER_Reweight

    lengthpartition = length(partition_tmp);
    
    f = sqrt(abs(u));
    result.f{iter} = f;
    result.D{iter} = f;
    result.u{iter} = u;
    
    if ~QUIET
        fprintf('Iteration %d / total %d iterations.\n', iter, MAX_ITER_Reweight);
    end
       
    if strcmp(parameter.solver,'admm')
        % =================Solver: Weighted l1-minimization==============
        % Cache A Ainv U wl1
        % check dimension use woodbury identity
        t2 = tic;
        [w, history] = wgl1(A,  y, f, partition_tmp, parameter);
              
        result.w{iter}  = w;
        result.history{iter}  = history;
        result.rho{iter} = history.rho;
        %     result.F{iter} = history.F;
        result.t_wgl1{iter} = toc(t2);
        % =======================Solver End==========================
        
    elseif strcmp(parameter.solver,'cvx')
        
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
        minimize    (lambda*sum(f.*norms( v, 2 ,1 )')+ 1/2*sum_square(A* vec(v)-y) )
        %                      subject to
        
        cvx_end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        w = vec(v);
        result.w{iter}  = w;
        result.t_cvxwgl1{iter} = toc(t2);
        % = = = = = = =  = = = = = CVX end = = = = = =  = = = = = = = = = = = =
    end
    
    
    % ===========================Partition===========================
    
    start_ind = 1;
    for i = 1:length(partition_tmp)
        sel = start_ind:cum_part(i);   %norm(w(sel1),2)/(fx(i)*M)
        gammai(i,1) = norm(w(sel),2)/(sqrt(u(i))*length(sel));  %%gammai update !!!!!!!!!!!!
        Gammai(sel,1) = kron(gammai(i,1), ones(length(sel),1));
        start_ind = cum_part(i) + 1;
    end
        
    % = = = = = = = = = = = = = = = = Prune parameters = = = = = = = = = = = = = = = = = = = = =
        
    result.gammai{iter} = gammai;
    result.Gammai{iter} = Gammai;
        
    if strcmp(truncate,'on') == 1
        gammai(gammai<Gammai_threshold(iter)) = 0;
        Gammai(Gammai<Gammai_threshold(iter)) = 0;
        index_g = find(gammai);
        index_G = find(Gammai);
        partition_tmp = partition_tmp(index_g);
        gammai = gammai(index_g);
        Gammai = Gammai(index_G);
        result.sparsity_block(iter,1) = length(index_g);
    end
    
    
    %% ====================Update the re-weights========================
    
    % Prune size and generate new feature lighter feature matrix
    
    if strcmp(truncate,'on') == 1
        A=A(:,index_G);
        Aa = Aa(:,index_g);
        group.Aa = Aa;
    end
    
    %  For Bw
    switch diaginv
        case 'woodbury'
            switch parameter.tmp
                case 'on'
                    VarDiag = woodbury2(A,  Gammai, lambda);
                case 'off'
                    sizegammai = size(gammai);
                    VarDiag = woodbury2(A,  Gammai, lambda);
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
    %     sizeudiag =  size(udiag)
    %     partition_tmp = partition_tmp
    cum_part = cumsum(partition_tmp);
    start_ind = 1;
    p_tmp = length(partition_tmp);
    u = zeros(p_tmp,1);
    for i = 1:p_tmp
        sel = start_ind:cum_part(i);
        u(i,1) = sum(udiag(sel))/length(sel);
        start_ind = cum_part(i) + 1;
    end
    
    result.uiag{iter} = udiag;
    result.VarDiag{iter} = VarDiag;
    
    % ======= Truncate ==========
    
    if strcmp(truncate,'on') == 1     
        Index{iter}=index_G;       
        ind1=Index{iter};
        for j=1:iter-1
            ind2=Index{iter-j};
            ind1=insideback(ind2,ind1); % check the function in the folder
        end
        
        x(INDEX(ind1),iter)=w(index_G);
        
    elseif strcmp(truncate,'off') == 1
        x (:,iter)= w;
    end
    
    result.Estimate(:,iter) = sparse(Anorm*x(:,iter))/scale_w;    
    result.nmse_y(iter,1)= norm(y - A0*x(:,iter))^2./norm(y)^2;
    
    %% ====================Report====================
    if ~QUIET
        %         reportnnz(iter) = length(find(x(:,iter)));
        reportnnz(iter) = length(find(gammai));
        fprintf('crwgl1 : %d nonzeros out of %d.\n',reportnnz(iter),p);
    end
    
end

%% Result
result.iter = iter;
x = sparse(x);
result.x = x;
result.estimate = result.Estimate(:,iter);
result.sparsity = find(result.estimate);
result.elapsedTime_crwgl1 =toc(t1);

cum_part = cumsum(partition);
cum_part = [0;cum_part(:)];
for i = 1:length(partition)
    grouping{i} = [cum_part(i)+1:cum_part(i+1)]';
    index_nnz(i) = norms(result.estimate(grouping{i}));
end
result.sparsepattern = find(index_nnz);

fprintf('\ncrwgl1: Elapsed time is  %f seconds.\n', result.elapsedTime_crwgl1);



