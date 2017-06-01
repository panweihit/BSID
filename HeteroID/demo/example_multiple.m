function data = example_multiple(select, parameter, which_var)
group.group_onoff = parameter.group_onoff;
group.group_type = parameter.group_type;
dynamicnoisestd = parameter.dynamicnoise;


switch select
    case 'generalrepressilator_discrete_even'
        %%
        
%         clear all; close all
%         parameter.multiple = 'on';
%         parameter.single = 'single';
%         parameter.replicate = 5;
%         parameter.plot_onoff = 'on';
%         group.group_onoff = 'on';
%         which_var = 1;
%         dynamicnoisestd = 0e-1;
%         parameter.measurementnoiseperc = 0.0;
%         parameter.derivative = 'eular'; % change this accordingly
%         %parameter.derivativekth = 4;
%         parameter.T = 500;

        t_start = tic;
        n_var = 8;
        start = 1;
        HILL = @(x,num,den)(x.^num./(1+x.^den));
        dt = .01;
        T = parameter.T;
        CUT= round(1*kron(T,ones(1, parameter.replicate))); % interval of interest and selected by you
        
        noise = dynamicnoisestd*randn(T,8);
        for replicate = 1:parameter.replicate
            p1 = 20; p5 = 2; p4 = .5;
            coefficient.p1 = (0.8+1.2*rand(1))*p1;
            coefficient.p4 = (0.8+1.2*rand(1))*p4;
            coefficient.p5 = (0.8+1.2*rand(1))*p5;
            p1 = coefficient.p1; p4 = coefficient.p4; p5 = coefficient.p5;
            initial = abs(randn(1, n_var));
            %             initial = 10*[0.3,0.04,0.1,0.4,0.3,0.3,0.3,1.0];
            
            x(1,:) = initial;
            Initial(replicate,:) = initial;
            Dic_sym1 = '[x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), x(k,7), x(k,8)]';
            Dic_sym2 = '[HILL(x(k,1),0,3), HILL(x(k,2),0,3), HILL(x(k,3),0,3), HILL(x(k,4),0,3),HILL(x(k,5),0,3), HILL(x(k,6),0,3), HILL(x(k,7),0,3), HILL(x(k,8),0,3)]';
            Dic_sym3 = '[HILL(x(k,1),3,3), HILL(x(k,2),3,3), HILL(x(k,3),3,3),HILL(x(k,4),3,3), HILL(x(k,5),3,3), HILL(x(k,6),3,3), HILL(x(k,7),3,3), HILL(x(k,8),3,3), 1]';
            %             Dic_sym = '[Dic_sym1, Dic_sym2, Dic_sym3]';
            X_true_mat1 = diag(-p5*ones(n_var,1));
            X_true_mat2 = [zeros(7,1), diag(p1*ones(7,1)); p1, zeros(1, 7)];
            X_true_mat3 = zeros(n_var,n_var);
            X_true_mat = [X_true_mat1; X_true_mat2; X_true_mat3; p4*ones(1, 8)];
                       
            % Iterate the model
            for k = 1:T
                A_mat1(k,:) = eval(Dic_sym1);A_mat2(k,:) = eval(Dic_sym2);A_mat3(k,:) = eval(Dic_sym3);
                A_mat(k,:) = [A_mat1(k,:), A_mat2(k,:), A_mat3(k,:)];
                for i = 1:1:n_var
                    x(k+1,i) = x(k, i) + dt*A_mat(k,:)*X_true_mat(:,i)+sqrt(dt)*noise(k,i);
                end
            end
            
            x = addwhite(x, parameter.measurementnoiseperc);
            t= 1:size(x,1);
            
            cut_i = CUT(replicate);
            Data_mat = x(start:start+cut_i,:);
            t_mat = dt*t(start:start+cut_i);
            A_mat = A_mat(start:start+cut_i-1,:);
            Y_mat = (Data_mat(2:end,:) - Data_mat(1:end-1,:))/dt;
            Y{replicate} = Y_mat(:,which_var);
            Dic{replicate} = A_mat;
            W_true{replicate} = X_true_mat(:,which_var);
            
            Dic{replicate} =   Dic{replicate}(4:end,:); % ignore the first few estimate
            Y{replicate} = Y{replicate}(4:end,:);
            DY{replicate} = Dic{replicate}*W_true{replicate};
            error_norm(replicate) = norm(Y{replicate} - DY{replicate});
            
            if strcmp(parameter.plot_onoff,  'on')
                figure;
                plot(t, x); hold on;
                plot(t, x(:, which_var),'-*'); hold on;
                xlabel('time'); ylabel('quantity')
                legend('the selected gene for identification');
                title('State Dynamics')
%                 figure;
%                 plot([Y{replicate}]);hold on;
%                 plot(DY{replicate},'r');
%                 legend('approximate','true')
%                 title('Derivative Dynamics')
            end
        end
        
        error_norm
        
        if strcmp(parameter.multiple,'on')
            [y,A, x_true] = multidic(Dic, Y, W_true);
            N = length(W_true{1});
            error_norm = norm(y-A*x_true)
            partition = kron(replicate, ones(N,1));
        elseif strcmp(parameter.multiple,'off')
            if strcmp(parameter.single,'on')
                y = Y{replicate};
                A = Dic{replicate};
                x_true = W_true{replicate};
                partition =[];
            elseif strcmp(parameter.single,'off')
                y = [];
                A = [];
                for j =1: 1:replicate
                    y = vertcat(Y{:});
                    A =vertcat(Dic{:});
                    x_true = W_true{replicate};
                    partition =[];
                end
            end
            error_norm = norm(y-A*x_true)
        end
        
        
    case 'generalrepressilator_discrete_odd'
        %%
%         
%         clear all; close all
%         parameter.multiple = 'on';
%         parameter.single = 'single';
%         parameter.replicate = 5;
%         parameter.plot_onoff = 'on';
%         group.group_onoff = 'on';
%         which_var = 1;
%         dynamicnoisestd = 0e-1;
%         parameter.measurementnoiseperc = 0.0;
%         parameter.derivative = 'eular'; % change this accordingly
%         %parameter.derivativekth = 4;
%         parameter.T = 500;
%         
        t_start = tic;
        n_var = 7;
        start = 1;
        HILL = @(x,num,den)(x.^num./(1+x.^den));
        dt = .01;
        T = parameter.T;
        CUT= round(1*kron(T,ones(1, parameter.replicate))); % interval of interest and selected by you
        
        noise = dynamicnoisestd*randn(T,n_var);
        for replicate = 1:parameter.replicate
            p1 = 20; p5 = 2; p4 = .5;
            coefficient.p1 = (0.8+1.2*rand(1))*p1;
            coefficient.p4 = (0.8+1.2*rand(1))*p4;
            coefficient.p5 = (0.8+1.2*rand(1))*p5;
            p1 = coefficient.p1; p4 = coefficient.p4; p5 = coefficient.p5;
            initial = abs(randn(1, n_var));
            %             initial = 10*[0.3,0.04,0.1,0.4,0.3,0.3,0.3,1.0];
            
            x(1,:) = initial;
            Initial(replicate,:) = initial;
            Dic_sym1 = '[x(k,1), x(k,2), x(k,3), x(k,4), x(k,5), x(k,6), x(k,7)]';
            Dic_sym2 = '[HILL(x(k,1),0,3), HILL(x(k,2),0,3), HILL(x(k,3),0,3), HILL(x(k,4),0,3),HILL(x(k,5),0,3), HILL(x(k,6),0,3), HILL(x(k,7),0,3)]';
            Dic_sym3 = '[HILL(x(k,1),3,3), HILL(x(k,2),3,3), HILL(x(k,3),3,3),HILL(x(k,4),3,3), HILL(x(k,5),3,3), HILL(x(k,6),3,3), HILL(x(k,7),3,3), 1]';
            %             Dic_sym = '[Dic_sym1, Dic_sym2, Dic_sym3]';
            X_true_mat1 = diag(-p5*ones(n_var,1));
            X_true_mat2 = [zeros(6,1), diag(p1*ones(6,1)); p1, zeros(1, 6)];
            X_true_mat3 = zeros(n_var,n_var);
            X_true_mat = [X_true_mat1; X_true_mat2; X_true_mat3; p4*ones(1, 7)];
            
            % Iterate the model
            for k = 1:T
                A_mat1(k,:) = eval(Dic_sym1);A_mat2(k,:) = eval(Dic_sym2);A_mat3(k,:) = eval(Dic_sym3);
                A_mat(k,:) = [A_mat1(k,:), A_mat2(k,:), A_mat3(k,:)];
                for i = 1:1:n_var
                    x(k+1,i) = x(k, i) + dt*A_mat(k,:)*X_true_mat(:,i)+sqrt(dt)*noise(k,i);
                end
            end
            
            x = addwhite(x, parameter.measurementnoiseperc);
            t= 1:size(x,1);
            
            cut_i = CUT(replicate);
            Data_mat = x(start:start+cut_i,:);
            t_mat = dt*t(start:start+cut_i);
            A_mat = A_mat(start:start+cut_i-1,:);
            Y_mat = (Data_mat(2:end,:) - Data_mat(1:end-1,:))/dt;
            Y{replicate} = Y_mat(:,which_var);
            Dic{replicate} = A_mat;
            W_true{replicate} = X_true_mat(:,which_var);
            
            Dic{replicate} =   Dic{replicate}(4:end,:); % ignore the first few estimate
            Y{replicate} = Y{replicate}(4:end,:);
            DY{replicate} = Dic{replicate}*W_true{replicate};
            error_norm(replicate) = norm(Y{replicate} - DY{replicate});
            
            if strcmp(parameter.plot_onoff,  'on')
                figure;
                plot(t, x); hold on;
                plot(t, x(:, which_var),'-*'); hold on;
                xlabel('time'); ylabel('quantity')
                legend('the selected gene for identification');
                title('State Dynamics')
%                 figure;
%                 plot([Y{replicate}]);hold on;
%                 plot(DY{replicate},'r');
%                 legend('approximate','true')
%                 title('Derivative Dynamics')
            end
        end
        
        error_norm
        
        if strcmp(parameter.multiple,'on')
            [y,A, x_true] = multidic(Dic, Y, W_true);
            N = length(W_true{1});
            error_norm = norm(y-A*x_true)
            partition = kron(replicate, ones(N,1));
        elseif strcmp(parameter.multiple,'off')
            if strcmp(parameter.single,'on')
                y = Y{replicate};
                A = Dic{replicate};
                x_true = W_true{replicate};
                partition =[];
            elseif strcmp(parameter.single,'off')
                y = [];
                A = [];
                for j =1: 1:replicate
                    y = vertcat(Y{:});
                    A =vertcat(Dic{:});
                    x_true = W_true{replicate};
                    partition =[];
                end
            end
            error_norm = norm(y-A*x_true)
        end
        
        toc(t_start)
    case 'generalrepressilator_continuous_even'
        %%
        
%         clear all; close all
%         select = 'generalrepressilator_continuous_even';
%         parameter.dynamicnoise = 0;
%         parameter.multiple = 'on';
%         parameter.single = 'single';
%         parameter.replicate = 1;
%         parameter.plot_onoff = 'on';
%         group.group_onoff = 'on';
%         which_var = 2;
%         parameter.measurementnoiseperc = 50;
%         parameter.prop = 10/10;
%         dynamicnoisestd = parameter.dynamicnoise;
        
        t_start = tic;
        n_var = 8;
        HILL = @(x,num,den)(x.^num./(1+x.^den));
        DT = 0.1*ones(parameter.replicate,1);
        T = 20;
        initial = abs(10*randn(1, n_var));
        
        for replicate = 1:parameter.replicate
            dt = DT(replicate);
            p1 = 40; p5 = 1; p4 = .5;
            coefficient.p1 = (0.2+1.0*rand(1))*p1;
            coefficient.p4 = (0.2+1.0*rand(1))*p4;
            coefficient.p5 = (0.02+1.0*rand(1))*p5;
            p1 = coefficient.p1; p4 = coefficient.p4; p5 = coefficient.p5;
            
            [t, data_matrix] = ode45(@(t,x) example_gene_continuous...
                (t,x,[], select, dynamicnoisestd, coefficient), 0:dt:T, initial);
            
            x_full = data_matrix;
            x_full_measnoise = addwhite(x_full, parameter.measurementnoiseperc);
            x = x_full_measnoise(1:1:end,:);
            
            Dic_sym1 = '[x(:,1), x(:,2), x(:,3), x(:,4), x(:,5), x(:,6), x(:,7), x(:,8)]';
            Dic_sym2 = '[HILL(x(:,1),0,3), HILL(x(:,2),0,3), HILL(x(:,3),0,3), HILL(x(:,4),0,3),HILL(x(:,5),0,3), HILL(x(:,6),0,3), HILL(x(:,7),0,3), HILL(x(:,8),0,3)]';
            Dic_sym3 = '[HILL(x(:,1),3,3), HILL(x(:,2),3,3), HILL(x(:,3),3,3),HILL(x(:,4),3,3), HILL(x(:,5),3,3), HILL(x(:,6),3,3), HILL(x(:,7),3,3), HILL(x(:,8),3,3), ones(size(x,1),1)]';
            %            Dic_sym = '[Dic_sym1, Dic_sym2, Dic_sym3]';
            X_true_mat1 = diag(-p5*ones(n_var,1));
            X_true_mat2 = [zeros(7,1), diag(p1*ones(7,1)); p1, zeros(1, 7)];
            X_true_mat3 = zeros(n_var,n_var);
            X_true_mat = [X_true_mat1; X_true_mat2; X_true_mat3; p4*ones(1, 8)];
            
            A_mat1 = eval(Dic_sym1);A_mat2 = eval(Dic_sym2);A_mat3 = eval(Dic_sym3);
            A_mat = [A_mat1, A_mat2, A_mat3];
            
            start = 1;
            cut_i = round(parameter.prop*(size(x,1)-start));
            cut(replicate) = cut_i;
            Data_mat = x(start:start+cut_i,:);
            t_mat = t(start:start+cut_i);
            A_mat = A_mat(start:start+cut_i-1,:);
            
            % estimate derivative
            [data_new, derivative, index] = estimatediff(Data_mat, t_mat, 'solver', 1, []);
            
            Y{replicate} = derivative(:,which_var);
            Dic{replicate} = A_mat(index, :);
            W_true{replicate} = X_true_mat(:,which_var);
            
            Dic{replicate} = Dic{replicate} (10:end,:); % ignore the first 10 estimate,
            Y{replicate} = Y{replicate}(10:end,:);
            DY{replicate} = Dic{replicate}*W_true{replicate};
            error_norm(replicate) = norm(Y{replicate} - DY{replicate});
            simulate_error{replicate} = Y{replicate} - DY{replicate};
            
            if strcmp(parameter.plot_onoff,  'on')
                figure;
                plot(t, x); hold on;
                plot(t, x(:, which_var),'-*'); hold on;
                xlabel('time'); ylabel('quantity')
                legend('the selected gene for identification');
                title('State Dynamics')
                figure;
                plot([Y{replicate}]);hold on;
                plot(DY{replicate},'r');
                legend('approximate','true')
                title('Derivative Dynamics')
            end
        end
        
        error_norm
        
        if strcmp(parameter.multiple,'on')
            [y,A, x_true] = multidic(Dic, Y, W_true);
            N = length(W_true{1});
            err_dic = norm(y-A*x_true);
            partition = kron(replicate, ones(N,1));
        elseif strcmp(parameter.multiple,'off')
            if strcmp(parameter.single,'on')
                y = Y{replicate};
                A = Dic{replicate};
                x_true = W_true{replicate};
                partition =[];
            elseif strcmp(parameter.single,'off')
                y = [];
                A = [];
                for j =1: 1:replicate
                    y = vertcat(Y{:});
                    A =vertcat(Dic{:});
                    x_true = W_true{replicate};
                    partition =[];
                end
            end
            error_norm = norm(y-A*x_true)
        end
        
        t_ode45 = toc(t_start)
    case 'generalrepressilator_continuous_odd'
        %%
        
%         clear all; close all
%         select = 'generalrepressilator_continuous_odd';
%         parameter.dynamicnoise = 0;
%         parameter.multiple = 'on';
%         parameter.single = 'single';
%         parameter.replicate = 1;
%         parameter.plot_onoff = 'on';
%         group.group_onoff = 'on';
%         which_var = 1;
%         parameter.measurementnoiseperc = 0;
%         parameter.prop = 10/10;
%         dynamicnoisestd = parameter.dynamicnoise;
        
        t_start = tic;
        n_var = 7;
        HILL = @(x,num,den)(x.^num./(1+x.^den));
        %         DT = [.1,.01, 0.001];
        DT = 0.1*ones(parameter.replicate,1);
        T = 50;
        initial = abs(10*randn(1, n_var));
        
        for replicate = 1:parameter.replicate
            dt = DT(replicate);
            p1 = 40; p5 = 1; p4 = .5;
            coefficient.p1 = (0.2+1.0*rand(1))*p1;
            coefficient.p4 = (0.2+1.0*rand(1))*p4;
            coefficient.p5 = (0.02+1.0*rand(1))*p5;
            p1 = coefficient.p1; p4 = coefficient.p4; p5 = coefficient.p5;
            
            %             initial = 10*[0.3,0.04,0.1,0.4,0.3,0.3,0.3,1.0];
            [t, data_matrix] = ode45(@(t,x) example_gene_continuous...
                (t,x,[], select, dynamicnoisestd, coefficient), 0:dt:T, initial);
            
            x_full = data_matrix;
            x_full_measnoise = addwhite(x_full, parameter.measurementnoiseperc);
            x = x_full_measnoise(1:1:end,:);
            
            Dic_sym1 = '[x(:,1), x(:,2), x(:,3), x(:,4), x(:,5), x(:,6), x(:,7)]';
            Dic_sym2 = '[HILL(x(:,1),0,3), HILL(x(:,2),0,3), HILL(x(:,3),0,3), HILL(x(:,4),0,3),HILL(x(:,5),0,3), HILL(x(:,6),0,3), HILL(x(:,7),0,3)]';
            Dic_sym3 = '[HILL(x(:,1),3,3), HILL(x(:,2),3,3), HILL(x(:,3),3,3),HILL(x(:,4),3,3), HILL(x(:,5),3,3), HILL(x(:,6),3,3), HILL(x(:,7),3,3), ones(size(x,1),1)]';
            X_true_mat1 = diag(-p5*ones(n_var ,1));
            X_true_mat2 = [zeros(6,1), diag(p1*ones(6,1)); p1, zeros(1, 6)];
            X_true_mat3 = zeros(n_var ,n_var );
            X_true_mat = [X_true_mat1; X_true_mat2; X_true_mat3; p4*ones(1, 7)];
            
            A_mat1 = eval(Dic_sym1);A_mat2 = eval(Dic_sym2);A_mat3 = eval(Dic_sym3);
            A_mat = [A_mat1, A_mat2, A_mat3];
            
            start = 1;
            cut_i = round(parameter.prop*(size(x,1)-start));
            cut(replicate) = cut_i;
            Data_mat = x(start:start+cut_i,:);
            t_mat = t(start:start+cut_i);
            A_mat = A_mat(start:start+cut_i-1,:);
            
            % estimate derivative
            [data_new, derivative, index] = estimatediff(Data_mat, t_mat, 'solver', 1, []);
            
            Y{replicate} = derivative(:,which_var);
            Dic{replicate} = A_mat(index, :);
            W_true{replicate} = X_true_mat(:,which_var);
            
            Dic{replicate} = Dic{replicate} (10:end,:); % ignore the first 10 estimate,
            Y{replicate} = Y{replicate}(10:end,:);
            DY{replicate} = Dic{replicate}*W_true{replicate};
            error_norm(replicate) = norm(Y{replicate} - DY{replicate});
            simulate_error{replicate} = Y{replicate} - DY{replicate};
            
            if strcmp(parameter.plot_onoff,  'on')
                figure;
                plot(t, x); hold on;
                plot(t, x(:, which_var),'-*'); hold on;
                xlabel('time'); ylabel('quantity')
                legend('the selected gene for identification');
                title('State Dynamics')
                figure;
                plot([Y{replicate}]);hold on;
                plot(DY{replicate},'r');
                legend('approximate','true')
                title('Derivative Dynamics')
            end
        end
        
        error_norm
        
        if strcmp(parameter.multiple,'on')
            [y,A, x_true] = multidic(Dic, Y, W_true);
            N = length(W_true{1});
            err_dic = norm(y-A*x_true);
            partition = kron(replicate, ones(N,1));
        elseif strcmp(parameter.multiple,'off')
            if strcmp(parameter.single,'on')
                y = Y{replicate};
                A = Dic{replicate};
                x_true = W_true{replicate};
                partition =[];
            elseif strcmp(parameter.single,'off')
                y = [];
                A = [];
                for j =1: 1:replicate
                    y = vertcat(Y{:});
                    A =vertcat(Dic{:});
                    x_true = W_true{replicate};
                    partition =[];
                end
            end
            error_norm = norm(y-A*x_true)
        end
        
        t_ode45 = toc(t_start)
end

%%
if strcmp(group.group_onoff, 'on')
    cum_part = cumsum(partition);
    cum_part = [0;cum_part];
    for i = 1:length(partition)
        grouping{i} = [cum_part(i)+1:cum_part(i+1)]';
        index_nnz{i} = find(norms(A(:,grouping{i}),[],1));
        partition_new(i) = length(index_nnz{i});
    end
else
    grouping = [];
end

index_partition = find(partition_new);
partition = partition_new(index_partition);
scale = 1;
y = scale*y;
index_nnz = find(norms(A,[],1));
A = A(:,index_nnz);
x_true = scale*x_true(index_nnz,:);

data.partition = partition;
data.grouping = grouping;
data.y = y;
data.Y = Y;
data.Dic = Dic;
data.A = A;
data.x_true = x_true;

