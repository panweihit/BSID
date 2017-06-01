function dx = example_gene_continuous(t, x, Z, example, dynamicnoise, coefficient)

% clear all;
% close all;
% which_var = 1;
% dynamicnoise = 0;
% % example = 'Lorenz2'
% % example = 'Ecology'
% % example = 'Mackey-Glass'
% example = 'Rossler'

% coefficient.p1 = 40; coefficient.p5 = 1;coefficient.p4 = 0.5;
% example =  'generalrepressilator_continuous'; which_var = 1; dynamicnoise = 0;
% [t,data_matrix] = ode45(@(t,x) example_gene_continuous(t,x,[], example, dynamicnoise, coefficient), [0 10], abs(randn(1, 8)));
% figure; plot(t,y);


switch example
    case 'generalrepressilator_continuous_even'
        %%


        n = 8;
        dx = zeros(n,1);
        noise = dynamicnoise*randn(n,1);
                
        HILL = @(x,num,den)(x.^num./(1+x.^den));
        p1 = coefficient.p1; p4 = coefficient.p4; p5 = coefficient.p5;
        
        Dic_sym1 = '[x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8)]';
        Dic_sym2 = '[HILL(x(1),0,3), HILL(x(2),0,3), HILL(x(3),0,3), HILL(x(4),0,3),HILL(x(5),0,3), HILL(x(6),0,3), HILL(x(7),0,3), HILL(x(8),0,3), 1]';
        Dic_sym3 = '[HILL(x(1),3,3), HILL(x(2),3,3), HILL(x(3),3,3),HILL(x(4),3,3), HILL(x(5),3,3), HILL(x(6),3,3), HILL(x(7),3,3), HILL(x(8),3,3)]';
        %             Dic_sym = '[Dic_sym1, Dic_sym2, Dic_sym3]';
        X_true_mat1 = diag(-p5*ones(8,1));
        X_true_mat2 = [zeros(7,1), diag(p1*ones(7,1)); p1, zeros(1, 7); p4*ones(1, 8)];
        X_true_mat3 = zeros(8,8);
        X_true_mat = [X_true_mat1; X_true_mat2; X_true_mat3];
        
        %[m1,m2] = size(X_true_mat);
        %X_true_mat = (0.5+1.0*rand(m1,m2)).*X_true_mat;
        
        
        % Iterate the model
        
        A_mat1 = eval(Dic_sym1);A_mat2 = eval(Dic_sym2);A_mat3 = eval(Dic_sym3);
        A_mat = [A_mat1, A_mat2, A_mat3];
        for i = 1:1:n
            dx(i) = A_mat*X_true_mat(:,i)+noise(i);
        end
        
    case 'generalrepressilator_continuous_odd'
        %%


        n = 7;
        dx = zeros(n,1);
        noise = dynamicnoise*randn(n,1);
                
        HILL = @(x,num,den)(x.^num./(1+x.^den));
        p1 = coefficient.p1; p4 = coefficient.p4; p5 = coefficient.p5;
        
        Dic_sym1 = '[x(1), x(2), x(3), x(4), x(5), x(6), x(7)]';
        Dic_sym2 = '[HILL(x(1),0,3), HILL(x(2),0,3), HILL(x(3),0,3), HILL(x(4),0,3),HILL(x(5),0,3), HILL(x(6),0,3), HILL(x(7),0,3), 1]';
        Dic_sym3 = '[HILL(x(1),3,3), HILL(x(2),3,3), HILL(x(3),3,3),HILL(x(4),3,3), HILL(x(5),3,3), HILL(x(6),3,3), HILL(x(7),3,3)]';
        %             Dic_sym = '[Dic_sym1, Dic_sym2, Dic_sym3]';
        X_true_mat1 = diag(-p5*ones(7,1));
        X_true_mat2 = [zeros(6,1), diag(p1*ones(6,1)); p1, zeros(1, 6); p4*ones(1, 7)];
        X_true_mat3 = zeros(7,7);
        X_true_mat = [X_true_mat1; X_true_mat2; X_true_mat3];
        
        %[m1,m2] = size(X_true_mat);
        %X_true_mat = (0.5+1.0*rand(m1,m2)).*X_true_mat;
        
        
        % Iterate the model
        
        A_mat1 = eval(Dic_sym1);A_mat2 = eval(Dic_sym2);A_mat3 = eval(Dic_sym3);
        A_mat = [A_mat1, A_mat2, A_mat3];
        for i = 1:1:n
            dx(i) = A_mat*X_true_mat(:,i)+noise(i);
        end
                

end