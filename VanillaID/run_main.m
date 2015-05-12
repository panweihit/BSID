
clear all;
close all;

which = 1;
[y, A, w_true] = GRN_dis(which);

w_lasso = lasso(A,y);
w_lasso = w_lasso(:,end);
w_ls = lscov(A,y); 
lambda = 0.01; MAXITER = 5;
w_ours =  tac_reconstruction(y, A, lambda,MAXITER);
compare = [w_true, w_ours(:,end), w_lasso, w_ls]
