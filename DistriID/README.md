This folder implements the following paper

W. Pan, A. Sootla, and G.-B. Stan. Distributed Reconstruction of Nonlinear Networks: An ADMM Approach. The International Federation of Automatic Control. Cape Town, South Africa, 2014.

==================== Quick Guide to Start ===========================

1. Download CVX toolbox from www.cvxr.com and run 'cvx_setup.m' to install the toolbox.
2. Run 'start_up.m' to the toolbox to your directory
3. RUN 'RUN_DistriID.m' and have fun!

============ Guidience on the Use of Parallel Computing Toolbox in Matlab =============== 

In Matlab 2015a, use     
parpool(numberofcluster); or delete(gcp);
to start or delete multicluster

In previous version, check it on MathWorks Website, typically, use matlabpool(numberofcluster)

=============== Acknowledgement ===============

1. The four files in folder 'inf' are implenmted by Dr Hannes Nickisch
http://hannes.nickisch.org/
http://hannes.nickisch.org/code/glm-ie/doc/index.html

2. dwl1.m and wl1.m are implmemetated based on 
http://web.stanford.edu/~boyd/papers/admm/lasso/lasso.html
and 
http://web.stanford.edu/~boyd/papers/admm/group_lasso/group_lasso_feat_split.html



