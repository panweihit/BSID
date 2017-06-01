clear all

% Approximates diag( B*inv(A)*B' ) with A = X'*R*X + B'*P*B through
 
X = rand(3,5);
[m,n] = size(X);
for i = 1:m
    XX{1,i} = X(i,:);
end
Xblk=blkdiag(XX{:});
XblkT = Xblk'*Xblk

B = eye(n);
b = speye(n);
R = eye(m);
P = eye(n);
P1 =abs(ones(1, n));
P2 = repmat(P1,n,1);
% B = P2.*B;
% P =1;

A = X'*R*X + B'*P*B;
tic;c =diag( B*inv(A)*B' );toc
tic;[c1,cu1,ldA,Q,T] = diaginv_lanczos(X,R,B,P,1001);toc
invA1 = Q*inv(T)*Q';
fsd1 = invA1*A;
tic;[c2,cu2] = diaginv_woodbury(X,R,P); toc
tic; [c3,cu3] = diaginv_sample(X,R,B,P,100,500); toc
% tic;[c3,cu3,~,~] =  diaginv_factorial(X,R,B,P); toc
% tic;c4 = diaginv_full(X,R,B,P); toc
% tic;c5 = linsolve_lcg(X,R,B,P,b); toc