function [y, A, x_true, noise, partition]= expandmatrix( y, A, x_true, noise, group)
% clear all
% A_d = randn(3,4);
% x_true = randn(4,1)
% onoff = 'on';
% [A_d, x_true, partition]= expandmatrix( A_d, x_true, onoff)


group_onoff = group.group_onoff;
group_type = group.group_type ;

 
[m,n] = size(A);
partition = ones(n,1);

if strcmp(group_onoff,'on')
    
    if strcmp(group_type,'SMV')
        
        mn = m*n;
        
        for i =1 :n
            AA(:,(i-1)*m+1:i*m) = diag(A(:,i));
            partition(i,1) = m;
        end
        
        A = AA; 
        
        x_true = full(kron(x_true,ones(m,1)));
        
        
    elseif strcmp(group_type,'MMV')
        
        [p,L] = size(y);
        y = vec(y');
        A = kron(A,eye(L));
        x_true = vec(x_true');
        noise = vec(noise');
        partition = L*ones(n,1); 
        
        
    end
    
elseif strcmp(group_onoff,'off')
    
    partition = [];
    return;
    
end



% Random Orthogonalization! funny, maybe some good effect :-)
m = size(A,1);
% tmp = RandOrthMat(m);
tmp = 1;
A = tmp*A;
y = tmp*y;
noise = tmp*noise;
