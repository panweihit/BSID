function [y_d,A_d] = diffdata(y, A, noise, diff_onoff)
%% plot the oringnal and transformed data


onoff= diff_onoff;

QUIET = 1;
switch onoff
    case 'on'
    
%     if length(y)>15
%         if adftest(y) == 1
%             y_d = y;
%             A_d = A;
%             noise_d = noise;
%         else
%             y_d = diff(y);
%             A_d = diff(A);
%             noise_d = diff(noise);
%         end
%     else
%         y_d = diff(y);
%         A_d = diff(A);
%         noise_d = diff(noise);
%     end
        y_d = diff(y);
        A_d = diff(A);
        noise_d = diff(noise);
    case 'off'
    
    y_d = y;
    A_d = A;
    noise_d = noise;
    
end

% if norm(y_d-A_d*x_true-noise_d)>1e-4
%     error('transformation is wrong');
% end


% Random Orthogonalization! funny, maybe some good effect :-)
m = size(A_d,1);
% tmp = RandOrthMat(m);
tmp = 1;
A_d = tmp*A_d;
y_d = tmp*y_d;
noise_d = tmp*noise_d;

% index_nnz = find(norms(A_d,[],1))
% A_d = A_d(:,index_nnz);
% x_true = x_true(index_nnz,:);

if ~QUIET
    figure;
    plot(y,'k');hold on;plot(y_d,'r-o');hold on;legend('oringnal','diff');
end
% [m,n] = size(A);
% Anorm_tmp = spdiags(1./sqrt(sum(A.^2))',0,n,n);
% A = A*Anorm_tmp;
% x_true = x_true./diag(Anorm_tmp);
