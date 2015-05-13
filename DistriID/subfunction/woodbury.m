 function Abarinv = woodbury(A, f, rho)
%     Formula:     inv(A*R*A' + P) = inv(P) - inv(P)*A*inv( inv(R)+A'*inv(P)*A )*A'*inv(P)
%                             where P = rho*f.^2; R =1
%     Abar=A'*A + rho*sparse(diag(f.*f)); 
%     Abarinv=inv(Abar);

% clear all
%      A = randn(3,4);
%      rho =1;
%      u = abs(randn(4,1));

p = 1e-12;

[m,n] = size(A);

if m>n
    AtA = A'*A;
    AtA(1:n+1:end) = diag(AtA) + rho*f.^2;
    ep =  p*max(AtA(:)); 
    Abarinv=inv(AtA+ep*eye(n));
%     AtA = A'*A;
%     AtA(1:n+1:end) = diag(AtA)/rho + u.^2;
%     ep =  p*max(AtA(:)); 
%     Abarinv=inv(AtA+ep*eye(n))/rho;
else
%      inv(A'*A + P) = Pinv - Pinv*A'*inv(eye(m)+A*Pinv*A' )*A*Pinv;
%      fprintf(' woodbury \n');
     pinv = 1./(rho*f.^2);
     Pinv = sparse(diag(pinv));
     P_tmp = eye(m)+A*Pinv*A';
     ep =  p*max(P_tmp(:)); 
%      iP = inv( P_tmp +ep*eye(m)); % stabilise inversion of P
     iP = P_tmp +ep*eye(m);
     tmp2 = A'*(iP\A);

%      Abarinv = Pinv - Pinv*tmp2*Pinv;
     Abarinv = Pinv - (pinv*pinv').*tmp2;

  
%     tic
%      Abarinv2 = Pinv - Pinv*tmp2*Pinv;
%     toc 
end
    
