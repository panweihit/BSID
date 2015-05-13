 function Dic0inv = woodbury2(A, Gammai, lambda)
%     Formula:     inv(A*R*A' + P) = inv(P) - inv(P)*A*inv( inv(R)+A'*inv(P)*A )*A'*inv(P)
%
%     Dic0inv=inv(lambda*eye(m)+A*diag(G)*A');

    [m_A,n_A] = size(A);
    [m_G,n_G] = size(Gammai);
    
    p = 1e-12;
    
    AG = (A* sparse(diag(sqrt(Gammai))))';
    % check dimension use woodbury identity
if m_A<n_A
    AGtAG = AG'*AG;
    AGtAG(1:m_A+1:end) = diag(AGtAG) + lambda*ones(m_A,1);
    ep =  p*max(AGtAG(:)); 
    Dic0inv=inv(AGtAG+ep*eye(m_A));
%     AGtAG = AG'*AG;
%     AGtAG(1:m+1:end) = diag(AGtAG)/lambda + ones(m,1);
%     ep =  p*max(AGtAG(:)); 
%     Dic0inv=inv(AGtAG+ep*eye(m))/lambda;
else
%      inv(A'*A + P) = Pinv - Pinv*A'*inv(eye(m_tmp)+A*Pinv*A' )*A*Pinv;
     DicPinv = diag(1./lambda*ones(m_A,1));
     m_A;
     n_A;
     DicP_tmp = eye(n_A)+AG*DicPinv*AG';
     ep = p*max(DicP_tmp(:));    
     DiciP =(DicP_tmp+ep*eye(n_A))\eye(n_A);
     DicPinv_tmp = AG*DicPinv;
     Dic0inv = DicPinv - DicPinv_tmp'*DiciP*DicPinv_tmp;
end
