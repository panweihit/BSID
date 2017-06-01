function [W, threshold_perc, threshold_abs]=prunedic(A, y, W,epsilon, group, partition_tmp)

% clear all
% A = randn(3,4);
% W = [1 0.0001 0 0 ]';
% y = A*W;
% epsilon = 1e-4;
% group.group_onoff='on';
% partition_tmp = [1,1,1,1];

 
if strcmp(group.group_onoff, 'off')
    Anorm = eye(size(A,2));
% %Q measure
%     S=A*diag(W);
%     AtS=abs(diag(A'*S));
%     Aty=abs(A'*y);
%     threshold=normc(AtS./Aty).^1
%     threshold(threshold<epsilon)=0;
%     index_s=find(~threshold);
%     W(index_s)=0;

% P measure     
    Anorminv=diag(1./diag(Anorm));
    A1=A*Anorminv;
    W1=Anorm*W;
    
    S=A1*diag(W1);
    AtS=abs(diag(A1'*S));
    Aty=abs(A1'*y);
    threshold_abs = AtS./Aty;
    threshold=normc(AtS./Aty).^2;
    threshold_perc=threshold*100;%percentage
%     threshold(threshold<epsilon)=0;
    index_s=find(~threshold);
    W(index_s)=0;
    

    % % row normarlization
    % Anorm1 = spdiags(1./sqrt(sum(A1.^2,2)')',0,m,m); % normalize columns
    % A1=Anorm1*A1;
    % x_true1 = x_true1;
    % y1=Anorm1*y1;
elseif strcmp(group.group_onoff, 'on')
    
    cum_part = cumsum(partition_tmp);
    start_ind = 1;
    for i = 1:length(partition_tmp)
        sel = start_ind:cum_part(i);   %norm(w(sel1),2)/(fx(i)*M)
        threshold_tmp(i) = norm(A(:,sel)'*A(:,sel)*W(sel))/norm(A(:,sel)'*y);
        start_ind = cum_part(i) + 1;
    end
    threshold_abs = threshold_tmp;
    threshold = threshold_tmp.^2/norm(threshold_tmp)^2;
    threshold_perc=threshold*100;%percentage
%     threshold(threshold<epsilon)=0;
    index_s=find(~threshold);
    
    start_ind = 1;
    for i = 1:length(partition_tmp)
        sel = start_ind:cum_part(i);   %norm(w(sel1),2)/(fx(i)*M)
        if ismember(i, index_s)
        W(sel)= 0;
        end
        start_ind = cum_part(i) + 1;
    end
    
   
end
    