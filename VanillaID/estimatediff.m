function [data_new, derivative, index] = estimatediff(data, t, type, kth, order)

% Estimate the derivative
% version 1.0 : we only estimate the first derivative, i.e. order = 1.
%            we will add higher order estimation in the future version
% please cite Kris De Brabanter's JMLR paper: Derivative Estimation with Local Polynomial Fitting



%%
% clear all; close all
% 
% type = 'solver';
% %                 type = 'central';
% %                 type = 'eularforward';
% %                 type = 'eularbackward';
% tt= 0:pi/10:24;
% y1 = tt.*sin(tt/2)-tt;
% y1 = y1(:);
% y2 = cos(tt);
% y2 = y2(:);
% data = [y1, y2];
% noise =0*randn(size(data));
% data = data+noise;
% derivative_true1 = sin(tt/2) + tt.*cos(tt/2)/2-1;
% derivative_true2 = -sin(tt);
% derivative_true = [derivative_true1(:), derivative_true2(:)];
% t = tt;
% which_var  = 1;
% kth =2; % please specify the k value


%% ================== solver =======================
plot_onoff='off';

n = size(data,1);

if ~exist('t') || isempty(t)
    t = 1:n;
    t = t(:);
end

if ~isempty(t)
    t = t(:);
    if length(t)~=n
        error('dimension of t is not consistent');
    end
end

if ~exist('type') || isempty(type)
    type = 'central';
end

if strcmp(type, 'central')
    type = 'solver';
    kth = 1;
end

if ~exist('order') || isempty(order)
    order = 1;
end


switch type
    case 'solver'
        %%
        index = kth+1:n-kth;
        
        for j  = 1:1:kth
            w(j) = (6*j^2)/(kth*(kth+1)*(2*kth+1));
        end
        
        t_new = t(index,:);
        
        % for i = 1: size(data_new,1)
        for i = index
            for j = 1:kth
                diff_i(j,:) = w(j) * ...
                    (data(i+j,:) - data(i-j,:))./(t(i+j,:) - t(i-j,:));
            end
            derivative(i,:) = sum(diff_i,1);
            d_forward(i,:) = (data(i+1,:) - data(i,:))./(t(i+1,:) - t(i,:));
            d_backward(i,:) = data(i,:) - data(i-1,:)./(t(i,:) - t(i-1,:));
        end
        derivative(1:kth,:) = [];
        d_forward(1:kth,:) = [];
        d_backward(1:kth,:) = [];
        
    case 'eularforward'
        %%
        index  = 1: n-1;
        for i = index
            d_forward(i,:) = (data(i+1,:) - data(i,:))./(t(i+1,:) - t(i,:));
        end
        derivative = d_forward(index,:);
        
    case 'eularbackward'
        %%
        index  = 2: n;
        for i = index
            d_backward(i,:) = data(i,:) - data(i-1,:)./(t(i,:) - t(i-1,:));
        end
        derivative = d_backward(index,:);
        index  = index - 1;
end

index = index(:);
data_new = data(index,:);

if size(data_new,1) ~= size(derivative,1)
    error('please check dimensions')
end


%% to comment
if strcmp(plot_onoff,'on')
    derivative_true = derivative_true(index,:);
    % derivative_true = normalize_var(derivative_true, -1, 1);
    % derivative = normalize_var(derivative, -1, 1);
    % d_forward = normalize_var(d_forward, -1, 1);
    
    
    figure; plot(data_new(:,which_var),'r');
    figure; plot(derivative(:,which_var),'b-*'); hold on; plot(derivative_true(:,which_var),'r-o');
    hold on; plot(d_forward(:,which_var),'g')
    legend('estimated derivative','true derivative','eular')
end
