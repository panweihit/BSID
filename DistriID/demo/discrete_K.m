%%  "Kuramoto"
function [y, A, w,snr]=discrete_K(N, delta, tspan,sigma, initial,candidate_num,scale,select)

% clear all
% select=2;
% N=100; %Define number of nodes
% delta = 0.1;  %Define sample interval
% tspan = 25;  %Define time span
% sigma = 0e-4;  %Define intensity of noise
% initial=pi*(1+randn(N,1));    %  initial value: here the initial value is randomly chosen
% candidate_num=1;
% scale=1e0;


interval_ini=-15;
interval_end=15;
sparsity=0.1;
%uniform distributed
weight0= sprand(N,N,sparsity);
weight=interval_ini*spones(weight0)+(interval_end-interval_ini)*weight0;


        sparsity2 = round(N*sparsity);
       
        x_true = zeros(N,1);
        x_true_nnz = interval_ini+(interval_end-interval_ini)*rand(sparsity2,1);
        itr = randperm(N);
        x_true(itr(1:sparsity2)) = x_true_nnz;
        weight(select,:) = x_true';



for j=1:1:N
    weight(j,j)=0;
end


omega=1*rand(N,1);


T=tspan/delta; %Sample number
x=zeros(N,T); % Initialization of state variable

x(:,1)=initial;%  initial value

%sigma=0.1;
NOISE=sigma*randn(N,T-1);


for t=1:1:T-1 
    for i=1:1:N
        x(i,t+1)=x(i,t)+delta*(omega(i)+weight(i,:)*sin(x(:,t)-ones(N,1)*x(i,t)))+NOISE(i,t);
    end
end
% % plot((0:T-1)*delta,x(1:5,:)); 
% plot((0:T-1)*delta,x(1,:)); 
% title('Phase Dynamics of Kuramoto Model'); 
% xlabel('time t');
% ylabel('Phase \phi');
% % legend('\phi_1', '\phi_2','\phi_3', '\phi_4','\phi_5');
for i=1:1:N
    xx(i,:)=x(i,:)-floor(x(i,:)./(2*pi))*2*pi;
end
plot((0:T-1)*delta,xx(1,:),'-*'); 
axis([0 tspan 0 2*pi]);
title('Phase Dynamics of Kuramoto Model'); 
xlabel('time t');
ylabel('Phase \phi');


%% 
% scale=1e1;
% y_matrix=zeros(T-1,N);
% w_matrix=zeros(candidate_num*N+1,N);
% phi_matrix=zeros(T-1,candidate_num*N*(N+1));
%for j=1:1:n
%j=27;
for j=1:1:N 
%     y=scale^2*(x(j,2:T)-x(j,1:T-1))'/delta;
    y=scale^2*(x(j,2:T)-x(j,1:T-1))';
    A=scale*delta*[sin(x(:,1:T-1)-ones(N,1)*x(j,1:T-1));...
        x(:,1:T-1)-ones(N,1)*x(j,1:T-1); ...
        (sin(x(:,1:T-1)-ones(N,1)*x(j,1:T-1))).^2;...
        (cos(x(:,1:T-1)-ones(N,1)*x(j,1:T-1))).^2;...
        cos(x(:,1:T-1)-ones(N,1)*x(j,1:T-1));...
%         (sin(x(:,1:T-1)-ones(n,1)*x(j,1:T-1))).* (cos(x(:,1:T-1)-ones(n,1)*x(j,1:T-1)));...
        ones(1,T-1)]';
        %]';
%     for q=1:1:candidate_num
%         phi(:,j+(q-1)*N)=randn(T-1,1);
%     end
A = A(:,[1: N*candidate_num,end]);
    w=scale*[weight(j,:), zeros(1, (candidate_num-1)*N), omega(j)]';
    y_matrix{j}=y;
    A_matrix{j}=A;
    w_matrix{j}=w;
end
e=sigma;


y=y_matrix{j};
A=A_matrix{j};
w=w_matrix{j};   

indzero =find(~norms(A, 2,1));
A(:,indzero) = [];
w(indzero) =[];



noise=scale^2*NOISE(select,:)';
err = norm(y-A*w-noise)
snr = 20*log10(norm(A*w)/norm(noise))
