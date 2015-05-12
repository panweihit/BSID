function [y, A,w_true] = GRN_dis(which,noise)
ga1=0.3;ga2=0.4;ga3=0.5;ga4=0.2;ga5=0.4;ga6=0.6;
be1=1.4;be2=1.5;be3=1.6;
al1=4;al2=3;al3=5;
h=4;
% w=rand(6+h*12,6);
% w_tru=[diag([-ga1,-ga2,-ga3,-ga4,-ga5,-ga6])+[zeros(3,3) diag([be1,be2,be3]);zeros(3,6)];[zeros(45,3);0 al2 0;0 0 al3;al1 0 0] zeros(48,3)];


w_tru=[diag([-ga1,-ga2,-ga3,-ga4,-ga5,-ga6])+...
    [zeros(3,3) diag([be1,be2,be3]);zeros(3,6)];...
    [zeros(15,3);0 al2 0;0 0 al3;al1 0 0;zeros(30,3);] zeros(48,3)];

x=zeros(6,1);
T=101;M=T-1;  % number of time points
X=zeros(T,6); %records of states: T rows
X(1,:)=rand(1,6); %set initial state
Phi=zeros(M,6+h*12);
Y=zeros(M,6);
samplerate = 1;
for k=2:T %real time from 0:1:50
    %call hill function:H
    Phi(k-1,:)=[X(k-1,:) hill(X(k-1,:),1) hill(X(k-1,:),2) hill(X(k-1,:),3) hill(X(k-1,:),4)];
    X(k,:)=X(k-1,:) + samplerate*Phi(k-1,:)*w_tru + noise*randn(1,6);
    Y(k-1,:)=  (X(k,:)-X(k-1,:))/samplerate; %M rows, derivative of x
end

% which state variable you are interested in
%which = 1; % 1 ~ 6
y = Y(:,which);
A = Phi;
w_true = w_tru(:,which);
% check construction accuracy
err = norm(y-A*w_true)


figure; 

plot(X); 
