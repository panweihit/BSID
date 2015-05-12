function [y, A,w_true] = GRN_con(which)
ga1=0.3;ga2=0.4;ga3=0.5;ga4=0.2;ga5=0.4;ga6=0.6;
be1=1.4;be2=1.5;be3=1.6;
al1=4;al2=3;al3=5;
h=4;

% w_tru=[diag([-ga1,-ga2,-ga3,-ga4,-ga5,-ga6])+[zeros(3,3) diag([be1,be2,be3]);zeros(3,6)];...
%     [zeros(39,3);0 al2 0;0 0 al3;al1 0 0;zeros(6,3)] zeros(48,3)];
w_tru=[diag([-ga1,-ga2,-ga3,-ga4,-ga5,-ga6])+...
    [zeros(3,3) diag([be1,be2,be3]);zeros(3,6)];...
    [zeros(15,3);0 al2 0;0 0 al3;al1 0 0;zeros(30,3);] zeros(48,3)];
[t,x]=ode45(@(t,x) odefun(t,x,w_tru),[0:1:50],abs(rand(6,1)));

for k=1:50
    if k==1
        h=[ hill(x(k,:),1) hill(x(k,:),2) hill(x(k,:),3) hill(x(k,:),4)];
    else
        h=[h;hill(x(k,:),1) hill(x(k,:),2) hill(x(k,:),3) hill(x(k,:),4)];
    end
end


k = 1;
[~, derivative]=estimatediff(x, t, 'solver', k, []);
y = derivative(:,which);
Phi =[ x(1:(end-1),:) h];
Phi = Phi(k+1:end-k+1,:);
w_true = w_tru(:,which);

A = Phi(5:end,:);
y = y(5:end,:);

norm(y-A*w_true)
figure;
plot(y,'r');hold on;
plot(A*w_true)

figure;
plot(t,x)
title('Solutions ');
xlabel('t');
ylabel('x');
grid on
% subplot(2,1,2)
% plot(x(:,1),x(:,2))
% title('Phase-plane portrait');
% xlabel('x1');
% ylabel('x2');
% grid on
