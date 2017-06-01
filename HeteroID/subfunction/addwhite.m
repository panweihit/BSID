function y = addwhite(y, perc)
% perc=input('% of additive noise=');

dev=std(y)*perc/100;
noise=randn(size(y)).*kron(dev,ones(size(y,1),1));
y=y+noise;

% save corrupted y perc;
 