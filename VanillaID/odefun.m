function [f ] = odefun(t,x,w, noise)
f=([x' hill(x,1) hill(x,2) hill(x,3) hill(x,4)]*w)'+noise*randn(6,1);
% x(1) x(2) x(3) x(4) x(5) x(6)
end
