function y = inner(a,b);
% This is a MatLab function to compute the inner product of
% two vectors a and b.  
% By Ralph Howard on 1/12/98
% Call syntax: y = inner(a,b) or inner(a,b)
% Input: The two vectors a and b 
% Output: The value of the inner product of a and b.
c=0;			% intialize the variable c
n= length(a);		% get the lenght of the vector a
for k=1:n		% start the loop
	c=c+a(k)*b(k);	% update c by the k-th product in inner product
end			% end loop
y = c;			% print value of c = inner product
