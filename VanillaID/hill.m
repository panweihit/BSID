function [ hillfunction ] = hill( x,h )
hillfunction=zeros(1,12);
for i=1:12
    if i<=6
        hillfunction(i)=1/(1+x(i)^h);
    else
        hillfunction(i)=x(i-6)^h/(1+x(i-6)^h);
    end
end
end

