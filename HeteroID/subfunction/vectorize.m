function p=vectorize(x)
% clear all
% x{1} =[1,2,3]';
% x{2} = [4,5];
    [m,n]=size(x);
        
    if  strcmp(class(x),'double') 
        p=reshape(x, m*n,1);
    elseif strcmp(class(x),'cell')
        for i =1:n
            ni(i) = max(size(x{i}));
        end
        cum_part = cumsum(ni);
        start_ind = 1;
        for i = 1:n,
            p(start_ind:cum_part(i),1) =  vec(x{i});
            start_ind = cum_part(i)+1;
        end
    end
    
 