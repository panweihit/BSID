function [y,A, x_true, partition] = multidic(Dic, Y, W)
% Dic and Y are cells!

% clear all
% Dic{1} = randn(2,3);
% Dic{2} = magic(3);
% Y{1}=[1;2];
% Y{2} =[1 2 3]';
% W{1}  = [1 2 3]';
% W{2} = [4 5 6]';

%% 
blocknumber = length(Dic);
for i = 1:1:blocknumber
    n(i) = size(Dic{i},2);
end
N = mean(n);
n(n==N)=0;

if norm(n)~=0
    error('N of each block is not the same')
end

if length(Y)~=blocknumber
    error('block not compatible');
end

for i = 1:1:blocknumber
    if length(Y{i})~=size(Dic{i}, 1)
        error('dimension of Y{i} and Dic{i} not same');
    end
    blocksize(i) = length(Y{i});
end

blockcumsum =  cumsum(blocksize);

for i = 1:1:blocknumber
    A_block(blockcumsum(i)-blocksize(i)+1:blockcumsum(i),(i-1)*N+1: i*N) = Dic{i};
    W_block((i-1)*N+1: i*N,1) = W{i}; 
    y(blockcumsum(i)-blocksize(i)+1:blockcumsum(i),1) = Y{i};
end


sortcumsum = cumsum(kron(N, ones(blocknumber,1)));

A = [];
x_true = [];
for j = N-1:-1:0
    A_tmp = A_block(:, sortcumsum-j);
    W_tmp = W_block(sortcumsum-j,1);
    x_true = [x_true; W_tmp];
    A = [A, A_tmp];
end

partition = blocknumber*ones(N,1);
% x_true = full(x_true);
