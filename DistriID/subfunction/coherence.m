function miu = coherence(A)
% clear all;A = RandOrthMat(5);
[m,n] = size(A);

MIU = zeros(n,n);
for i = 1:n
    Ai = A(:,i);
    for j = (i+1):n
        Aj = A(:,j);
        MIU(i,j) = abs(Ai'*Aj) / (norm(Ai)*norm(Aj));
    end
end

miu = max(max(MIU));
end
