function out = fcls(HIM,M)
% Fully Constrained Linear Spectral Unmixing
% Perform a Linear least squares with nonnegativity constraints.


[ns,nb] = size(HIM);
[l,p] = size(M);

Delta = 1/1000; 

N = zeros(l+1,p);
N(1:l,1:p) = Delta*M;
N(l+1,:) = ones(1,p);
s = zeros(l+1,1);

OutputImage = zeros(ns,p);

go=0;
for i = 1:ns
    s(1:l) = Delta*HIM(i,:)';
    s(l+1) = 1;
    if go==0
        [Abundances,options] = lsqnonneg_fast(N,s); % solve NNLS model (E.q. 19)
        go=1;
    else
        Abundances = lsqnonneg_fast(N,s,options);
    end
    OutputImage(i,:) = Abundances;
end


out = OutputImage;