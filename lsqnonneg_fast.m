function [x,options] = lsqnonneg_fast(C,d,options)

% Reference:
%  Lawson and Hanson, "Solving Least Squares Problems", Prentice-Hall, 1974.

if nargin<3
    %disp('Calculating optimset');
    
    defaultopt = optimset('display','final','TolX','10*eps*norm(C,1)*length(C)');
    if ~isreal(C) || ~isreal(d), error('C and d must be real.'); end
    options = [];
    options = optimset(defaultopt,options);
    printtype = optimget(options,'display');
    tol = optimget(options,'tolx');
    
    % In case the defaults were gathered from calling: optimset('fminsearch'):
    c=C;
    if ischar(tol)
        tol = eval(tol);
    end
    
    switch printtype
        case {'none','off'}
            verbosity = 0;
        case 'iter'
            warning('''iter'' value not valid for ''Display'' parameter for lsqnonneg');
            verbosity = 2;
        case 'final'
            verbosity = 1;
        otherwise
            error('Bad value for options parameter: ''Display''');
    end
end
tol = optimget(options,'tolx');
c=C;
if ischar(tol)
    tol = eval(tol);
end

[m,n] = size(C);
P = zeros(1,n);
Z = 1:n;
x = P';

ZZ=Z;
resid = d-C*x;
w = C'*(resid);

% set up iteration criterion
outeriter = 0;
iter = 0;
itmax = 3*n;
exitflag = 1;

% outer loop to put variables into set to hold positive coefficients
while any(Z) && any(w(ZZ) > tol)
   outeriter = outeriter + 1;
   [wt,t] = max(w(ZZ));
   t = ZZ(t);
   P(1,t) = t;
   Z(t) = 0;
   PP = find(P);
   ZZ = find(Z);
   nzz = size(ZZ);
   CP(1:m,PP) = C(:,PP);
   CP(:,ZZ) = zeros(m,nzz(2));
   z = pinv(CP)*d;
   z(ZZ) = zeros(nzz(2),nzz(1));
   % inner loop to remove elements from the positive set which no longer belong
   while any((z(PP) <= tol))
      iter = iter + 1;
      if iter > itmax
         if verbosity 
            warnstr = sprintf('Exiting: Iteration count is exceeded, exiting LSQNONNEG.', ...
               '\n','Try raising the tolerance (OPTIONS.TolX).');
            disp(warnstr);
         end
         exitflag = 0;
         output.iterations = outeriter;
         resnorm = sum(resid.*resid);
         x = z;
         lambda = w;
         return
      end
      QQ = find((z <= tol) & P');
      alpha = min(x(QQ)./(x(QQ) - z(QQ)));
      x = x + alpha*(z - x);
      ij = find(abs(x) < tol & P' ~= 0);
      Z(ij)=ij';
      P(ij)=zeros(1,length(ij));
      PP = find(P);
      ZZ = find(Z);
      nzz = size(ZZ);
      CP(1:m,PP) = C(:,PP);
      CP(:,ZZ) = zeros(m,nzz(2));
      z = pinv(CP)*d;
      z(ZZ) = zeros(nzz(2),nzz(1));
   end
   x = z;
   resid = d-C*x;
   w = C'*(resid);
end

lambda = w;
resnorm = sum(resid.*resid);
output.iterations = outeriter;
output.algorithm = 'active-set using svd';

verbosity =0;

if verbosity > 0
   disp('Optimization terminated successfully.');   
end

