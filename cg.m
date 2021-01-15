% This program performs the conjugate gradient method,
% 
function [conv, x,xhist] = cg(A,b,x,ITMAX,TOL,DECEL,N)
xhist=zeros(N,ITMAX);
r = b - A*x;
p = r;
j = 0;
conv(1) = norm(r);
if b == 0
   err(1) = norm(x,inf);
end

while (norm(r)/norm(b) > TOL) && (j < ITMAX)
    j = j + 1;
    xhist(:,j)=x;
    alpha = DECEL*(r'*r)/(p'*A*p);
    x = x + alpha*p;    
    rnew = r - alpha*A*p;
    beta = (rnew'*rnew)/(r'*r);
    p = rnew + beta*p;
    r = rnew;
    conv(j) = norm(r)/norm(b);
    if b == 0
       err(j) = norm(x,inf);
    end
end

return
