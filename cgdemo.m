% A simple demonstration of the conjugate gradient algorithm
% and the method of steepest descents for quadratic minimization

A = gallery('wathen',5,9); 
% a random sparse SPD finite element mass matrix
n=length(A);
b=ones(n,1);

% initial guess
x=zeros(n,1);
% number of iterations
niter = 2*n;
% store norm of residual in res
res=zeros(niter+1,1);

% conjugate gradients
r = b-A*x;
s = r;
r2 = r'*r;
res(1) = sqrt(r2);
for i = 1:niter
  t = A*s;
  s2 = s'*t;
  lambda = r2/s2;
  x = x + lambda*s;
  r = r-lambda*t;
  r2old = r2;
  r2 = r'*r;
  s=r+(r2/r2old)*s;
  res(i+1)=sqrt(r2);
end

% plot norms of residuals: linear convergence means a straight line
semilogy(res)

hold on

% for comparison show steepest descents

% initial guess
x=zeros(n,1);
% number of iterations
niter = 2*n;
% store norm of residual in res
ressd=zeros(niter+1,1);

% steepest descents
for i = 1:niter
  s = b-A*x;
  s2 = s'*s;
  l = s2/(s'*A*s);
  x = x + l*s;
  ressd(i)=sqrt(s2);
end
s = b-A*x;
ressd(niter+1)=sqrt(s'*s);

% plot norms of residuals: linear convergence means a straight line
semilogy(ressd,'r')
grid on

xlabel('iterations')
ylabel('norm of residual')
legend('conjugate gradient','steepest descents')
disp('hit any key to close window')
pause
close
