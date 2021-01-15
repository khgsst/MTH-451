%power iteration Algorithm 27.1
clear, close all
N=50;
itmax = 400;
TOL = 100*eps;
vnew = zeros(N,itmax);
ivenw=vnew;
lamb=zeros(itmax,1);
ilamb=lamb;
x=linspace(0,1,N);
D=diag(exp(-x)/N);
evmin=exp(-1); evmax = 1;

U = randn(N);
A = U*D*U';

v=randn(N,1);
v=v/norm(v);
ivnew(:,1)=v;
ilamb(1,1) = ivnew(:,1)'*A*ivnew(:,1);
icnt = 1;
ierr(icnt,1) = norm(A*ivnew(:,1)-ilamb(1,1)*ivnew(:,1),inf);
for k = 2:itmax
    w = A*ivnew(:,k-1);
    ivnew(:,k)= w/norm(w);
    ilamb(k,1) = ivnew(:,k)'*A*ivnew(:,k);
    ierr(icnt,1) = norm(A*ivnew(:,k)-ilamb(k,1)*ivnew(:,k),inf);
    if ierr(icnt,1)< TOL, break, end;
    icnt = icnt+1;
end
subplot(2,1,1), plot(1:icnt,ilamb(1:icnt,1))
ylabel('eigenvalue')
xlabel('iterate')
subplot(2,1,2), plot(ivnew(:,icnt))
ylabel('eigenvector')
xlabel('N')
pause
close
%v=randn(N,1);
%v=v/norm(v);
vnew(:,1)=v;
lamb(1,1) = vnew(:,1)'*A*vnew(:,1);
cnt = 1;
err(1,1) = norm(A*vnew(:,1)-lamb(1,1)*vnew(:,1),inf);
for k = 2:itmax    
    B=(A-lamb(k-1,1)*eye(N));
    w=B\vnew(:,k-1);    
    vnew(:,k)= w/norm(w);
    lamb(k,1) = vnew(:,k)'*A*vnew(:,k);
     err(cnt,1) = norm(A*vnew(:,k)-lamb(k,1)*vnew(:,k),inf);
    if err(cnt,1)< TOL, break, end;
    cnt = cnt+1;
end
subplot(2,1,1), plot(1:cnt,lamb(1:cnt,1))
ylabel('eigenvalue')
xlabel('iterate')
subplot(2,1,2), plot(vnew(:,cnt))
ylabel('eigenvector')
xlabel('N')
pause
close
loglog(1:icnt,ierr(1:icnt),'-r*')
hold
loglog(1:cnt,err(1:cnt),'-ko')
grid
xlabel('iterate')
ylabel('error, L_\infty')
legend('Rayleigh','Iterated Rayleigh')
fprintf('power method estimate of eigenvalue=%12.8f \n',ilamb(icnt))
fprintf('iterated method estimate of eigenvalue=%12.8f \n',lamb(cnt))
disp('Matlab estimated Eigenvalues:')
evs=eig(A)

