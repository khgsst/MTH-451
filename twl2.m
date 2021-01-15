clear,close all
for n=2:30
    m=2*n-1;
    x=linspace(-1,1,n)';
    y=linspace(-1,1,m)';
    V=fliplr(vander(x));
    W=fliplr(vander(y));
    W=W(:,1:n);
    A=W/V;
    nrm(n)=norm(A,inf);
end
N=[2:30];
est=(2.^N)./(exp(1)*log(N).*(N-1));
semilogy(N,nrm(2:30),'+',N,est,'-');
grid
set(gca,'fontsize',16);
xlabel('n')
ylabel('||A||_\infty')
title('Infinity norm vs dimension problem 12.2')
legend('Computed','Estimated')