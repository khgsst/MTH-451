% This is a demo of Conjugate Gradient that
% visualizes the iterative process.

clear, close all
ITMAX=10; TOL=1e-10; DECEL=0.98; N=2;
fprintf('DECEL should be set to 1 for regular CG\n')
fprintf('DECEL set to %f\n',DECEL);

xhist=zeros(2,ITMAX);
p0=zeros(2,1);
eigvals=[exp(-randn),exp(-randn)];
dd=diag(eigvals);
U=randn(2);
[Q,R]=qr(U);
A=Q*dd*Q';
[vv,DD]=eig(A);

xs=randn(2,1);
b = A*xs;
phis=0.5*xs'*A*xs - xs'*b;

vectarrow(p0,DD(1,1)*vv(:,1))
hold
vectarrow(p0,DD(2,2)*vv(:,2))
axis('equal')
pause
close

xx=linspace(-3,3,100);
yy=xx;
[XX,YY]=meshgrid(xx,yy');
for jj = 1:100
    for ii = 1:100
        X=[xx(1,ii);yy(1,jj)];
        phi(ii,jj) = 0.5*X'*A*X- X'*b;
    end
end
colormap('cool')
contourf(YY,XX,phi)
axis('equal')
hold
x0 = xs(1);
y0 = xs(2);
x1 = DD(1,1)*vv(1,1)+x0;
y1 = DD(1,1)*vv(2,1)+y0;
plot([x0;x1],[y0;y1],'linewidth',2);
x1 = DD(2,2)*vv(1,2)+x0;
y1 = DD(2,2)*vv(2,2)+y0;
plot([x0;x1],[y0;y1],'linewidth',2);          
p0=[3;-3];
[conv,xp,xhist]=cg(A,b,p0,ITMAX,TOL,DECEL,N);   
for ii = 2:ITMAX      
    x0=xhist(1,ii-1);x1=xhist(1,ii);
    y0=xhist(2,ii-1);y1=xhist(2,ii);
    plot([x0;x1],[y0;y1],'linewidth',2);          
    pause
end

figure(2)
semilogy((1:ITMAX),conv(1:ITMAX),'o',...
    (1:ITMAX),conv(1:ITMAX),'--');
xlabel('iteration, n','fontsize',18)
ylabel('Inf Norm(r(n))/Norm(b)','fontsize',18)
grid
