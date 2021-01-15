clear,close

%  Comparison of the interpolation polynomial  Pn
%  for sin(x) in the interval [1 2 pi] with sin(x)
%  Pn(x) = c_1 x^n + ..+ c_n x + c_{n+1}
%  The goal is to show that this interpolation
%  strategy is highly ill-conditioned due to the
%  Vandermonde matrix generated when Nx is large.


  for ii = 1:4
	     Nx = 2^(2*ii);
x = (linspace(0,2*pi,Nx))';
     dx = x(2)-x(1);
% Find the interpolation polynomial using
% the Vandermunde matrix
     V = vander(x);
     C = V\sin(x); 
     rr = cond(V);
     fprintf('Nx=%d, condition number=%e \n',Nx,rr);
% Plot the results at grid values other than the 
% ones used in the interpolation.
     y = x+0.1*dx*rand(size(x));     
     figure(1)
     subplot(2,2,ii), plot(y,polyval(C,y)-sin(y),'linewidth',2)
     xlabel('x','fontsize',18)
     ylabel('Pn(x)-f(x)','fontsize',18)
     title(['Nx = ',num2str(Nx)])     
     figure(2)
     plot(y,polyval(C,y),y,sin(y),'linewidth',2)
     legend('interpolated','sin(x)')
     xlabel('x','fontsize',18)
     dis('hit any key to continue')
     pause
  end

