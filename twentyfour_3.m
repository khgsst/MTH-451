clear, close all
nn=100; %number of cases
tt=linspace(0,20,nn); % time
expnorms = tt; %initialization

for jj = 1:10
   A = randn(10)-2*eye(10); % random matrices with a shift of -2
                            % so some negative values get generated.
                            
  [ww,dvs]=eig(A);
  sigma = max(real(diag(dvs))); % y axis
  for k = 1:nn
      expnorms(k) = norm(expm(tt(k)*A));
  end
  % plot the norms and compare to the exp(t sigma) and kappa exp(tsigma)
  wcond = cond(ww);
  semilogy(tt,exp(tt*sigma),'-r', ...
      tt,wcond*exp(tt*sigma),'-k','linewidth',1);
   hold on
  semilogy(tt,expnorms,'b','linewidth',2); 
end
hold off
xlabel('t')
ylabel('norm(e^{At})')
title('10-dim random matrices, norm of expm(At)')
grid
text(2,10^10,'exp(t*sigma) in red')
text(2,10^12,'exp(tA) norm in blue')
text(2,10^14,'w*exp(t*sigma) in black')

  