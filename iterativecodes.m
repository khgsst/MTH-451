clear; close all;
format long;
N=10; % size of matrix
% Set iteration tolerance and max num of iterations here:
%TOL = input('enter tolerance:  ');
TOL = 1e-12;
Nmax = 5000;

% load matrix problem 
%       Ax = B
%
U = randn(N);
A=U*U'+150*eye(N); % matrix is diagonally dominant for small N
rr=cond(A);
fprintf('Condition number of matrix = %f \n \n',rr);
x=randn(N,1); % solution
B = A*x;      % RHS
disp('PROGRAM WILL FIRST SOLVE USING JACOBI');
disp('THEN IT WILL SOLVE USING GAUSS-SEIDEL OR SOR');
counts=0;

xp=zeros(N,1);xo=xp;
%           SOR AND GAUSS SEIDEL:
% initial guess assumed to be the null vector:
clear error;

%
%	Successive Over-Relaxation  Technique:
%       Gauss-Seidel: same as SOR but w = 1
%       set w=1 in SOR for Gauss-Seidel:
xp = zeros(N,1);
counts = 0;
w = input('Enter 1 for G-S. Enter SOR factor 1 <= w < 2:   ')
% uncomment to test:
%w = 1.2;
tic;
for k=1:Nmax
    counts = counts + 1;
    for i=1:N
        s =0; 
        for j=1:N
            if i ~= j, s = s + A(i,j)*xp(j,1); end
        end
        xp(i,1) = w*(-s + B(i,1))/A(i,i) + (1-w)*xp(i,1);
    end
    error(k) = norm(A*xp-B,inf);
%      fprintf(' It. Num. = %3.0f, error = %15.13e\n', ...
%                k, error(k))
    if error(k) < TOL, break; end
    if k == Nmax, 
         fprintf('\n MAXIMUM NUMBER OF ITERATIONS EXCEEDED\n'), 
    end;
end
tgs = toc;
      fprintf(' \n \n time=%f, It. Num. = %3.0f, error = %15.13e\n\n\n', ...
                tgs, k, error(counts))
norm(x-xp)
semilogy(1:counts,error(1:counts));
grid on
xlabel('iteration number')
ylabel('error')
if w==1, title('GAUSS-SEIDEL ITERATION TECHNIQUE'), end;
if w~=1
   omeg=num2str(w);
   title(['SOR ITERATION WITH OMEGA =', omeg])
end;
disp('PRESS ANY KEY TO CONTINUE')
pause

clear error;

%
%	Successive Over-Relaxation  Technique:
%       Gauss-Seidel: same as SOR but w = 1
%       set w=1 in SOR for Gauss-Seidel:
xp = zeros(N,1);
counts = 0;
w = input('Enter 1 for G-S. Enter SOR factor 1 <= w < 2:   ')
% uncomment to test:
%w = 1.2;
tic;
for k=1:Nmax
    counts = counts + 1;
    for i=1:N
        s =0; 
        for j=1:N
            if i ~= j, s = s + A(i,j)*xp(j,1); end
        end
        xp(i,1) = w*(-s + B(i,1))/A(i,i) + (1-w)*xp(i,1);
    end
    error(k) = norm(A*xp-B,inf);
%      fprintf(' It. Num. = %3.0f, error = %15.13e\n', ...
%                k, error(k))
    if error(k) < TOL, break; end
    if k == Nmax, 
         fprintf('\n MAXIMUM NUMBER OF ITERATIONS EXCEEDED\n'), 
    end;
end
tgs = toc;
      fprintf(' \n \n time=%f, It. Num. = %3.0f, error = %15.13e\n\n\n', ...
                tgs, k, error(counts))
norm(x-xp)
semilogy(1:counts,error(1:counts));
grid on
xlabel('iteration number')
ylabel('error')
if w==1, title('GAUSS-SEIDEL ITERATION TECHNIQUE'), end;
if w~=1
   omeg=num2str(w);
   title(['SOR ITERATION WITH OMEGA =', omeg])
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%5
