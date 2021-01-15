function svdexample(percent)
% To run, type svdexample(X), where
% X is 0 <= X <= 100, the percentage of
% singular values that you want to use in 
% the compressed image. On output, a png file
% will be created and displayed. It will be the 
% compressed image.

% Note: to find out the actual compression achieved in the
% png image, query the size of the original file png and the 
% one produced by this function.

% This function requires you obtain the png file
% lena.png to run.

 I = imread('lena.png');
 figure(1)
 imshow(I);
 
 fprintf('FIGURE 1 shows uncompressed image \n')
 [m n c] = size(I);
 I = double(I);
 R = I(:, :, 1);
 G = I(:, :, 2);
 B = I(:, :, 3);
 [RU, RS, RV] = svd(R);
 % TIPS: svd(R) returns diag(RS) because it is the most important
 [GU, GS, GV] = svd(G);
 [BU, BS, BV] = svd(B);
 RdS = diag(RS);
 GdS = diag(GS);
 BdS = diag(BS);
 D = length(RdS);
 DIM = floor(D*(percent*0.01));
 RdS(DIM+1:end) = 0; % reduce dimensions
 GdS(DIM+1:end) = 0;
 BdS(DIM+1:end) = 0;
 RSR = sparse(1:D, 1:D, RdS, m, n); % reconstruct into matrix
 GSR = sparse(1:D, 1:D, GdS, m, n);
 BSR = sparse(1:D, 1:D, BdS, m, n);
 RR = RU * RSR * RV'; % recover
 GR = GU * GSR * GV';
 BR = BU * BSR * BV';
 IR = cat(3, RR, GR, BR);
 IR = uint8(IR);
 figure(2)
 imshow(IR);
 fprintf('FIGURE 2 shows compressed image. \n')
 imwrite(IR, sprintf('lenasvd%02d.png', percent), 'PNG');
 disp('hit any key to close off displayed image')
 pause
close all

end