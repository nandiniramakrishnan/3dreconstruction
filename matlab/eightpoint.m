function [ F ] = eightpoint( pts1, pts2, M )

% eightpoint:
%   pts1 - Nx2 matrix of (x,y) coordinates
%   pts2 - Nx2 matrix of (x,y) coordinates
%   M    - max (imwidth, imheight)

% Q2.1 - Todo:
%     Implement the eightpoint algorithm
%     Generate a matrix F from some '../data/some_corresp.mat'
%     Save F, M, pts1, pts2 to q2_1.mat

%     Write F and display the output of displayEpipolarF in your writeup

nCorr = size(pts1, 1);

% Scale the coordinates for normalization
T = [2/M, 0, -1; 0, 2/M, -1; 0, 0, 1];
pts1_homogeneous = cat(2, pts1, ones(nCorr, 1))';
pts2_homogeneous = cat(2, pts2, ones(nCorr, 1))';
pts1_norm = (T*pts1_homogeneous)';
pts2_norm = (T*pts2_homogeneous)';

x_l = pts1_norm(:,1);
y_l = pts1_norm(:,2);
x_r = pts2_norm(:,1);
y_r = pts2_norm(:,2);

X = [x_r.*x_l x_r.*y_l x_r y_r.*x_l y_r.*y_l y_r x_l y_l];
X = cat(2, X, ones(nCorr,1));
[~,~,V] = svd(X);
last_eiegnvector = V(:, end);
F = reshape(last_eiegnvector, 3, 3)';

% enforce rank-2 constraint
[UU,DD,VV] = svd(F);
DD(3,3) = 0;
F = UU*DD*VV';

% unscale
F = T'*F*T;

save('q2_1.mat','F','M','pts1','pts2');
end
