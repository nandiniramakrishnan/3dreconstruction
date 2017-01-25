function [ F ] = sevenpoint( pts1, pts2, M )
% sevenpoint:
%   pts1 - Nx2 matrix of (x,y) coordinates
%   pts2 - Nx2 matrix of (x,y) coordinates
%   M    - max (imwidth, imheight)

% Q2.2 - Todo:
%     Implement the eightpoint algorithm
%     Generate a matrix F from some '../data/some_corresp.mat'
%     Save recovered F (either 1 or 3 in cell), M, pts1, pts2 to q2_2.mat

%     Write recovered F and display the output of displayEpipolarF in your writeup

% Load correspondence data
load('../data/some_corresp.mat');

% Scale the coordinates for normalization
T = [2/M, 0, -1; 0, 2/M, -1; 0, 0, 1];
seven_points = 52:58;
seven_pts_1 = [];
seven_pts_1 = cat(2, seven_pts_1, pts1(seven_points,1));
seven_pts_1 = cat(2, seven_pts_1, pts1(seven_points,2));
seven_pts_1 = cat(2, seven_pts_1, ones(7,1));
seven_pts_2 = [];
seven_pts_2 = cat(2, seven_pts_2, pts2(seven_points,1));
seven_pts_2 = cat(2, seven_pts_2, pts2(seven_points,2));
seven_pts_2 = cat(2, seven_pts_2, ones(7,1));

pts1_norm = (T*seven_pts_1')';
pts2_norm = (T*seven_pts_2')';

x_l = pts1_norm(:,1);
y_l = pts1_norm(:,2);
x_r = pts2_norm(:,1);
y_r = pts2_norm(:,2);

X = [x_r.*x_l x_r.*y_l x_r y_r.*x_l y_r.*y_l y_r x_l y_l];
X = cat(2, X, ones(7,1));
[~,~,V] = svd(X);

F1 = V(:, end);
F1 = reshape(F1, 3, 3)';

% enforce rank-2 constraint
[UU,DD,VV] = svd(F1);
DD(3,3) = 0;
F1 = UU*DD*VV';

% unscale
F1 = T'*F1*T;

F2 = V(:, end - 1);
F2 = reshape(F2, 3, 3)';
% enforce rank-2 constraint
[UU,DD,VV] = svd(F2);
DD(3,3) = 0;
F2 = UU*DD*VV';
% unscale
F2 = T'*F2*T;

syms lam_var;
F_formula = (1-lam_var)*F1 + lam_var*F2;
lam = vpa(solve(det(F_formula) == 0, lam_var));
for i = 1:length(lam)
    F{i} = subs(F_formula, lam(i));
end

save('q2_2.mat','F','M','pts1','pts2');
end
