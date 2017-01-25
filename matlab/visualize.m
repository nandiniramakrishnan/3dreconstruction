% Q2.7 - Todo:
% Integrating everything together.
% Loads necessary files from ../data/ and visualizes 3D reconstruction
% using scatter3

load('../data/templeCoords.mat');       % x1 and y1
load('../data/some_corresp.mat');       % pts1 and pts2
load('../data/intrinsics.mat');         % K1 and K2
im1 = imread('../data/im1.png');
im2 = imread('../data/im2.png');

x2 = zeros(length(x1), 1);
y2 = zeros(length(y1), 1);
M = max(size(im1, 1), size(im1, 2));    % Scaling factor
[ F ] = eightpoint( pts1, pts2, M );
for i = 1:length(x1)
[ x2(i), y2(i) ] = epipolarCorrespondence( im1, im2, F, x1(i), y1(i) );
end

% Triangulate to compute the 3d locations
M1 = [eye(3), zeros(3, 1)];
p1 = [x1, y1];
load('q2_5.mat','M2');                       % M2
p2 = [x2, y2];
[ P, error ] = triangulate( K1*M1, p1, K2*M2, p2 );
scatter3(P(:,1),P(:,2),P(:,3));
save('q2_7.mat','F','M1','M2');