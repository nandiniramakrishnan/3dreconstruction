% Q2.5 - Todo:
%       1. Load point correspondences
%       2. Obtain the correct M2
%       3. Save the correct M2, p1, p2, R and P to q2_5.mat

load('../data/some_corresp.mat');   % pts1 and pts2 = point correspondences
load('../data/intrinsics.mat');     % K1 and K2 intrinsic matrices
im1 = imread('../data/im1.png');    % Read image to get the scale factor
M = max(size(im1,1), size(im1,2));  % M = Scale factor

% Get fundamental matrix using 8-point algorithm
[ F ] = eightpoint( pts1, pts2, M );
% Get essential matrix
[ E ] = essentialMatrix( F, K1, K2 );
[M2s] = camera2(E);

M1 = [eye(3), zeros(3, 1)];

% Iterate through all 4 camera matrices
for i = 1:size(M2s, 3)
    % Obtain 3d point locations
    [ Pi, error ] = triangulate( K1*M1, pts1, K2*M2s(:, :, i), pts2 );
    
    % Check if 3d point locations all have positive z-values
    if isempty(find(Pi(:,3) < 0, 1))
        M2 = M2s(:, :, i);
        P = Pi;
        break;
    end
end
p1 = pts1;
p2 = pts2;
save('q2_5.mat','M2','p1','p2','P');