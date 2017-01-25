function [ P, error ] = triangulate( M1, p1, M2, p2 )
% triangulate:
%       M1 - 3x4 Camera Matrix 1
%       p1 - Nx2 set of points
%       M2 - 3x4 Camera Matrix 2
%       p2 - Nx2 set of points

% Q2.4 - Todo:
%       Implement a triangulation algorithm to compute the 3d locations
%       See Szeliski Chapter 7 for ideas
%

nCorr = size(p1, 1);        % Number of correspondences
P = zeros(nCorr, 3);        % 3D Point locations
error = 0;                  % Reprojection error

% Iterate through all point correspondences
for i = 1:nCorr
    A = [p1(i,1)*M1(3,:) - M1(1,:); ...
    p1(i,2)*M1(3,:) - M1(2,:); ...
    p2(i,1)*M2(3,:) - M2(1,:);  ...
    p2(i,2)*M2(3,:) - M2(2,:)];
    [~, ~, V] = svd(A);
    P(i,:) = V(1:3,end)/V(4, end);
    
    % Calculate reprojection error
    p1_proj = (M1*[P(i,:), 1]');
    p2_proj = (M2*[P(i,:), 1]');
    p1_proj = p1_proj/p1_proj(3);       % z = 0?
    p2_proj = p2_proj/p2_proj(3);
    p1i = [p1(i,:), 1];
    p2i = [p2(i,:), 1];
    error = error + norm(p1i - p1_proj')^2 + norm(p2i - p2_proj')^2;
end
end