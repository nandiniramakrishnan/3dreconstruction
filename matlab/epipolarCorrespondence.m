function [ x2, y2 ] = epipolarCorrespondence( im1, im2, F, x1, y1 )
% epipolarCorrespondence:
%       im1 - Image 1
%       im2 - Image 2
%       F - Fundamental Matrix between im1 and im2
%       x1 - x coord in image 1
%       y1 - y coord in image 1

% Q2.6 - Todo:
%           Implement a method to compute (x2,y2) given (x1,y1)
%           Use F to only scan along the epipolar line
%           Experiment with different window sizes or weighting schemes
%           Save F, pts1, and pts2 used to generate view to q2_6.mat
%
%           Explain your methods and optimization in your writeup
x2 = 0;
y2 = 0;
imheight = size(im1, 1);
window_size = 16;
min_diff = 1000000;
T1 = y1;
B1 = y1 + window_size;
L1 = x1;
R1 = x1 + window_size;
win1 = im1(T1:B1, L1:R1);
X = [x1;y1;1];
L = F*X;
a = L(1);
b = L(2);
c = L(3);
for y2i = 1:(imheight - window_size)
    x2i = (-c-b*y2i)/a;
    T2 = y2i;
    B2 = y2i + window_size;
    L2 = x2i;
    R2 = x2i + window_size;
    win2 = im2(T2:B2, L2:R2);
    diff = norm(double(win2) - double(win1));
    if diff < min_diff
        min_diff = diff;
        x2 = x2i;
        y2 = y2i;
    end
end
end
