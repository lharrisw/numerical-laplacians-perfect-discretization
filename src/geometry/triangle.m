function R = triangle(x,y,p)
%TRIANGLE Triangular region inside a 2-D mesh.
%   R = TRIANGLE(X,Y) is a matrix consisting of ones inside the default
%   triangular region and zeros outside the triangular region. The default
%   region is an equilateral triangle at the center of the grid.
%                                                                               
%   R = TRIANGLE(X,Y,P), where P is a vector with 6 elements, returns the
%   triangular region defined by the points P = [a1 a2 a3 b1 b2 b3], where
%   the a-values correspond to x-values and the b-values correspond to
%   y-values. When paired, the points represent the vertices of the
%   triangle. The x-values and y-values are arranged in the following
%   manner: a1 <= a2 <= a3 and b1,b2 <= b3
%   
%   Note: Built-in MATLAB functions like FLIP, ROT90 can be used to achieve
%   upside down and 90 degree rotated triangular regions.
%
%   Example: Plot the following two-variable function inside a triangular
%            region with vertices (-.5,-.5), (.25,1.75), and (1,.5):
%            f(x,y) = x^2+y^2 for -1 <= x <= 1 and -2 <= y <= 2.
%
%       [X,Y] = meshgrid(-1:.1:1,-2:.2:2);
%       P = [-.5 .25 1 -.5 1.75 .5];
%       R = triangle(X,Y,P);
%       f = (X.^2 + Y.^2) .* R;
%       mesh(X,Y,f)
%
%   See also ELLIPSE, TETRAGON.

%   Copyright 2020 Latimer De'Shone Harris-Ward.

% Return the default region of an equilateral triangle at center of mesh
if nargin == 2 || isempty(p) 
    % TRIANGLE(X,Y)
    dx = (x(1,end) - x(1,1))/4;     % step-size for x-interval
    dy = (y(end,1) - y(1,1))/4;     % step-size for y-interval
    a1 = x(1,1) + dx;               % x-value of left vertex
    a2 = x(1,1) + 2*dx;             % x-value of top vertex
    a3 = x(1,1) + 3*dx;             % x-value of right vertex
    b1 = y(1,1) + dy;               % y-value of bottom vertex
    b2 = b1;                        % y-value of bottom vertex
    b3 = y(1,1) + 3*dy;             % y-value of top vertex
    m1 = (b3 - b1)/(a2 - a1);       % slope of left leg
    m2 = (b2 - b3)/(a3 - a2);       % slope of right leg
    m3 = (b2 - b1)/(a3 - a1);       % slope of bottom leg
    leg1 = y <= m1*(x - a1) + b1;   % left leg
    leg2 = y <= m2*(x - a2) + b3;   % right leg
    leg3 = y >= m3*(x - a1) + b1;   % bottom leg
    R = leg1 & leg2 & leg3;         % region
else
    % TRIANGLE(X,Y,P)
    if length(p) ~= 6
        errorp = ['You need a vector of 6 values to specify the three'...
                  ' vertices of the triangle.'];
        error(errorp)
    end
    a1 = p(1);
    a2 = p(2); 
    a3 = p(3);
    b1 = p(4);
    b2 = p(5);
    b3 = p(6);
    % making sure vertices are within the bounds
    for k = 1:numel(p)
        if p(k) > x(1,end) || p(k) > y(end,1) || p(k) < x(1,1) ||...
                p(k) < y(1,1)
            error('The vertices must lie within the grid.')
        end
    end
    % making sure the order of the x-values is preserved
    if a3 < a2
        error1 = ['The x-value of the right-most vertex must be greater'...
                  ' than or equal to the x-value of the top-most vertex.'];
        error(error1)
    end
    if a1 > a2
        error2 = ['The x-value of the left-most vertex must be less' ...
                  ' than or equal to the x-value of the top-most vertex.'];
        error(error2)
    end
    % making sure order of y-values are preserved
    if (b1 > b2 && b3 < b1 && b3 > b2) || (b1 > b2 && b3 < b2) || ...
            (b2 > b1 && b3 < b2 && b3 > b1) || (b2 > b1 && b3 < b1)
        error3 = ['The y-value of the top-most vertex must be greater'...
                  ' than both the y-value of the left and right-most' ...
                  ' vertices.'];
        error(error3)
    end
    m1 = (b3 - b1)/(a2 - a1); 
    m2 = (b2 - b3)/(a3 - a2); 
    m3 = (b2 - b1)/(a3 - a1); 
    leg1 = y <= m1*(x - a1) + b1; 
    leg2 = y <= m2*(x - a2) + b3; 
    leg3 = y >= m3*(x - a1) + b1; 
    R = leg1 & leg2 & leg3; 
end
end
