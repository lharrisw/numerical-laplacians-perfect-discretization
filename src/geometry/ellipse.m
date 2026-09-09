function R = ellipse(x,y,center,axes,angle)
%ELLIPSE Elliptical region inside a 2-D mesh.
%   R = ELLIPSE(X,Y) is a matrix of ones inside the default region and
%   zeros lying outside the region. The default region is a circle with
%   AXES = [L/4 H/4], centered at CENTER = [L/2 H/2], and ANGLE = 0, where
%   L is the length of the x-axis and H is the length of the y-axis.
%   
%   R = ELLIPSE(X,Y,CENTER,AXES,ANGLE) describes the ellipse specified by
%   the values of CENTER, AXES, and ANGLE. CENTER is a two-element vector
%   of the form CENTER = [H K] where H specifies the x-value and K
%   specifies the y-value of the center of the ellipse. AXES is a
%   two-element vector AXES = [MAJOR MINOR] is MAJOR specifies the major
%   axis of the ellipse and MINOR is the minor axis of the ellipse. ANGLE
%   is the angle between the major axis and the positive x-axis in radians.
%   
%   Note: Ideally, inputs for AXES should be nonnegative real numbers. For
%   negative inputs, the aboslute value is used by definition.           
%
%   Example: Plot the following two-variable function inside an elliptical
%            region with major axis of .5, minor axis of .25, a center of
%            (h,k) = (1/3,0), and rotated +45 degrees.
%            f(x,y) = x^2+y^2 for -2 <= x <= 2 and -2 <= y <= 2
%
%       [X,Y] = meshgrid(-2:.1:2,-2:.2:2);
%       CENTER = [1/3 0];
%       AXES = [.5 .25];
%       ANGLE = pi/4;
%       R = ellipse(X,Y,CENTER,AXES,ANGLE);
%       f = (X.^2 + Y.^2) .* R;
%       mesh(X,Y,f)
%
%   See also TETRAGON, TRIANGLE.

%   Copyright 2020  Latimer De'Shone Harris-Ward

if nargin < 3 || (isempty(center) && isempty(axes) && isempty(angle))
    dx = (x(1,end) - x(1,1))/4; % step-size of .25*(length of x-interval) 
    dy = (y(end,1) - y(1,1))/4; % step-size of .25*(length of y-interval)
    h = x(1,1) + 2*dx;          % x-value of center
    k = y(1,1) + 2*dy;          % y-value of center
    a = dx;                     % major axis
    b = dy;                     % minor axis
    A = 0;                      % angle
    R = (((x - h)*cos(A) + (y - k)*sin(A)).^2/a^2 + ...
         ((x - h)*sin(A) - (y - k)*cos(A)).^2/b^2) <= 1;
else
    if length(center) ~= 2
        error('Two values are needed to specify the center of the ellipse.')
    end
    if length(axes) ~= 2
        error('Two values are needed to specify the major and minor axes.')
    end
    h = center(1);
    k = center(2);
    a = axes(1);
    b = axes(2);
    % ensure no division by zero
    if a == 0 || b == 0
        warning('Divide by zero.')
        error('The major and minor axes must be nonzero.')
    end
    % making sure the center, minor, and major axes in region
    if h < x(1,1) || h > x(1,end) || k < y(1,1) || k > y(end,1)
        error('The center must lie within the region.')
    end
    A = angle; 
    R = (((x - h)*cos(A) + (y - k)*sin(A)).^2/a^2 + ...
         ((x - h)*sin(A) - (y - k)*cos(A)).^2/b^2) <= 1;
end
end
