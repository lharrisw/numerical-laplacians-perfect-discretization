function R = tetragon(x,y,varargin)
%TETRAGON Convex rectangular region inside a 2-D mesh.
%   R = TETRAGON(X,Y) returns the default region of a square positioned at
%   the center of the grid with horizontal side length of L/2 and vertical
%   side length of H/2, where L is the length of the x-interval and H is
%   the length of the y-interval. 
%
%   R = TETRAGON(X,Y,...) returns the default square, the default rhombus,
%   or a general region. TETRAGON accepts these parameters: 0, 1, 2,
%   and p. The parameter p is a vector with eight elements used to specify
%   the vertices of the desired quadrilateral. The x-values and y-values
%   of the vertices are arranged in the following manner:
%   a1 < a2,a3 < a4, b1 < b2,b3 < b4. 
%
%   Note: Built-in MATLAB functions like FLIP, ROT90 can be used to achieve
%   upside down and 90 degree rotated regions.
%   
%       Option value:       Output:
%               0               default square
%               1               default square
%               2               default rhombus
%               p               quadrilateral with vertices 
%                                   p = [a1 a2 a3 a4 b1 b2 b3 b4]
%
%   Example: Plot the function f(x,y) = 1 inside a rectangular region given
%            by the vertices (-1.5,0), (-pi/4,-2), (1.1,1.75), and (2,-1),
%            where -2 <= x,y <= 2.
%       
%       [X,Y] = meshgrid(linspace(-2,2));
%       P = [-1.5 1.1 -pi/4 2 -2 0 -1 1.75];
%       R = tetragon(X,Y,P);
%       f = ones(size(X)).*R;
%       mesh(X,Y,f)
%
%   See also ELLIPSE, TRIANGLE.

% Copyright 2020 Latimer De'Shone Harris-Ward.

% default parameters
dx = (x(1,end) - x(1,1))/4; % step-size for x-interval
dy = (y(end,1) - y(1,1))/4; % step-size for y-interval
a1 = x(1,1) + dx;           % x-value for bottom left vertex
a3 = x(1,1) + 3*dx;         % x-value for bottom right vertex
b1 = y(1,1) + dy;           % y-value for bottom vertices
b3 = y(1,1) + 3*dy;         % y-value for top vertices
s1 = y >= b1;               % this and the next four lines define the
s2 = y <= b3;               % rectangular region
s3 = x >= a1;
s4 = x <= a3;

if nargin == 2
    % TETRAGON(X,Y)
    % Generate the default square
    R = s1 & s2 & s3 & s4;
elseif nargin == 3
    % TETRAGON(X,Y,options)
    if numel(varargin{1}) == 1
        if varargin{1} == 0 || varargin{1} == 1
            % Generate the default square, the same as above
            R = s1 & s2 & s3 & s4;
        elseif varargin{1} == 2
            % Generate default rhombus/diamond
            a2 = x(1,1) + 2*dx;
            a3 = a2;
            a4 = x(1,1) + 3*dx;
            b2 = y(1,1) + 2*dy;
            b3 = b2;
            b4 = y(1,1) + 3*dy;
            m1 = (b4 - b2)/(a2 - a1);
            m2 = (b3 - b4)/(a4 - a2);
            m3 = (b3 - b1)/(a4 - a3);
            m4 = (b1 - b2)/(a3 - a1);
            s1 = y <= m1*(x - a1) + b2;
            s2 = y <= m2*(x - a4) + b3;
            s3 = y >= m3*(x - a3) + b1;
            s4 = y >= m4*(x - a3) + b1;
            R = s1 & s2 & s3 & s4;
        else
            errornum = ['The only values tetragon accepts for the third'...
                        ' argument are 0, 1, and 2. Otherwise, pass a'...
                       ' vector with 8 elements to specify the vertices.'];
            error(errornum)
        end
    else
        % Generate general quadrilateral
        p = varargin{1};
        % p must be a vector with 8 elements
        if numel(p) ~= 8
            errorp = ['There must be eight values to specify the four'...
                      ' vertices of the rectangle.'];
            error(errorp)
        end
        % returns matrix of zeros if p contains all zeros
        if nnz(p) == 0
            msg = ['The vector consists of all zeros and the region is'...
                   ' zero.'];
            warning(msg)
        end
        % make sure vertices lie within the region
        for k = 1:numel(p)
            if p(k) > x(1,end) || p(k) > y(end,1) || p(k) < x(1,1) ||...
                    p(k) < y(1,1)
                error('The vertices must lie within the grid.')
            end
        end
        a1 = p(1);
        a2 = p(2);
        a3 = p(3);
        a4 = p(4);
        b1 = p(5);
        b2 = p(6);
        b3 = p(7);
        b4 = p(8);
        % preserve the order for x and y values of vertices
        if a1 > a4
            errora1 = ['The first value in the vector must be less'...
                       ' than the fourth value.'];
            error(errora1)
        end
        if b1 > b4
            errorb1 = ['The fifth value in the vector must be less'...
                       ' than the eighth value.'];
            error(errorb1)
        end
        if a2 < a1 || a3 < a1 || a2 > a4 || a3 > a4
            errora = ['The second and third values of the vector must'...
                      ' lie between the first and the fourth values.'];
            error(errora)
        end
        if b2 < b1 || b3 < b1 || b2 > b4 || b3 > b4
            errorb = ['The sixth and seventh values of the vector must'...
                      ' lie between the fifth and the eigth values.'];
            error(errorb)
        end
        m1 = (b4 - b2)/(a2 - a1);
        m2 = (b3 - b4)/(a4 - a2);
        m3 = (b3 - b1)/(a4 - a3);
        m4 = (b1 - b2)/(a3 - a1);
        s1 = y <= m1*(x - a1) + b2;
        s2 = y <= m2*(x - a4) + b3;
        s3 = y >= m3*(x - a3) + b1;
        s4 = y >= m4*(x - a3) + b1;
        R = s1 & s2 & s3 & s4;
    end
end
end
