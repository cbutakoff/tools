function M = xxyy_to_xyxy ( x, y );

% M = xxyy_to_xyxy ( x, y );;
% 
% Converts from 'x' and 'y' to 'M' where
% M is a matrix with EACH ROW containing a 2D 
% shape points x1, y1, x2, y2.... 
% The function will separate x's and y's
% always ordering the shapes on ecah row
% of the matrices
%

M = [x x];
M(:, 1 : 2 : end) = x;
M(:, 2 : 2 : end) = y;
