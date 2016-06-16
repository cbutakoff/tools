function [x, y] = xyxy_to_xxyy ( M );

% [x, y] = xyxy_to_xxyy ( M );
% 
% Converts from 'M' to 'x' and 'y' where
% M is a matrix with EACH ROW containing a 2D 
% shape points x1, y1, x2, y2.... 
% The function will separate x's and y's
% always ordering the shapes on ecah row
% of the matrices
%

x = M (:, 1 : 2 : end);
y = M (:, 2 : 2 : end);
