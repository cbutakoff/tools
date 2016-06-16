function PlotShape (shapeSet, groups, varargin);

% PlotShape (shapeSet, groups, charOpt);
% Plots a landmark-based 2D shape
% 'shape' is in format x1, y1, x2, y2, ...., xN, yN
% 'groups' is a 4-colomn vector, defining 1st landmark, last landmark
% (of the group numbered from 1 to N), open or closed shape (0 or 1)
% and normal direction (1 or -1) of the contour
% 'charOpt' is an option to pass to the plot function

nLandmarks = floor (size (shapeSet, 2) / 2);
nGroups = size (groups, 1);
hold on

for nS = 1 : size (shapeSet, 1)
    shape = shapeSet (nS, :);
    x = shape (1 : 2 : 2 * nLandmarks - 1);
    y = shape (2 : 2 : 2 * nLandmarks);
    
    for g = 1 : nGroups
        xp = x(groups(g, 1) : groups(g, 2));
        yp = y(groups(g, 1) : groups(g, 2));
        if groups (g, 3) == 1
            xp = [xp, xp(1)];
            yp = [yp, yp(1)];
        end
        plot ( xp, yp, varargin{:});
    end
end
hold off
