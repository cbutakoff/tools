function c = ShapeCentroid (shape, landmarks)

% c = ShapeCentroid (shape, landmarks)
%

idx = landmarks * 2 - 1;
c(1) = mean (shape (idx));
c(2) = mean (shape (idx + 1));



