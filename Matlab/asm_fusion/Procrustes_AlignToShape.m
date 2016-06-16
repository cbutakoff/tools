function [alignedShape, T] = Procrustes_AlignToShape (inputShape, referenceShape)

% [alignedShape, T] = Procrustes_AlignToShape (inputShape, referenceShape)
%
% Aligns the inputShape to the referenceShape by means of
% Procrustes Analysis
%                                     a  -b  tx
% T is a 3z3 matrix defined as        b   a  ty
%                                     0   0   1
%
% such that we can use it as a projective transformation
% to obtain alignedShape_projective = T * inputShape_projective
% Notice, however, that the way we define shapes does not stand
% for projective representation. Shape2proj and proj2shape functions
% does the appropiate conversions. 
%

% Center both shapes to the origin
% -------------------------------------------------
nL = floor (size (inputShape, 2) / 2);
cnt_input = ShapeCentroid (inputShape, (1 : nL));
cnt_ref = ShapeCentroid (referenceShape, (1 : nL));

iShape = inputShape - xxyy_to_xyxy (cnt_input(1) * ones(1, nL), cnt_input(2) * ones(1, nL));
rShape = referenceShape - xxyy_to_xyxy (cnt_ref(1) * ones(1, nL), cnt_ref(2) * ones(1, nL));

% Center both shapes to the origin
% -------------------------------------------------  
modulo2 = sum (iShape .* iShape);
aa = sum (iShape .* rShape) / modulo2;

ix = iShape (1 : 2 : end);
iy = iShape (2 : 2 : end);
rx = rShape (1 : 2 : end);
ry = rShape (2 : 2 : end);

bb = sum (ix .* ry - iy .* rx) / modulo2;
out_x = aa * ix - bb * iy;
out_y = bb * ix + aa * iy;
alignedShape = xxyy_to_xyxy (out_x, out_y);

% Compute transformation parameters
% -------------------------------------------------  
tx = -aa * cnt_input(1) + bb * cnt_input(2) + cnt_ref(1);
ty = -bb * cnt_input(1) - aa * cnt_input(2) + cnt_ref(2);

T = [aa, -bb, tx; bb, aa, ty; 0, 0, 1];