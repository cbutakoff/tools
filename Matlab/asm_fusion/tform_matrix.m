function res = tform_matrix(matr, T)
%tranforms columns of a matrix

[r,c] = size( matr );

res = zeros(r,c);
for i=1:c
    res(:,i) = tform_vector( matr(:,i), T );
end;
