function res = tform_vector( vec, T )
% T - 2x2 matrix
% vec - column vector divisible by 2

res = vec;
for i=1:2:length(vec)
    res(i:i+1) = T*vec(i:i+1);
end;