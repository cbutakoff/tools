function [fused, new_means] = fuse_asms( eigenspaces )
% input: 
% eignespaces - cell array of structures
% .V - eigenvectors (in columns)
% .D - eigenvalues
% .N - number of observations
% .m - mean
% .w - weight
% Returns the same structure

% check that weights sum to 1
sum_of_probab = 0;
sum_of_weights = 0;
fused.N = 0;
fused.m = zeros( size( eigenspaces{1}.m ) );

for i=1:length(eigenspaces)
    sum_of_weights = sum_of_weights + eigenspaces{i}.w;
    sum_of_probab = sum_of_probab + eigenspaces{i}.w * eigenspaces{i}.N;
    fused.N = fused.N + eigenspaces{i}.N;
end;

%fused.m = fused.m / sum_of_probab;

if abs( sum_of_weights - 1.0 )>eps
    'incorrect weights. Must sum up to 1'
    return
end;

% form an array of modified weights (weight per sample, not per model as given)
weights = [];
weights( length(eigenspaces) ) = 0;
for i=1:length(eigenspaces)
    weights(i) = eigenspaces{i}.w / sum_of_probab;
end;

% first estimate of mean
% for i=1:length(eigenspaces)
%     fused.m = fused.m + eigenspaces{i}.m * (eigenspaces{i}.N * weights(i));
% end;

%fused.m = fused.m / norm(fused.m,2);

fused.m = eigenspaces{i}.m;

%make a copy of eigenspaces
eigenspaces1 = eigenspaces;

% align all the means to the estimate
for i=1:length(eigenspaces)
    centroid = [ mean(eigenspaces1{i}.m(1:2:end)); mean(eigenspaces1{i}.m(2:2:end)) ];
    eigenspaces1{i}.m = eigenspaces1{i}.m - ...
        repmat( centroid,  length(eigenspaces1{i}.m)/2, 1 );
end;
            
procrustes_finished = false;
while ~procrustes_finished
    new_mean = zeros( size( eigenspaces1{1}.m ) );
    for i=1:length(eigenspaces1)
        [eigenspaces1{i}.m, T] = Procrustes_AlignToShape(eigenspaces1{i}.m', fused.m');
        eigenspaces1{i}.m = eigenspaces1{i}.m';
        new_mean = new_mean + eigenspaces1{i}.m * (eigenspaces1{i}.N * weights(i));
    end;
    new_mean = new_mean/norm(new_mean,2);
     
    if sum( (new_mean-fused.m).^2 ) < 1e-20
        procrustes_finished = true;
    end;
    
   fused.m = new_mean;
end;


clear eigenspaces1;

%get the transforms and transform all the eigenspaces
tforms={};
new_means = {};
for i=1:length(eigenspaces)
    [new_means{i}, T] = Procrustes_AlignToShape(eigenspaces{i}.m', fused.m');
    new_means{i} = new_means{i}';
    tforms{i} = T(1:2,1:2);
    eigenspaces{i}.V = tform_matrix( eigenspaces{i}.V, tforms{i} );
end;


%put all eigenvectors together
H = [];
for i=1:length(eigenspaces)
    if eigenspaces{i}.w > eps
        H = [H, eigenspaces{i}.V];
    end;
end;

%add differences between means
for i=1:(length(eigenspaces)-1)
    for j=i+1:length(eigenspaces)
        if ( eigenspaces{i}.w ~= 0 ) & ( eigenspaces{j}.w ~= 0 )
            %diff = eigenspaces{i}.w * eigenspaces{i}.m - eigenspaces{j}.w * eigenspaces{j}.m;
            diff = new_means{i} - new_means{j};
            H = [H, diff];
        end;
    end;
end;

%H_trans = Orthonormalize( H )';
H_trans = orth( H )';

% Create matrix A
[rows, cols]= size(H_trans);
A = zeros( rows, rows );

for i=1:length(eigenspaces)
    matrix = H_trans * eigenspaces{i}.V;
    A = A + weights(i) * (eigenspaces{i}.N-1) * matrix * eigenspaces{i}.D * matrix';
end;

%create temporary matrix set 2
for i=1:(length(eigenspaces)-1)
    for j=i+1:length(eigenspaces)
        diff = eigenspaces{i}.m - eigenspaces{j}.m;
        matrix = H_trans * diff;
        matrix = ( matrix * matrix' )*weights(i)*weights(j)*eigenspaces{i}.N*eigenspaces{j}.N ;
        A = A + matrix;
    end;
end;

%get unbiased estimate
s = 0;
for i=1:length(eigenspaces)
    s = s + weights(i)*weights(i)*eigenspaces{i}.N;
end;
A = A./(1-s);

%do the rest of things
[R fused.D] = eig(A);

fused.V = H_trans' * R; 

[sorted_val eval_indexes] = sort( diag(fused.D), 'descend' );
fused.D = diag(sorted_val);
fused.V = fused.V(:,eval_indexes);

