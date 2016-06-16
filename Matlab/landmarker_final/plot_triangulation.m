function plot_triangulation( varargin )

x = varargin{1}; 
y = varargin{2}; 
tri = varargin{3};
plot_data = varargin{4:end};

tx = ones(4,1);
ty = ones(4,1);
hold on; 
for t = 1:size(tri,1)
    tx(1:3) = x( tri(t,:) );
    tx(4) = tx(1);
    ty(1:3) = y( tri(t,:) );
    ty(4) = ty(1);
    plot(tx,ty,plot_data);
end;

