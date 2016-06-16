function shape = LoadShape( filename )

number_of_points = 98;

data = textread(filename, '%s' );

data = data( 21:((98*3)+21-1) );

shape= [];
for i=1:length(data)
    shape1(i) = str2num( data{i} );
end;

shape = zeros(98*2,1);
shape(1:2:end) = shape1(1:3:end);
shape(2:2:end) = shape1(2:3:end);