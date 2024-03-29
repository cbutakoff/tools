function SaveShape( filename, shape )

start = {...
'<?xml version="1.0"?>',...
'<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">',...
'  <PolyData>',...
'    <Piece NumberOfPoints="98" NumberOfVerts="0" NumberOfLines="2" NumberOfStrips="0" NumberOfPolys="6">',...
'      <Points NumberOfDimensions="2">',...
'        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'};

ending = {...
'   </DataArray>',...
'      </Points>',...
'      <PointData>',...
'      </PointData>',...
'      <CellData>',...
'      </CellData>',...
'      <Verts>',...
'        <DataArray type="Int32" Name="connectivity" format="ascii">',...
'        </DataArray>',...
'        <DataArray type="Int32" Name="offsets" format="ascii">',...
'        </DataArray>',...
'      </Verts>',...
'      <Lines>',...
'        <DataArray type="Int32" Name="connectivity" format="ascii">',...
'          0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 85 86 87 88 89 90 91 92 93 94 95 96 97 ',...
'        </DataArray>',...
'        <DataArray type="Int32" Name="offsets" format="ascii">',...
'          25 38 ',...
'        </DataArray>',...
'      </Lines>',...
'      <Strips>',...
'        <DataArray type="Int32" Name="connectivity" format="ascii">',...
'        </DataArray>',...
'        <DataArray type="Int32" Name="offsets" format="ascii">',...
'        </DataArray>',...
'      </Strips>',...
'      <Polys>',...
'        <DataArray type="Int32" Name="connectivity" format="ascii">',...
'          25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 ',...
'        </DataArray>',...
'        <DataArray type="Int32" Name="offsets" format="ascii">',...
'          12 24 34 44 52 60 ',...
'        </DataArray>',...
'      </Polys>',...
'    </Piece>',...
'  </PolyData>',...
'</VTKFile>'};

x = shape(1:2:end);
y = shape(2:2:end);

data = {};
for i=1:(length(shape)/2)
    data{i} = [ '          ',num2str(x(i)),' ', num2str(y(i)), ' 0' ];
end;

save_cells( [start, data, ending], filename );
        