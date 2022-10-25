/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include <vtkSmartPointer.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkShortArray.h>
#include <vtkCell.h>
#include <vtkPointLocator.h>
#include <vtkType.h>
#include <vtkPointData.h>
#include <vtkCellDataToPointData.h>
#include <vtkCellLocator.h>
#include <vtkPolyDataReader.h>
#include <vtkCellData.h>
#include <vtkCell.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkXMLPolyDataReader.h>

#include <vtkCallbackCommand.h>
#include <VTKCommonTools.h>
#include <vtkType.h>


#include <Eigen/Dense> 

#include <algorithm>
#include <fstream>
#include <vector>
#include <typeinfo>
#include <set>


using namespace std;

typedef struct __FaceInfo {
    vtkIdType boundary_id;
    vector<vtkIdType> pt_id;
    vtkIdType vol_cell_id; //id of the volumetric mesh cell touching the face
} FaceInfo;


template <typename T>
using MatrixRX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using MatrixRXd = MatrixRX<double>; //these are the types for MPIO, do not change
using MatrixRXi = MatrixRX<int64_t>;

template <class T>
void WriteMPIOMatrix( MatrixRX<T>& m, const char* problem_name, const char *variable, const char* association, const char* dtype );


template<class T>
bool Contains(const std::vector<T> &list, T x)
{
	return std::find(list.begin(), list.end(), x) != list.end();
}
    

inline std::vector<vtkIdType> intersection(std::vector<vtkIdType> &v1,
                                    std::vector<vtkIdType> &v2){
    std::vector<vtkIdType> v3;

    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    std::set_intersection(v1.begin(),v1.end(),
                          v2.begin(),v2.end(),
                          back_inserter(v3));
    return v3;
}


template <typename T>
std::ostream& operator<<(std::ostream& output, std::vector<T> const& values)
{
    for (auto const& value : values)
    {
        output << value << std::endl;
    }
    return output;
}


bool FileExists( const char* filename, bool no_exception=false )
{
	bool result = true;

	FILE *fid = fopen(filename,"r");
	if( fid == NULL )
	{
		std::cerr<<"File not found: "<<filename;

		char error[500];
		sprintf(error,"File not found: %s",filename);

		if( !no_exception ) throw std::runtime_error( error );

		result = false;
	}
	else fclose(fid);

	return result;
}


int main(int argc, char** argv)
{
    std::cout<<"Parameters: volmesh.vtk|vtu surfmesh.vtk|vtp arrayname volmesh_out_prefix -f {VTK,ALYATXT,ALYAMPIO} -s scale -o correct_orientation(0|1)"<<std::endl;
    std::cout<<"Scale - rescale mesh by this factor, !!! for now works only for ALYA, for everything else put 1"<<std::endl;
    std::cout<<"correct_orientation - 1 or 0 whether to correct cell orientation or not (normally 0, default 0)"<<std::endl;
    std::cout<<"Labels - boundary ids stored as celldata in surfmesh.vtk. Ids MUST BE > 0"<<std::endl;
    std::cout<<"-f can be repeated to save more formats. Default ALYAMPIO"<<std::endl;         

    enum enumFileFormats  { VTK, ALYATXT, ALYAMPIO, ELMER };

    {
        double test=0;
        if( sizeof(test) != 8 ) 
        {
            cout<<"Your compiler's double type is not 8 bytes. Fix the sources"<<endl;
            throw;
        }
    }


    if(argc<4) return -1;
    
    int c =1;
    const char* volmesh_file = argv[c++];
    const char* surfmesh_file = argv[c++];
    const char* array_name = argv[c++];
    const char* outfile_prefix = argv[c++];
    vector<enumFileFormats> save_formats;
    float scale = 1;
    bool correct_orientation = false;


    while( c<argc )
    {

        if( strcmp( argv[c], "-f" ) == 0 ){
            c++;
            if( strcmp(argv[c], "VTK") == 0 ){
                save_formats.push_back(VTK);
            }
            else if( strcmp( argv[c], "ALYATXT" ) == 0 ){
                save_formats.push_back(ALYATXT);
            }
            else if( strcmp( argv[c], "ALYAMPIO" ) == 0 ){
                save_formats.push_back(ALYAMPIO);
            }
            else if( strcmp( argv[c], "ELMER" ) == 0 ){
                save_formats.push_back(ELMER);
            }
        }
        else if( strcmp( argv[c], "-s" ) == 0 ){
            scale = atof(argv[++c]);
        }
        else if( strcmp( argv[c], "-o" ) == 0 ){
            correct_orientation = atoi(argv[++c])==1;
        }

        c++;
    }


    
    
    std::cout<<"Volumetric mesh: "<<volmesh_file<<std::endl;
    std::cout<<"Surface mesh: "<<surfmesh_file<<std::endl;
    std::cout<<"Label array: "<<array_name<<std::endl;
    std::cout<<"Output mesh prefix: "<<outfile_prefix<<std::endl;
    std::cout<<"Scale: "<<scale<<std::endl;
    if(correct_orientation)
        std::cout<<"Correcting orientation: ON"<<std::endl;
    else
        std::cout<<"Correcting orientation: OFF"<<std::endl;
    
    FileExists(volmesh_file);
    FileExists(surfmesh_file);

    if (scale<0) {
        cout<<"Scale factor must be > 0"<<endl;
        throw;
    }

    if (save_formats.size()==0) {
        save_formats.push_back(ALYAMPIO);
    }


    /*===================================================================================
    *            
    *            read and process data
    *
    *
    * ===================================================================================*/

    vtkSmartPointer<vtkUnstructuredGrid> volmesh =     vtkSmartPointer<vtkUnstructuredGrid>::New();
    if (volmesh_file[ strlen(volmesh_file)-1 ] == 'k') //vtk
    {

        vtkSmartPointer<vtkDataSetReader> vol_rdr = vtkSmartPointer<vtkDataSetReader>::New();
        CommonTools::AssociateProgressFunction(vol_rdr);

        vol_rdr->SetFileName(volmesh_file);
        vol_rdr->Update();
        volmesh->DeepCopy( (vtkUnstructuredGrid*)vol_rdr->GetOutput() );
    }    
    else
    {
        vtkSmartPointer<vtkXMLUnstructuredGridReader> vol_rdr = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        CommonTools::AssociateProgressFunction(vol_rdr);

        vol_rdr->SetFileName(volmesh_file);
        vol_rdr->Update();
        volmesh->DeepCopy( (vtkUnstructuredGrid*)vol_rdr->GetOutput() );
    }

    //check just in case tha the mesh is not surface
    //really a silly check for one face being a triangle
    if(volmesh->GetCell(0)->GetNumberOfPoints()<4)
    {
        std::cout<<"Supplied volumetric mesh has cells with : "<<volmesh->GetCell(0)->GetNumberOfPoints()<<" vertices"<<std::endl;
        std::cout<<"Make sure your mesh is actually volumetric"<<std::endl;
        exit(-1);
    }

    
    
    vtkSmartPointer<vtkPolyData> surfmesh = vtkSmartPointer<vtkPolyData>::New();
    if(surfmesh_file[strlen(surfmesh_file)-1] == 'k') {
        vtkSmartPointer<vtkPolyDataReader> poly_rdr = vtkSmartPointer<vtkPolyDataReader>::New();
        poly_rdr->SetFileName(surfmesh_file);
        poly_rdr->Update();
        surfmesh->DeepCopy( poly_rdr->GetOutput() );
    }
    else
    {
        vtkSmartPointer<vtkXMLPolyDataReader> poly_rdr = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        poly_rdr->SetFileName(surfmesh_file);
        poly_rdr->Update();
        surfmesh->DeepCopy( poly_rdr->GetOutput() );
    }
    

    //read cell scalars from the surface mesh (labels)
    auto scalars_boundary_id = surfmesh->GetCellData()->GetArray(array_name);
    if ( scalars_boundary_id == nullptr )
    {
        cout << "Surafce mesh has no cellarray with name "<<array_name<<endl;
        throw;
    }
    
    
    vtkSmartPointer<vtkPointLocator> ptloc = vtkSmartPointer<vtkPointLocator>::New();
    ptloc->SetDataSet(volmesh);
    ptloc->BuildLocator();
    
    std::vector<FaceInfo> faceData( surfmesh->GetNumberOfCells() );
    
    //create std::vector with information about the cell and point ids for every label
    //for every cell of the surface mesh
    std::cout<<"Looking for boundaries"<<std::endl;
    int max_nodes_per_face = 3;
    double closestPoint[3];
    int subId;
    double dist2;
    vtkIdType cellid;  
    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();


    for(vtkIdType i=0; i<surfmesh->GetNumberOfCells(); i++)
    {
        if( i%10000 == 0 )
            std::cout<<"Cell "<<i<<"/"<<surfmesh->GetNumberOfCells()<<"\r"<<std::flush;
        
        vtkCell* cell = surfmesh->GetCell(i);
        FaceInfo entry;
        

        entry.boundary_id = scalars_boundary_id->GetTuple1(i);
        const int npoints = cell->GetNumberOfPoints();
        entry.pt_id.resize( npoints );

        if( npoints > max_nodes_per_face ) max_nodes_per_face = npoints;

        vector<vtkIdType> cell_toucing_face; //this will keep an intersection af all the cells, in the end it will be 1 element
        vector<vtkIdType> cells_toucing_face; 

        for( int jj=0; jj<npoints; jj++){
            const vtkIdType ptid = ptloc->FindClosestPoint(cell->GetPoints()->GetPoint(jj));
            if(ptid<0){
                cout<<"Point locator failed for point "<<jj<<" coords "<<cell->GetPoints()->GetPoint(jj)[0]<<" "<<cell->GetPoints()->GetPoint(jj)[1]<<" "<<cell->GetPoints()->GetPoint(jj)[2]<<endl;
                throw;
            }

            entry.pt_id[jj] = ptid;

            //start idenitying the only veolumetric cell touching the face
            volmesh->GetPointCells(ptid, cellIds);
            cells_toucing_face.clear();
            for(int cellId=0; cellId<cellIds->GetNumberOfIds(); cellId++)
                cells_toucing_face.push_back(cellIds->GetId(cellId));


            if( cell_toucing_face.size()==0 ) 
                cell_toucing_face = cells_toucing_face;
            else
                cell_toucing_face = intersection(cell_toucing_face, cells_toucing_face);

        }
                
        if(cell_toucing_face.size()!=1)
        {
            std::cout<<"The number of cells ("<<cell_toucing_face.size()<<"): "<< cell_toucing_face <<", touching face "<<i<<" is not 1. LELBO will be incorrect. Revise the mesh."<<std::endl;
            exit(1);
        }
        entry.vol_cell_id = cell_toucing_face[0];
        
        faceData[i] = entry;
    }
    std::cout<<std::endl;
    

    //cout<<"Correspondence:"<<endl;
    //for(auto &face: faceData ){
    //    cout<<face.boundary_id;
    //    for(auto id: face.pt_id)
    //        cout<<" "<<id;
    //    cout<<face.vol_cell_id<<endl;
    //}


    /*===================================================================================
    *            
    *            write data
    *
    *
    * ===================================================================================*/
    if( Contains(save_formats, VTK) ){
        cout<<"================================================================="<<endl;
        cout<<"Writing VTK"<<endl;
        cout<<"================================================================="<<endl;

        //Save the boundary ids at volmesh elements to verify correct boundary identification
        vtkSmartPointer<vtkShortArray> volmesh_boundaries = vtkSmartPointer<vtkShortArray>::New();
        volmesh_boundaries->SetName(array_name);
        volmesh_boundaries->SetNumberOfComponents(1);
        volmesh_boundaries->SetNumberOfValues(volmesh->GetNumberOfCells());


        //for every cell of the surface mesh
        for( auto& face: faceData ) {
            volmesh_boundaries->SetTuple1( face.vol_cell_id, face.boundary_id );
        }
        
        volmesh->GetCellData()->AddArray(volmesh_boundaries);
        
        string filename = string(outfile_prefix)+".vtu";

        vtkSmartPointer<vtkTransform> T = vtkSmartPointer<vtkTransform>::New();
        vtkSmartPointer<vtkTransformFilter> tf = vtkSmartPointer<vtkTransformFilter>::New();
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> wrwr = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        CommonTools::AssociateProgressFunction(wrwr);

        wrwr->SetInputData(volmesh);

        if(scale!=1) {
            cout<<"Rescaling mesh"<<endl;
            T->Scale(scale, scale, scale);

            tf->SetTransform(T);
            tf->SetInputData(volmesh);
            tf->Update();

            wrwr->SetInputData(tf->GetOutput());
        }

        wrwr->SetDataModeToAppended ();
        wrwr->EncodeAppendedDataOff();
        wrwr->SetFileName(filename.c_str());
        wrwr->Write();
    }

    if( Contains(save_formats, ALYATXT) ){//needs checking. Looks fishy
        cout<<"================================================================="<<endl;
        cout<<"Writing ALYA ASCII"<<endl;
        cout<<"================================================================="<<endl;

        string filename_BC = string(outfile_prefix)+".BC_NBOUN.in";
        string filename_LELBO = string(outfile_prefix)+".LELBO.in";
        string filename_BOUNDARIES = string(outfile_prefix)+".BOUNDARIES.in";

        //store the labels
        std::ofstream file_BC(filename_BC.c_str());
        std::ofstream file_LELBO(filename_LELBO.c_str());
        std::ofstream file_BOUNDARIES(filename_BOUNDARIES.c_str());

        cout<<"Writing boundaries, LELBO, BC"<<endl;
        for( vtkIdType face_id = 0; face_id < faceData.size(); face_id++ )
        {            
            const auto& face = faceData.at(face_id);

            file_BOUNDARIES << face_id+1;
            for( int jj=0; jj<face.pt_id.size(); jj++ ){
                file_BOUNDARIES<< " " << face.pt_id[jj]+1;
            }
            file_BOUNDARIES << endl;


            file_LELBO << face_id+1 <<" "<< face.vol_cell_id+1<<endl;
            

            file_BC <<  face_id+1 << " " << face.boundary_id << endl;                    
        }        

        std::cout<<"Writing alya ASCII mesh"<<std::endl;

        CommonTools::SaveVolMeshBSC(volmesh, outfile_prefix, scale, correct_orientation);      
    }

    if( Contains(save_formats, ALYAMPIO) ){
        cout<<"================================================================="<<endl;
        cout<<"Writing ALYA MPIO"<<endl;
        cout<<"================================================================="<<endl;

        {
            cout<<"Saving "<<volmesh->GetNumberOfPoints()<<" nodes"<<endl;
            MatrixRXd nodes(volmesh->GetNumberOfPoints(),3);
            for( vtkIdType ptid =0; ptid<nodes.rows(); ptid++ ) {
                double pt[3];
                volmesh->GetPoint(ptid, pt);
                nodes( ptid, 0 ) = pt[0]*scale;
                nodes( ptid, 1 ) = pt[1]*scale;
                nodes( ptid, 2 ) = pt[2]*scale;
            }
            WriteMPIOMatrix( nodes, outfile_prefix, "COORD", "NPOIN", "REAL" );
        }

        {
            cout<<"Saving LNODS and LTYPE"<<endl;
            //# TET04 30
            //# PYR05 32
            //# PEN06 34
            //# HEX08 37

            cout<<"Finding cell with largest number of nodes"<<endl;
            int max_nodes_per_solid_cell = 4; //default tetra
            MatrixRXi nodetypes(volmesh->GetNumberOfCells(),1);

            for( vtkIdType i=0; i<volmesh->GetNumberOfCells(); i++ ){
                const auto cell = volmesh->GetCell(i);


                if(cell->GetNumberOfPoints()==4) { //TET04
                    nodetypes(i,0) = 30; //TET04
                }
                else if(cell->GetNumberOfPoints()==8){ //HEX08
                    max_nodes_per_solid_cell = max(8, max_nodes_per_solid_cell);
                    nodetypes(i,0) = 37; //HEX08
                }
                else if(cell->GetNumberOfPoints()==5) { //PYR05
                    max_nodes_per_solid_cell = max(5, max_nodes_per_solid_cell);
                    nodetypes(i,0) = 32; //PYR05
                }
                else if(cell->GetNumberOfPoints()==6) {//PEN06
                    max_nodes_per_solid_cell = max(6, max_nodes_per_solid_cell);
                    nodetypes(i,0) = 34; //PEN06
                    max_nodes_per_face = 4; //ALya works in misterios ways, it expects 4 columns always if there are pen06 elements
                }
                else
                {
                    cout<<"Volume element "<<i<<" has "<<cell->GetNumberOfPoints()<<" vertices, unsupported"<<endl;
                    throw;
                }
                
            }

            cout<<"Extracting elements"<<endl;
            MatrixRXi nodes = MatrixRXi::Zero(volmesh->GetNumberOfCells(), max_nodes_per_solid_cell);
            vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

            set<int> unique_elements;

            for( vtkIdType i=0; i<volmesh->GetNumberOfCells(); i++ ){
                volmesh->GetCellPoints(i, ptIds);	
                for( int jj=0; jj<ptIds->GetNumberOfIds(); jj++ ){
                    nodes(i, jj) = ptIds->GetId(jj)+1;
                }

                unique_elements.insert(ptIds->GetNumberOfIds());
            }

            cout<<"Writing LTYPE (Element types) on elems"<<endl;
            WriteMPIOMatrix( nodes, outfile_prefix, "LNODS", "NELEM", "INTEG" );
            cout<<"Writing LNODS (Elements connectivity) on elems"<<endl;
            WriteMPIOMatrix( nodetypes, outfile_prefix, "LTYPE", "NELEM", "INTEG" );


            //write info
            string info_filename = string(outfile_prefix) + "-INFO.in";
            ofstream info_file(info_filename);
            info_file << "  TYPES_OF_ELEMENTS = ";
            for( auto value: unique_elements ){
                switch (value) {
                    case 4:
                        info_file<<"TET04, ";
                        break; 
                    case 5:
                        info_file<<"PYR05, ";
                        break;
                    case 6: 
                        info_file<<"PEN06, ";
                        break;
                    case 8:
                        info_file<<"HEX08, ";
                        break;
                }        
            }
            info_file<<endl;

            info_file << "  BOUNDARIES   = "<< faceData.size() <<endl;
            info_file << "  ELEMENTS     = "<< volmesh->GetNumberOfCells() <<endl;
            info_file << "  NODAL_POINTS = "<< volmesh->GetNumberOfPoints() <<endl;

        }


        {
            //# ====================================
            //#
            //#          save LNODB, LELBO, CODBO on NBOUN and LBSET
            //#
            //# ====================================

            MatrixRXi nodes = MatrixRXi::Zero(faceData.size(), max_nodes_per_face); 
            MatrixRXi boundary_ids( faceData.size(), 1 );
            MatrixRXi LELBO( faceData.size(), 1 );

            cout<<"Extracting LNODB, LELBO, CODBO"<<endl;
            for( vtkIdType face_id=0; face_id<faceData.size(); face_id++ ){
                const auto& face = faceData.at(face_id);

                for( int jj=0; jj<face.pt_id.size(); jj++ ){
                    nodes(face_id,jj) = face.pt_id[jj]+1;
                }

                LELBO(face_id,0) = face.vol_cell_id+1;
                boundary_ids(face_id,0) = face.boundary_id;
            }

            cout<<"Writing LNODB (boundary faces)"<<endl;
            WriteMPIOMatrix( nodes, outfile_prefix, "LNODB", "NBOUN", "INTEG" );

            cout<<"Writing LELBO (Boundary elements) volume cell ids, NBOUN, LNODB, scalar"<<endl;
            WriteMPIOMatrix( LELBO, outfile_prefix, "LELBO", "NBOUN", "INTEG" );
            
            cout<<"Writing CODBO on NBOUN"<<endl;
            WriteMPIOMatrix( boundary_ids, outfile_prefix, "CODBO", "NBOUN", "INTEG" );
            cout<<"Writing LBSET on NBOUN"<<endl;
            WriteMPIOMatrix( boundary_ids, outfile_prefix, "LBSET", "NBOUN", "INTEG" );
        }
    }

    if( Contains(save_formats, ELMER) ){
        //needs checking. Looks fishy
        cout<<"Elmer format needs verification. Ignoring."<<endl;
        /*
        std::cout<<"Writing elmer file"<<std::endl;
        std::string header_filename = "mesh.header";
        std::ofstream header_file(header_filename.c_str());
        header_file<<volmesh->GetNumberOfPoints()<<" "<<volmesh->GetNumberOfCells()<<" "<<faceData.size()<<std::endl;
        header_file<<"2"<<std::endl;
        header_file<<"303 "<<faceData.size()<<std::endl;
        header_file<<"504 "<<volmesh->GetNumberOfCells()<<std::endl;

        std::string node_filename = "mesh.nodes";
        std::ofstream node_file(node_filename.c_str());
        for(vtkIdType i=0; i<volmesh->GetNumberOfPoints(); i++)
        { MatrixRX<T>& m, const char* problem_name, const char *variable, const char* association, const char* dtype
            double* pt = volmesh->GetPoint(i);
            node_file<<i+1<<" -1 "<<pt[0]*scale<<" "<<pt[1]*scale<<" "<<pt[2]*scale<<std::endl;
        }

        std::string ele_filename = "mesh.elements";
        std::ofstream ele_file(ele_filename.c_str());
        for(vtkIdType i=0; i<volmesh->GetNumberOfCells(); i++)
        {
            vtkCell* cell = volmesh->GetCell(i);
            ele_file<<i+1<<" 1 504";
            for(int k=0; k<cell->GetNumberOfPointsMatrix<int, 3, 4, ColMajor> Acolmajor;(); k++)                
                ele_file<<" "<<cell->GetPointId(k)+1;
            ele_file<<std::endl;
        }
        
        
        std::string bound_filename =  "mesh.boundary";
        std::ofstream bound_file(bound_filename.c_str());
        for(vtkIdType i=0; i<labeldata.size(); i++)
        {
            BscEntry entry = faceData[i];
            bound_file<<i+1<<" "<<
                    entry.id<<" "<<
                    entry.vol_cell_id+1<<" 0 303 "<<
                    entry.pt1+1<<" "<<
                    entry.pt2+1<<" "<<
                    entry.pt3+1;
            if( entry.nvertices==4 )
                bound_file<<" "<<entry.pt4+1;

            bound_file<<std::endl;
        }        
        */
    }

    return 0;
}




template <class T>
void WriteMPIOMatrix( MatrixRX<T>& m, const char* problem_name, const char *variable, const char* association, const char* dtype )
{
    //m is Npoints x ndim 
    //variable = COORD,  LNODS, LTYPE, LTYPE, LNODS, LMATE, LNODB, LELBO
    //association = NPOIN, NELEM, NBOUN
    //datatype = REAL, INTEG

    int64_t ndim = m.cols();
    int64_t npts = m.rows();

    char mpio_var[8] = "0000000";
    char mpio_assoc[8] = "0000000";
    char mpio_dtype[8] = "0000000";
    for( int i=0; i<min(strlen(variable), (size_t)7); i++ ) mpio_var[i] = variable[i];
    for( int i=0; i<min(strlen(association), (size_t)7); i++ ) mpio_assoc[i] = association[i];
    for( int i=0; i<min(strlen(dtype), (size_t)7); i++ ) mpio_dtype[i] = dtype[i];


    string filename = string(problem_name)+"-"+ string(variable) + ".mpio.bin";

    std::ofstream file(filename.c_str(), ios::binary);
    int64_t file_id = 27093;
    file.write( reinterpret_cast<const char *>(&file_id), sizeof(file_id));
    file.write("MPIAL00", 8);
    file.write("V000400", 8);
    file.write( mpio_var, 8);

    if( ndim == 1 ) {
        file.write("SCALA00", 8);
    }
    else {
        file.write("VECTO00", 8);
    }

    file.write( mpio_assoc, 8 );
    file.write( mpio_dtype, 8 );
        
    file.write("8BYTE00", 8);
    file.write("SEQUE00", 8);
    file.write("NOFIL00", 8);
    file.write("ASCEN00", 8);
    file.write("NOID000", 8);
    file.write("0000000", 8);

    file.write( reinterpret_cast<const char *>(&ndim), sizeof(ndim));
    file.write( reinterpret_cast<const char *>(&npts), sizeof(npts));

    int64_t num = 0;
    file.write( reinterpret_cast<const char *>(&num), sizeof(num)); // Time step number (signed int 64 bits):    ittim
    num = 1;
    file.write( reinterpret_cast<const char *>(&num), sizeof(num)); // n of subdomains (signed int 64 bits):nsubd (1=SEQUENTIAL)
    num = 0;
    file.write( reinterpret_cast<const char *>(&num), sizeof(num)); // Mesh division (signed int 64 bits):       divi
    file.write( reinterpret_cast<const char *>(&num), sizeof(num)); // Tag 1 (signed int 64 bits):               tag1
    file.write( reinterpret_cast<const char *>(&num), sizeof(num)); // Tag 2 (signed int 64 bits):               tag2
    double num1 = 0;
    file.write( reinterpret_cast<const char *>(&num1), sizeof(num1)); // Time (real 64 bits):                      time


    file.write("0000000", 8);
    file.write("NONE000", 8); //#1
    file.write("NONE000", 8); //#2
    file.write("NONE000", 8); //#3
    file.write("NONE000", 8); //#4
    file.write("NONE000", 8); //#5
    file.write("NONE000", 8); //#6
    file.write("NONE000", 8); //#7
    file.write("NONE000", 8); //#8
    file.write("NONE000", 8); //#9
    file.write("NONE000", 8); //#10

    file.write((char *) m.data(), m.rows() * m.cols() * sizeof(m(0,0)));
}
