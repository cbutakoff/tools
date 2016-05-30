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
#include <vtkCell.h>
#include <vector>
#include <vtkCellArray.h>
#include <vtkType.h>
#include <vtkCellData.h>
#include <vtkTetra.h>

int main(int argc, char **argv)
{
    std::cout<<"Removes elements that are not tetrahedra v1.0"<<std::endl;
    if(argc<2)
    {
        std::cout<<"Usage: CheckTetraMesh -i inputmesh.vtk -o outputmessh.vtk -clean -flip(p|n)"<<std::endl;
        std::cout<<"-clean -- removes non-tetrahedra"<<std::endl;
        std::cout<<"-flipp -- flips tetra (ABCD) with positive dot product (ABxAC).AD"<<std::endl;
        std::cout<<"-flipn -- flips tetra (ABCD) with negative dot product (ABxAC).AD"<<std::endl;
        std::cout<<"-info -- diagnose orientation"<<std::endl;
        return -1;
    }
    
    if(argc<=1) exit(-1);
    
    char *input_volmesh_filename = NULL;
    char *output_volmesh_filename = NULL;
    bool check_orientation = false;
    bool flipp = false;
    bool flipn = false;
    bool clean = false;
    

    for (int c = 1; c < argc; c++) {
        if (strcmp(argv[c], "-i") == 0) 
            input_volmesh_filename = argv[++c];
        if (strcmp(argv[c], "-o") == 0) 
            output_volmesh_filename = argv[++c];    
        if (strcmp(argv[c], "-flipp") == 0) 
        {
            flipp = true;    
            flipn = false;    
        }
        if (strcmp(argv[c], "-flipn") == 0) 
        {
            flipp = false;    
            flipn = true;    
        }
        if (strcmp(argv[c], "-info") == 0) 
            check_orientation = true;
        if (strcmp(argv[c], "-clean") == 0) 
            clean = true;
    }

    std::cout<<"Input mesh: "<<input_volmesh_filename<<std::endl;
    if(output_volmesh_filename!=NULL)
        std::cout<<"Output mesh: "<<output_volmesh_filename<<std::endl;
    
    std::cout<<"Requested operations:"<<std::endl;
    if(clean)
        std::cout<<"Remove non-tetrahedra"<<std::endl;
    if(check_orientation)
        std::cout<<"Print orientation info"<<std::endl;
    if(flipp)
        std::cout<<"Flipping tetra with positive dot product (ABxAC).AD"<<std::endl;
    if(flipn)
        std::cout<<"Flipping tetra with negative dot product (ABxAC).AD"<<std::endl;

    std::cout<<std::endl;    
    
    //break if no file supplied
    if(input_volmesh_filename==NULL)
        return -1;
    
    //
    //
    //  Processing starts here
    //
    //
    
    if(flipp||flipn||clean)
        if(output_volmesh_filename==NULL)
        {
            std::cout<<"No output filename supplied."<<std::endl;
            return -1;
        }
    
    std::cout<<"Loading mesh"<<std::endl;
    vtkSmartPointer<vtkDataSetReader> rdr = vtkSmartPointer<vtkDataSetReader>::New();
    rdr->SetFileName(input_volmesh_filename);
    rdr->Update();
    
    vtkUnstructuredGrid *mesh = (vtkUnstructuredGrid *)rdr->GetOutput();
    
    if(clean)
    {
        std::cout<<"Cleaning mesh"<<std::endl;
        vtkSmartPointer<vtkCellArray> array =vtkSmartPointer<vtkCellArray>::New();

        //verify that all the cells are tetrahedra
        //and create an indicator array of point ids belonging to cells
        //std::vector<int> indic(mesh->GetNumberOfPoints,0);

        for(vtkIdType i=0; i<mesh->GetNumberOfCells(); i++)
        {
            if(mesh->GetCell(i)->GetNumberOfPoints()!=4)
            {
                std::cout<<"Cell "<<i<<" is not tetrahedron. Number of vertices: "<<mesh->GetCell(i)->GetNumberOfPoints()<<std::endl;
            }
            else
                array->InsertNextCell(mesh->GetCell(i));
        }


        mesh->SetCells(VTK_TETRA,array);
        mesh->GetCellData()->Initialize();
    }
     
    
    //
    //
    //
    //
    //
    //
    std::vector<char> positive_dp(mesh->GetNumberOfCells());
    
    if(check_orientation || flipp || flipn)
    {
        std::cout<<"Checking tetra orientation. For (ABCD) tetra it is (AB x AC).AD:"<<std::flush<<std::endl;
        //test tetra orientation
        //formula for (ABCD) tetra is (AB x AC).AD
        unsigned long dp_positive=0; 
        unsigned long dp_negative=0; 

        for(vtkIdType i=0; i<mesh->GetNumberOfCells(); i++){
            double A[3], B[3], C[3], D[3];
            mesh->GetCell(i)->GetPoints()->GetPoint(0,A);
            mesh->GetCell(i)->GetPoints()->GetPoint(1,B);
            mesh->GetCell(i)->GetPoints()->GetPoint(2,C);
            mesh->GetCell(i)->GetPoints()->GetPoint(3,D);


            double AB[3], AC[3], AD[3];
            for(int j=0; j<3; j++) 
            {
                AB[j]=B[j]-A[j];
                AC[j]=C[j]-A[j];
                AD[j]=D[j]-A[j];
            }

            //cross product
            double ABxAC[3]; 
            ABxAC[0] = AB[1]*AC[2] - AB[2]*AC[1];
            ABxAC[1] = AB[2]*AC[0] - AB[0]*AC[2];
            ABxAC[2] = AB[0]*AC[1] - AB[1]*AC[0];

            //dot product
            const double dp = ABxAC[0]*AD[0]+ABxAC[1]*AD[1]+ABxAC[2]*AD[2];

            if(dp>=0) 
            {
                positive_dp[i]=1;
                dp_positive++;
            }
            else 
            {
                positive_dp[i]=0;
                dp_negative++;
            }
        }

        std::cout<<"Number positive dot products: "<<dp_positive<<std::endl;
        std::cout<<"Number negative dot products: "<<dp_negative<<std::endl;
    }
    
    
    //
    //  Flip tetrahedra
    //
    //
    if(flipp || flipn)
    {
        std::cout<<"Flipping"<<std::flush<<std::endl;

        vtkSmartPointer<vtkCellArray> array =vtkSmartPointer<vtkCellArray>::New();
                
        for(vtkIdType i=0; i<mesh->GetNumberOfCells(); i++){
            
            if( (positive_dp[i]==1 && flipp) || (positive_dp[i]==0 && flipn))
            {
                const vtkIdType A_id = mesh->GetCell(i)->GetPointId(0);
                const vtkIdType B_id = mesh->GetCell(i)->GetPointId(1);
                const vtkIdType C_id = mesh->GetCell(i)->GetPointId(2);
                const vtkIdType D_id = mesh->GetCell(i)->GetPointId(3);
                //std::cout<<"Order original:"<<A_id<<" "<<B_id<<" "<<C_id<<" "<<D_id<<std::endl;
                array->InsertNextCell(4); //reverse orientation
                array->InsertCellPoint(D_id);
                array->InsertCellPoint(B_id);
                array->InsertCellPoint(C_id);
                array->InsertCellPoint(A_id);
            }
            else
                array->InsertNextCell(mesh->GetCell(i));
            
        }
        
        mesh->SetCells(VTK_TETRA,array);
        mesh->GetCellData()->Initialize();
    }
    
    if(output_volmesh_filename!=NULL)
    {
        std::cout<<"Saving mesh"<<std::endl;
        vtkSmartPointer<vtkDataSetWriter> wr = vtkSmartPointer<vtkDataSetWriter>::New();
        wr->SetFileName(output_volmesh_filename);
        wr->SetInputData(mesh);
        wr->Write();
    }
    
    
    return 0;
}
