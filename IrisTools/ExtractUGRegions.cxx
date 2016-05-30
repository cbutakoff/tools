#include "vtkCylinderSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGridReader.h" 
#include "vtkUnstructuredGrid.h" 
#include "vtkThreshold.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkCellData.h"
#include "vtkShortArray.h"

#include <iostream>

int main(int argc, char** argv) {
    vtkSmartPointer<vtkUnstructuredGridReader> rd = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    rd->SetFileName(argv[1]);
    rd->Update();

    std::cout << "Npoints :" << rd->GetOutput()->GetNumberOfPoints() << std::endl;
    std::cout << "Ncells :" << rd->GetOutput()->GetNumberOfCells() << std::endl;

    vtkUnstructuredGrid* mesh = rd->GetOutput();
    mesh->GetCellData()->SetActiveScalars("Material");
    
    vtkShortArray* ids = dynamic_cast<vtkShortArray*>(mesh->GetCellData()->GetArray("Material"));
    
    vtkSmartPointer<vtkThreshold> th = vtkSmartPointer<vtkThreshold>::New();
    th->SetInputData(mesh);
    th->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "Material");
    
    std::cout<<"Ids: "<<ids->GetNumberOfTuples()<<std::endl;
    double range[2];
    ids->GetRange(range);
    std::cout<<"Ids: "<<range[0]<<"-"<<range[1]<<std::endl;
    
    
    for(vtkIdType i=range[0]; i<=range[1]; i++)
    {
        std::cout<<"Id: "<<i<<std::endl;
        th->ThresholdBetween(i, i);
        th->Update();
        vtkUnstructuredGrid* submesh = th->GetOutput();
        
        
        vtkSmartPointer<vtkUnstructuredGridWriter> wr = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
        char filename[100];
        sprintf(filename,"region%03lld.vtk",i);
        wr->SetFileName(filename);
        wr->SetInputData(submesh);
        wr->SetFileTypeToBinary();
        wr->Write();
    }
    
    
    return 0;
}
