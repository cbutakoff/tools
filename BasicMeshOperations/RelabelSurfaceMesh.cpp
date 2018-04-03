/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   RelabelSurfaceMesh.cpp
 * Author: costa
 *
 * Created on April 2, 2018, 3:43 PM
 */

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>

#include <cstdlib>
#include <map>

using namespace std;


std::map<vtkIdType, vtkIdType> ReadLabelsFromCSV(const char* filename)
{
    std::map<vtkIdType, vtkIdType> mapping;
    
    string line;
    ifstream f(filename);
    //getline(f,line); //read header
    
    char c;
    vtkIdType idold, idnew;
    while (f >> idold >> c >> idnew )
    {
        cout<<"Convert:"<<idold<<" -> "<<idnew<<endl;
        mapping[idold]= idnew;
    }
    
    return mapping;
}

/*
 * 
 */
int main(int argc, char** argv) {

    if(argc<4)
    {
        cout<<"Usage: RelabelSurfaceMesh mesh.vtk cell_array_name label_correspondece.csv mesh_out.vtk"<<endl;
        cout<<"label_correspondece.csv is a comma separated file, without header:"<<endl;
        cout<<"each line is of the form: old_id, new_id"<<endl;
        return -1;
    }
    
    int c=1;
    const char* inmeshfilename = argv[c++];
    const char* cellarray = argv[c++];
    const char* label_filename = argv[c++];
    const char* outmeshfilename = argv[c++];
    
    cout<<"Input Mesh :"<<inmeshfilename<<endl;
    cout<<"Output mesh:"<<outmeshfilename<<endl;
    cout<<"Labels file:"<<label_filename<<endl;
    cout<<"Labels array:"<<cellarray<<endl;
    
    cout<<"Reading mesh"<<endl;
    vtkSmartPointer<vtkPolyDataReader> rd =     vtkSmartPointer<vtkPolyDataReader> ::New();
    rd->SetFileName(inmeshfilename);
    rd->Update();
    vtkPolyData* mesh = rd->GetOutput();
    
    auto array = mesh->GetCellData()->GetArray(cellarray);
    if(array==nullptr)
    {
        cout<<"NULL cell array"<<endl;
        return -1;
    }
    
    
    auto label_map = ReadLabelsFromCSV(label_filename);
    
    for(vtkIdType cellid = 0; cellid<mesh->GetNumberOfCells(); cellid++)
    {
        vtkIdType id = array->GetTuple1(cellid);
        
        auto it=  label_map.find(id);
        if(it!=label_map.end())
        {
            array->SetTuple1(cellid, it->second);
        }
    }
    
    vtkSmartPointer<vtkPolyDataWriter> wr =     vtkSmartPointer<vtkPolyDataWriter> ::New();
    wr->SetFileName(outmeshfilename);
    wr->SetInputData(mesh);
    wr->SetFileTypeToBinary();
    wr->Write();
    
    
    
    return 0;
}

