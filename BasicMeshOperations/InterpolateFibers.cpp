#include <vtkDataSet.h>
#include <vtkDataSetReader.h>
#include <vtkDataArray.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkProbeFilter.h>
#include <vtkShortArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>

#include <VTKCommonTools.h>


//https://stackoverflow.com/questions/236129/the-most-elegant-way-to-iterate-the-words-of-a-string
std::vector<std::string> split(const std::string& s, const std::string& delim, const bool keep_empty = true) {
    std::vector<std::string> result;
    if (delim.empty()) {
        result.push_back(s);
        return result;
    }
    std::string::const_iterator substart = s.begin(), subend;
    while (true) {
        subend = search(substart, s.end(), delim.begin(), delim.end());
        std::string temp(substart, subend);
        if (keep_empty || !temp.empty()) {
            result.push_back(temp);
        }
        if (subend == s.end()) {
            break;
        }
        substart = subend + delim.size();
    }
    return result;
}



int main(int argc, char **argv)
{    
    if(argc<6)
    {
        std::cout<<"Usage: InterpolateFibers unstructured_grid.vtk fibre_array labels.txt target_points.txt out_fibers.txt out_labels.txt out.vtk"<<std::endl;
        std::cout<<"Given: "<<std::endl<<
                "1) source tetrahedral mesh in vtk format"<<std::endl<<
                "2) fibers associated to the points of the source mesh in pointdata"<<std::endl<<
                "3) label associated to the points of the source mesh in a text file"<<std::endl<<
                "4) Points at which you need to calculate the fibers and labels in a text file"<<std::endl<<
                "Return:"<<std::endl<<
                "1) text file with interpolated fibers (linear interpolation)"<<std::endl<<
                "2) text file with interpolated labels (nearest neighbor interpolation)"<<std::endl<<std::endl;
        std::cout<<"Text files have format: id value(s). FOr example for fibers:"<<std::endl;
        std::cout<<"1 -0.714724 0.696612 -0.062452"<<std::endl<<
                   "2 -0.787582 0.616203 -0.003063"<<std::endl;
        return -1;
    }
    
    int c=1;
    const char *ug_filename = argv[c++];
    const char *fiberarray_name = argv[c++];
    const char *labels_filename = argv[c++];
    const char *points_filename = argv[c++];
    const char *out_fibers_filename = argv[c++];
    const char *out_labels_filename = argv[c++];
    const char *outvtk_filename = argv[c++];
    
    
    std::cout<<"Mesh: "<<ug_filename<<std::endl;
    std::cout<<"Fiber array: "<<fiberarray_name<<std::endl;
    std::cout<<"Labels filename: "<<labels_filename<<std::endl;
    std::cout<<"Target points: "<<points_filename<<std::endl;
    std::cout<<"Output for fibers: "<<out_fibers_filename<<std::endl;
    std::cout<<"Output for labels: "<<out_labels_filename<<std::endl;
    std::cout<<"VTK output with the result: "<<outvtk_filename<<std::endl;
    
    
    std::cout<<"Reading the mesh"<<std::endl;
    vtkSmartPointer<vtkDataSetReader> mesh_rdr =     vtkSmartPointer<vtkDataSetReader>::New();
    mesh_rdr->SetFileName(ug_filename);
    CommonTools::AssociateProgressFunction(mesh_rdr);
    mesh_rdr->Update();
    
    vtkUnstructuredGrid *mesh = vtkUnstructuredGrid::SafeDownCast(mesh_rdr->GetOutput());
    
    auto fibers = mesh->GetPointData()->GetArray(fiberarray_name);
    if(fibers==nullptr)
    {
        std::cout<<"FIber array "<<fiberarray_name<<" does not exist. Aborting!"<<std::endl;
        exit(-1);
    }   
        
    std::cout<<"Reading labels"<<std::endl;
    vtkSmartPointer<vtkShortArray> lbl =     vtkSmartPointer<vtkShortArray> ::New();
    lbl->SetName("Labels");
    lbl->SetNumberOfComponents(1);
    lbl->SetNumberOfTuples(mesh->GetNumberOfPoints());
    
    std::ifstream file;
    std::string line;
    file.open(labels_filename);
    if (!file.is_open())
    {
        std::cout<<"Cannot open file"<<std::endl;
        exit(-1);
    }

    for(vtkIdType i=0; i<mesh->GetNumberOfPoints(); i++)
    {
        if(i%10000==0)
            std::cout<<"Point "<<i<<"/"<<mesh->GetNumberOfPoints()<<"\r"<<std::flush;
        
        getline (file, line);
        auto values = split(line, " ", false);
        lbl->SetTuple1(i, atoll(values[1].c_str()) );
    }
    file.close();
    std::cout<<std::endl;
    
    mesh->GetPointData()->AddArray(lbl);
    
    //Interpolating fibers
    std::cout<<"Reading points"<<std::endl;
    vtkSmartPointer<vtkPoints> pts =     vtkSmartPointer<vtkPoints> ::New();
    file.open(points_filename);
    if (!file.is_open())
    {
        std::cout<<"Cannot open file"<<std::endl;
        exit(-1);
    }

    
    while ( getline (file,line) )
    {
        //std::cout<<"Line "<<line<<std::endl;
        
        auto values = split(line, " ", false);
        if(values.size()>0)
        {
            //std::cout<<"Values "<<values[0]<<"  "<<values[1]<<"  "<<values[2]<<"  "<<values[3]<<std::endl;

            auto id = pts->InsertNextPoint( atof(values[1].c_str()), atof(values[2].c_str()), atof(values[3].c_str()) );
            //std::cout<<"Id :"<<id<<std::endl;
            //std::cout<<"Inserted point"<<pts->GetPoint(id)[0]<<" "<<pts->GetPoint(id)[1]<<" "<<pts->GetPoint(id)[2]<<std::endl;

            if(id%10000==0)
                std::cout<<"Reading point "<<id<<"\r"<<std::flush;
        }
    }
    file.close();
    std::cout<<std::endl;

    vtkSmartPointer<vtkPolyData> pd =     vtkSmartPointer<vtkPolyData> ::New();
    pd->SetPoints(pts);
    
    //interpolation
    std::cout<<"Interpolating labels"<<std::endl;
    vtkSmartPointer<vtkProbeFilter> probe =     vtkSmartPointer<vtkProbeFilter> ::New();
    CommonTools::AssociateProgressFunction(probe);
    probe->SetInputData(pd);
    probe->SetSourceData(mesh);
    probe->Update();
    
    auto probed = probe->GetOutput();
    
    //saving
    std::cout<<"Saving fibers and labels"<<std::endl;
    vtkSmartPointer<vtkPolyDataWriter> pdwr =     vtkSmartPointer<vtkPolyDataWriter> ::New();
    pdwr->SetFileName(outvtk_filename);
    pdwr->SetInputData(probed);
    pdwr->Write();
    
    std::ofstream outfile_fibers;
    std::ofstream outfile_labels;
    outfile_fibers.open(out_fibers_filename);
    for(vtkIdType i=0; i<probed->GetNumberOfPoints(); i++)       
    {
        if(i%10000==0)
            std::cout<<"Saving id "<<i<<"/"<<probed->GetNumberOfPoints()<<"\r"<<std::flush;

        double *fiber = probed->GetPointData()->GetArray(fiberarray_name)->GetTuple(i);
        double length = sqrt(fiber[0]*fiber[0]+fiber[1]*fiber[1]+fiber[2]*fiber[2]);
        outfile_fibers<<i+1<<" "<<fiber[0]/length<<" "<<fiber[1]/length<<" "<<fiber[2]/length<<std::endl;
        
        auto label = probed->GetPointData()->GetArray("Labels")->GetTuple1(i);
        outfile_labels<<i+1<<" "<<label<<std::endl;
    }
    outfile_fibers.close();
    outfile_labels.close();
    std::cout<<std::endl;
    
    
    
    return 0;
}
