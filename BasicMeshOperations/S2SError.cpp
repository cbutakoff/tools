/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file 
    \brief Calculate S2S error.
    */

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>

#include "VTKCommonTools.h"
#include "CommonTools.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>



int main(int argc, char **argv)
{
	std::cout << "S2SError. Version 1.0"<< std::endl;

	if(argc<4)
	{
		std::cout << "S2SError shape1.vtk shape2.vtk outshape1.vtk outshape2.vtk"<< std::endl;		
		return 0;
	}

        int c=1;
	const char* inp1 = argv[c++];
	const char* inp2 = argv[c++];
	const char* outp1 = argv[c++];
	const char* outp2 = argv[c++];
	
        vtkSmartPointer<vtkPolyData> shape1 = vtkSmartPointer<vtkPolyData>::Take(CommonTools::LoadShapeFromFile(inp1));
        vtkSmartPointer<vtkPolyData> shape2 = vtkSmartPointer<vtkPolyData>::Take(CommonTools::LoadShapeFromFile(inp2));
        
        vtkSmartPointer<vtkCleanPolyData> cl1 = vtkSmartPointer<vtkCleanPolyData>::New();
        vtkSmartPointer<vtkCleanPolyData> cl2 = vtkSmartPointer<vtkCleanPolyData>::New();
        
        cl1->SetInputData(shape1);
        cl1->Update();
        
        cl2->SetInputData(shape2);
        cl2->Update();
        
        
        vtkSmartPointer<vtkPolyData> outShape1 = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPolyData> outShape2 = vtkSmartPointer<vtkPolyData>::New();
        
        double mean, stdev, maxerror;
        
        CommonTools::GetS2S(cl1->GetOutput(), cl2->GetOutput(), outShape1.GetPointer(), outShape2.GetPointer(), mean, stdev, maxerror);
 
        CommonTools::SaveShapeToFile(outShape1.GetPointer(), outp1);
        CommonTools::SaveShapeToFile(outShape2.GetPointer(), outp2);

        std::cout<<"S2S Error: "<<mean<<" +- "<<stdev<<"; max = "<<maxerror<<std::endl;

              
        //outShape1->Delete();
        //outShape2->Delete();
	return 0;
}

