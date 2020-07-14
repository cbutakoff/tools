/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
#ifndef __vtkTestProgressReportFilter_h
#define __vtkTestProgressReportFilter_h

#include <vtkPolyDataAlgorithm.h>

class vtkTestProgressReportFilter : public vtkPolyDataAlgorithm
{
public:
  static vtkTestProgressReportFilter *New();
  vtkTypeMacro(vtkTestProgressReportFilter,vtkAlgorithm);

protected:
  vtkTestProgressReportFilter(){}
  ~vtkTestProgressReportFilter() {}

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

private:
  vtkTestProgressReportFilter(const vtkTestProgressReportFilter&) = delete;
  void operator=(const vtkTestProgressReportFilter&) = delete;

};

#endif


