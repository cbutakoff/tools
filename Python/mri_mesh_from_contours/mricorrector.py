#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 09:58:49 2013
Needs: VTK, Numpy

@author: Constantine Butakoff
"""

import vtk
import numpy as np
import mytools
#


def CreateBottomText(displacement, saveflag=False):
    """
    Creates message at the bottom of the render window
    
    displacement float: slice displacement value to print
    saveflag bool: print message if save was requested or not

    Return nothing
    """
    s = 'e - finish; s - save; 1/2 - toggle visibility; P/; - change slice; [/\' - change displ; arrows - move slice; \nDisplacement: '+str(displacement)
    if saveflag:
        s = s+'; Save requested'
    return s
    

def Points2Shape(pts):
    """
    Creates vtkPolyData from points and adds array 'colors' for slice highlighting.
    Assigns 0 to every point
    
    pts Nx3 numpy array of points

    Return vtkPolyData
    """
    shape = mytools.VTKPoints2PolyData( pts )
    scalars = vtk.vtkShortArray()
    scalars.SetName('colors')
    scalars.SetNumberOfComponents(1)
    scalars.SetNumberOfTuples(shape.GetNumberOfPoints())
    
    for i in range(scalars.GetNumberOfTuples()):
        scalars.SetValue(i,0)
        
    shape.GetPointData().AddArray(scalars)
    shape.GetPointData().SetActiveScalars('colors')    
    return shape
    
def ChangeLayerColor(polydata, layerids, oldlayerid, layerid):
    """
    Finds points belonging to layerid using layerids indicator array,
    assigns them color 1
    
    polydata vtkPolyData: the shape
    layerids 1xN numpy array: stores the slice number for every point
    oldlayerid int: the layer highlighted before
    layerid int: the layer to be highlighted now

    Return nothing
    """

    scalars = polydata.GetPointData().GetArray('colors')

    #reset colors of the previous active slice
    if( oldlayerid!=layerid ):
        ind = (layerids==oldlayerid).nonzero()
        for i in range(len(ind[0])):
            scalars.SetValue(ind[0][i],0)

    #set colors of the active slice
    ind = (layerids==layerid).nonzero()
    for i in range(len(ind[0])):
        scalars.SetValue(ind[0][i],1)


def DisplaceLayer(polydata, layerids, layerid, displacement):
    """
    Displaces the highlighted layer by the vector displacement
    
    polydata vtkPolyData: the shape
    layerids 1xN numpy array: stores the slice number for every point
    layerid int: the layer highlighted now
    displacement 3x1 np.array: displacement vector

    Return nothing
    """    
    ind = (layerids==layerid).nonzero()
    
    pts = polydata.GetPoints()
    for i in range(len(ind[0])):
        pt = pts.GetPoint(ind[0][i])
        pt1 = (pt[0]  + displacement[0], pt[1] + displacement[1], pt[2] + displacement[2])
        pts.SetPoint(ind[0][i], pt1)
    
    pts.Modified()
        

class MyInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):
  
    def __init__(self,parent=None):
        self.AddObserver("KeyPressEvent",self.KeyPressEvent)
        self.ActiveSlice = 0 #initial slice marked as active
        self.SliceDisplacement = 1.0 #displacemeent increment in one of the axes
        self.SliceAxis = 0  #slices along x axis
        self.SaveRequested = False #save after the interactor is closed

    def IncrActiveSlice(self, maxid, incr):
        oldslice = self.ActiveSlice
        self.ActiveSlice = self.ActiveSlice + incr
        
        if self.ActiveSlice > maxid:
            self.ActiveSlice = 0
        elif self.ActiveSlice < 0:
            self.ActiveSlice = maxid
        
        return oldslice

            
 
    def KeyPressEvent(self,obj,event):
        key = self.GetInteractor().GetKeySym()
        rnd = self.GetCurrentRenderer()
        actors = rnd.GetActors()
        actors.InitTraversal()

        actors2d = rnd.GetActors2D()
        actors2d.InitTraversal()

        slice_highlight = False
        slice_move = False
        update_text = False
        epi_toggle = False
        endo_toggle = False

        if key=='semicolon':  #highlight slice above
            oldslice = self.IncrActiveSlice(sliceids_endo.max(), 1)
            slice_highlight = True
        elif key=='p':      #highlight slice below
            oldslice = self.IncrActiveSlice(sliceids_endo.max(), -1)
            slice_highlight = True
        elif key=='Up':  #displace slice
            displacement = np.roll( np.array([0, self.SliceDisplacement, 0]), self.SliceAxis )
            slice_move = True
        elif key=='Down': #displace slice
            displacement = np.roll( np.array([0, -self.SliceDisplacement, 0]), self.SliceAxis )
            slice_move = True
        elif key=='Left': #displace slice
            displacement = np.roll( np.array([0, 0, self.SliceDisplacement]), self.SliceAxis )
            slice_move = True
        elif key=='Right': #displace slice
            displacement = np.roll( np.array([0, 0, -self.SliceDisplacement]), self.SliceAxis )
            slice_move = True
        elif key=='bracketleft': #increase displaceemnt step
            self.SliceDisplacement = self.SliceDisplacement*2
            update_text = True
        elif key=='apostrophe': #decrease displaceemnt step
            self.SliceDisplacement = self.SliceDisplacement/2
            update_text = True
        elif key=='1': #toggle apicardium visibility
            epi_toggle = True
        elif key=='2': #toggle endocardium visibility
            endo_toggle = True
        elif key=='s': #toggle SAVE flag
            self.SaveRequested = not self.SaveRequested
            update_text = True


        while epi_toggle:
            actor = actors.GetNextProp()

            if actor==actor_epi:
                actor.SetVisibility( not actor.GetVisibility() )

            if actor==actors.GetLastProp():
                break


        while endo_toggle:
            actor = actors.GetNextProp()

            if actor==actor_endo:
                actor.SetVisibility( not actor.GetVisibility() )

            if actor==actors.GetLastProp():
                break



        
        while update_text:
            actor = actors2d.GetNextActor2D()
        
            if actor==textactor:
                actor.SetInput( CreateBottomText(self.SliceDisplacement, self.SaveRequested) )
                actor.Modified()

            if actor==actors2d.GetLastActor2D():
                break        
        
        while slice_highlight:
            actor = actors.GetNextProp()
            pd = actor.GetMapper().GetInputAsDataSet()

            if actor==actor_endo:
                ChangeLayerColor(pd, sliceids_endo, oldslice, self.ActiveSlice)
            elif actor==actor_epi:
                ChangeLayerColor(pd, sliceids_epi, oldslice, self.ActiveSlice )
                
            pd.Modified()            

            if actor==actors.GetLastProp():
                break


        while slice_move:
            actor = actors.GetNextProp()
            pd = actor.GetMapper().GetInputAsDataSet()

            if actor==actor_endo:
                DisplaceLayer(pd, sliceids_endo, self.ActiveSlice, displacement)
            elif actor==actor_epi:
                DisplaceLayer(pd, sliceids_epi, self.ActiveSlice, displacement)
                
            pd.Modified()            

            if actor==actors.GetLastProp():
                break


        #redraw
        self.GetInteractor().GetRenderWindow().Render()        
        self.OnKeyPress()   
        return
        
        
#------------------------------------------------------------
#
#
#
#
#     MAIN CODE
#
#
#
#------------------------------------------------------------
    
#
#  Parameters, adapt to your taste
#
input_filename = 'points.npz'  #numpy file with epi and endo
output_filename = 'points_corrected.npz'  #numpy file to save the corrections
window_size = [1200,800] #size for the render window

#
#
#  The remianing processing
#


npzfile = np.load(input_filename)
pts_endo = npzfile['pts_endo']
pts_epi = npzfile['pts_epi']
sliceids_endo = npzfile['sliceids_endo']
sliceids_epi = npzfile['sliceids_epi']


#allids = np.concatenate((sliceids_endo,sliceids_epi))
activeslice = 0
sliceaxis = 0 #slices along X axis

#endo points
shape_endo = Points2Shape( pts_endo )
shape_epi = Points2Shape( pts_epi )
ChangeLayerColor(shape_endo, sliceids_endo, 0, 0)
ChangeLayerColor(shape_epi, sliceids_epi, 0, 0)




# Visualize
mapper_endo = vtk.vtkPolyDataMapper()
mapper_endo.SetInput(shape_endo)
actor_endo = vtk.vtkActor()
actor_endo.SetMapper(mapper_endo)
actor_endo.GetProperty().SetPointSize(5)

mapper_epi = vtk.vtkPolyDataMapper()
mapper_epi.SetInput(shape_epi)
actor_epi = vtk.vtkActor()
actor_epi.SetMapper(mapper_epi)
actor_epi.GetProperty().SetPointSize(1)
 
 
 
renderer = vtk.vtkRenderer()
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

myinteractor = MyInteractorStyle()
myinteractor.ActiveSlice = activeslice

textactor = vtk.vtkTextActor()
textactor.GetTextProperty().SetFontSize(14)
textactor.SetPosition2(50,20)
textactor.SetInput( CreateBottomText(myinteractor.SliceDisplacement) )
textactor.GetTextProperty().SetColor((1,0,0))

 
renderer.AddActor(actor_endo)
renderer.AddActor(actor_epi)
renderer.AddActor(textactor)
renderWindow.SetSize(window_size[0],window_size[1])
renderWindow.Render()
 

renderWindowInteractor.SetInteractorStyle( myinteractor )
renderWindowInteractor.Start()

#save the points
if myinteractor.SaveRequested:
    pts_endo = mytools.ExtractVTKPoints(shape_endo)
    pts_epi = mytools.ExtractVTKPoints(shape_epi)
    sliceids_endo = npzfile['sliceids_endo']
    sliceids_epi = npzfile['sliceids_epi']
    np.savez(output_filename, pts_epi=pts_epi, pts_endo=pts_endo, sliceids_endo=sliceids_endo, sliceids_epi=sliceids_epi)