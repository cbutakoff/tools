import sys
import vtk
import numpy as np
from PyQt5 import QtCore, QtGui
from PyQt5 import Qt
from PyQt5.QtWidgets import  QPushButton
from PyQt5.QtCore import pyqtSlot

from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

class MainWindow(Qt.QMainWindow):

    def __init__(self, parent = None):
        Qt.QMainWindow.__init__(self, parent)

        self.setGeometry(0, 0, 800, 800)   

        self.frame = Qt.QFrame()
        self.vl = Qt.QVBoxLayout()
        self.vtkWidget = QVTKRenderWindowInteractor(self.frame)
        self.vl.addWidget(self.vtkWidget)

        self.ren = vtk.vtkRenderer()
        self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()

        mesh_filename = Qt.QApplication.instance().arguments()[1]
        rd = vtk.vtkDataSetReader()
        rd.SetFileName( mesh_filename )
        rd.Update()

        self.full_mesh = vtk.vtkPolyData()
        self.full_mesh.ShallowCopy( rd.GetOutput() )
        self.full_mesh.ComputeBounds()
        bounds = np.array(self.full_mesh.GetBounds()) #(xmin,xmax, ymin,ymax, zmin,zmax).

        self.plane = vtk.vtkPlane()
        self.plane.SetOrigin( (bounds[::2] + bounds[1::2])/2 )
        self.clip = vtk.vtkClipPolyData()
        self.clip.SetInputData( self.full_mesh )
        self.clip.SetClipFunction( self.plane )
        self.clip.GenerateClippedOutputOff ()
        self.clip.Update()



        # ImplicitPlaneWidget - outlines and plane
        self.planeWidget = vtk.vtkImplicitPlaneWidget()
        self.planeWidget.SetInteractor( self.iren )
        self.planeWidget.SetPlaceFactor( 1 )
        self.planeWidget.SetInputData( self.full_mesh )
        self.planeWidget.PlaceWidget()
        self.planeWidget.SetOrigin(self.plane.GetOrigin()) 
        self.planeWidget.OutlineTranslationOff ()
        self.planeWidget.AddObserver("InteractionEvent", self.myCallback)
        self.planeWidget.Off() # ImplicitPlaneWidget end
        self.planeWidget.GetPlane(self.plane)


        # Create a mapper
        self.main_mesh_mapper = vtk.vtkPolyDataMapper()
        self.main_mesh_mapper.SetInputData( self.full_mesh )

        # Create an actor
        actor = vtk.vtkActor()
        actor.SetMapper(self.main_mesh_mapper)

        self.iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())

        self.ren.AddActor(actor)

        self.ren.ResetCamera()

        self.frame.setLayout(self.vl)
        self.setCentralWidget(self.frame)

        self.ren.SetBackground(0,0,0)

        self.create_buttons_main()
        self.create_cut_toolbar()
        self.visibility_toolbar(self.cut_toolbar, False)

        self.show()
        self.iren.Initialize()
        self.iren.Start()

    def create_buttons_main(self):
        left = 0
        self.main_toolbar = {}
        self.main_toolbar['cut_button'] = QPushButton('Cut plane', self)
        self.main_toolbar['cut_button'].setToolTip('Define cutting plane')
        self.main_toolbar['cut_button'].move(left,0)
        self.main_toolbar['cut_button'].clicked.connect(self.on_click_cutplane)

        left += self.main_toolbar['cut_button'].frameGeometry().width()
        self.main_toolbar['lv_endo_seed'] = QPushButton('LV Endo', self)
        self.main_toolbar['lv_endo_seed'].setToolTip('Define LV endo')
        self.main_toolbar['lv_endo_seed'].move(left,0)
        self.main_toolbar['lv_endo_seed'].clicked.connect(self.on_click_lvendo)


    def create_cut_toolbar(self):
        left = 0
        self.cut_toolbar = {}
        self.cut_toolbar['flip'] = QPushButton('Flip normal', self)
        self.cut_toolbar['flip'].setToolTip('Flip plane normal')
        self.cut_toolbar['flip'].move(left,self.main_toolbar['cut_button'].frameGeometry().height())
        self.cut_toolbar['flip'].clicked.connect(self.on_click_flipcutnormal)

        left += self.cut_toolbar['flip'].frameGeometry().width()
        self.cut_toolbar['apply'] = QPushButton('Apply', self)
        self.cut_toolbar['apply'].setToolTip('Apply')
        self.cut_toolbar['apply'].move(left, self.main_toolbar['cut_button'].frameGeometry().height())
        self.cut_toolbar['apply'].clicked.connect(self.on_click_cut_apply)

        left += self.cut_toolbar['apply'].frameGeometry().width()
        self.cut_toolbar['cancel'] = QPushButton('Cancel', self)
        self.cut_toolbar['cancel'].setToolTip('Cancel')
        self.cut_toolbar['cancel'].move(left, self.main_toolbar['cut_button'].frameGeometry().height())
        self.cut_toolbar['cancel'].clicked.connect(self.on_click_cut_cancel)


    def visibility_toolbar(self, toolbar, visible):
        for key,value in toolbar.items():
            value.setVisible(visible)

    def enable_toolbar(self, toolbar, enabled):
        for key,value in toolbar.items():
            value.setEnabled(enabled)


    def myCallback(self, obj, event):
        obj.GetPlane(self.plane)
        #self.selectActor.VisibilityOn() 

    @pyqtSlot()
    def on_click_cutplane(self):
        self.planeWidget.SetEnabled( 1 - self.planeWidget.GetEnabled() )

        self.main_mesh_mapper.SetInputData( self.full_mesh )
        self.visibility_toolbar(self.cut_toolbar, True)
        self.enable_toolbar(self.main_toolbar, False)

        self.update_view()
        

    @pyqtSlot()
    def on_click_lvendo(self):
        print('LV endo click')


    @pyqtSlot()
    def on_click_flipcutnormal(self):
        self.planeWidget.SetNormal( -np.array(self.planeWidget.GetNormal()) ) 
        self.planeWidget.GetPlane(self.plane)
        self.update_view()

    @pyqtSlot()
    def on_click_cut_apply(self):
        self.planeWidget.SetEnabled( 1 - self.planeWidget.GetEnabled() )
        self.clip.Update()
        self.main_mesh_mapper.SetInputData( self.clip.GetOutput() )
        self.visibility_toolbar(self.cut_toolbar, False)
        self.enable_toolbar(self.main_toolbar, True)
        self.update_view()


    @pyqtSlot()
    def on_click_cut_cancel(self):
        self.planeWidget.SetEnabled( 1 - self.planeWidget.GetEnabled() )
        self.main_mesh_mapper.SetInputData( self.full_mesh )
        self.visibility_toolbar(self.cut_toolbar, False)
        self.enable_toolbar(self.main_toolbar, True)
        self.update_view()



    def update_view(self):
        self.vtkWidget.GetRenderWindow().Render()


if __name__ == "__main__":
    app = Qt.QApplication(sys.argv)
    window = MainWindow()
    sys.exit(app.exec_())
