import sys
import vtk
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

        self.mesh = vtk.vtkPolyData()
        self.mesh.ShallowCopy( rd.GetOutput() )

        self.plane = vtk.vtkPlane()
        self.clip = vtk.vtkClipPolyData()
        self.clip.SetInputData( self.mesh )
        self.clip.SetClipFunction( self.plane )
        self.clip.GenerateClippedOutputOff ()
        self.clip.Update()
        self.cut_mesh = self.clip.GetClippedOutput()



        # ImplicitPlaneWidget - outlines and plane
        self.planeWidget = vtk.vtkImplicitPlaneWidget()
        self.planeWidget.SetInteractor( self.iren )
        self.planeWidget.SetPlaceFactor( 1 )
        self.planeWidget.SetInputConnection( self.clip.GetOutputPort() )
        self.planeWidget.PlaceWidget()
        self.planeWidget.AddObserver("InteractionEvent", self.myCallback)
        self.planeWidget.Off() # ImplicitPlaneWidget end


        # Create a mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection( self.clip.GetOutputPort() )

        # Create an actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        self.iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())

        self.ren.AddActor(actor)

        self.ren.ResetCamera()

        self.frame.setLayout(self.vl)
        self.setCentralWidget(self.frame)

        self.ren.SetBackground(0,0,0)

        self.create_buttons()

        self.show()
        self.iren.Initialize()
        self.iren.Start()

    def create_buttons(self):
        self.cut_button = QPushButton('Cut plane', self)
        self.cut_button.setToolTip('Define cutting plane')
        self.cut_button.move(0,0)
        self.cut_button.clicked.connect(self.on_click_cutplane)

        button2 = QPushButton('LV Endo', self)
        button2.setToolTip('Define LV endo')
        button2.move(self.cut_button.frameGeometry().width(),0)
        button2.clicked.connect(self.on_click_lvendo)


    def myCallback(self, obj, event):
        obj.GetPlane(self.plane)
        #self.selectActor.VisibilityOn() 

    @pyqtSlot()
    def on_click_cutplane(self):
        if ( self.planeWidget.GetEnabled() == 1 ):
            self.clip.Update()
            self.planeWidget.SetEnabled( 0 )
            self.cut_button.setText("Cut plane")
        else:
            self.planeWidget.SetEnabled( 1 )
            self.cut_button.setText("Apply Cut")

        print('Cutting plane click')

    @pyqtSlot()
    def on_click_lvendo(self):
        print('LV endo click')



if __name__ == "__main__":
    app = Qt.QApplication(sys.argv)
    window = MainWindow()
    sys.exit(app.exec_())
