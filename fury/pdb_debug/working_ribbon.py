import vtk
from fury import window, ui, utils


reader = vtk.vtkPDBReader()
protein_name = '1pgb.pdb'
reader.SetFileName(protein_name)

# pdata = vtk.PolyData()
# ribbonFilter.SetInputData(pdata)


ribbonFilter = vtk.vtkProteinRibbonFilter()
ribbonFilter.SetInputConnection(reader.GetOutputPort())
#print(vars(reader))
#var_ = reader.GetOutput()
#print(dir(var_))
#print(var_.GetNumberOfStrips())
#mol = bonder.GetOutput()
#print(type(mol)
#print(reader.GetOutputPort().PrintSelf())

#print(reader.GetOutputPort)

cw = ribbonFilter.GetCoilWidth()
hw = ribbonFilter.GetHelixWidth()
#ribbonFilter.SetDrawSmallMoleculesAsSpheres(False)
sf = ribbonFilter.GetSubdivideFactor()

## improve boundary area of small molecules
ribbonFilter.SetSphereResolution(40)
#print(cw, hw)

mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(ribbonFilter.GetOutputPort())

actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetDiffuse(2)
actor.GetProperty().SetSpecular(0.5)
actor.GetProperty().SetSpecularPower(80.0)
# actor.GetProperty().SetRepresentationToWireframe()
# actor.GetProperty().SetRepresentationToSurface()
# actor.GetProperty().SetRepresentationToPoints()

scene = window.Scene()

#tb = ui.TextBlock2D(text=protein_name[:4], position=(220, 400), font_size=20, color=(1, 1, 1))
#scene.add(tb)
#scene.background((1, 1, 1))
scene.add(actor)
window.show(scene, size=(500, 500))
