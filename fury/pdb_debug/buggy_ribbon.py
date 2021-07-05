import numpy as np
from fury import window, actor, ui
from fury.shaders import load
import vtk

# from vtk.util.numpy_support import vtk_to_numpy

NumberOfAtoms = 0
Points = vtk.vtkPoints()
# Points.Allocate(500)

AtomType = vtk.vtkUnsignedIntArray()
AtomType.SetNumberOfComponents(1)
# AtomType.Allocate(500)

AtomTypeStrings = vtk.vtkStringArray()
# AtomTypeStrings.Allocate(500)

Model = vtk.vtkUnsignedIntArray()
# Model.Allocate(500)

Sheets = vtk.vtkIntArray()
# Sheets.Allocate(500)
Sheets.SetNumberOfComponents(4)

Helix = vtk.vtkIntArray()
# Helix.Allocate(500)
Helix.SetNumberOfComponents(4)

Residue = vtk.vtkUnsignedIntArray()
# Residue.Allocate(500)

Chain = vtk.vtkUnsignedCharArray()
# Chain.Allocate(500)

IsHetatm = vtk.vtkUnsignedCharArray()
# IsHetatm.Allocate(500)

SecondaryStructures = vtk.vtkUnsignedCharArray()
# SecondaryStructures.Allocate(500)

SecondaryStructuresBegin = vtk.vtkUnsignedCharArray()
# SecondaryStructuresBegin.Allocate(500)

SecondaryStructuresEnd = vtk.vtkUnsignedCharArray()
# SecondaryStructuresEnd.Allocate(500)

current_model_number = 1

table = vtk.vtkPeriodicTable()

# parsing the mmCIF file for information about coordinates and atoms
pdbx_file_name = '1pgb.pdb'
pdbxfile = open(pdbx_file_name, 'r')
pdbx_lines = pdbxfile.readlines()

atomic_num = []

for line in pdbx_lines:
        l = line.split()
        try:
            if l[0] == 'ATOM' or l[0] == 'HETATM':
                if l[-1]!='H':
                    # print(2)
                    coorX, coorY, coorZ = float(l[6]), float(l[7]), float(l[8])
                    resi = l[1]
                    chain = l[4]
                    Points.InsertNextPoint(coorX, coorY, coorZ)
                    Residue.InsertNextValue(resi)
                    Chain.InsertNextValue(chain)
                    # for later AtomType
                    # print(table.GetAtomicNumber(l[-1]))
                    AtomType.InsertNextValue(table.GetAtomicNumber(l[-1]))
                    atomic_num += table.GetAtomicNumber(table.GetAtomicNumber(l[-1]))
                    AtomTypeStrings.InsertNextvalue(l[2])
                    Model.InsertNextValue(current_model_number)
                    NumberOfAtoms += 1
                    if(l[0]=='HETATM'):
                        IsHetatm.InsertNexValue(True)
            if l[0] == 'SHEET':
                startChain = l[5]
                startResi = l[6]
                endChain = l[7]
                endResi = l[9]
                r = [startChain, startResi, endChain, endResi]
                Sheets.InsertNextTypedTuple(r)
            if l[0] == 'HELIX':
                startChain = l[4]
                startResi = l[5]
                endChain = l[7]
                endResi = l[8]
                r = [startChain, startResi, endChain, endResi]
                Helix.InsertNextTypedTuple(r)
        except:
            continue

# Points.Squeeze()
# AtomType.Squeeze()
# AtomTypeStrings.Squeeze()
# Residue.Squeeze()
# IsHetatm.Squeeze()
# Model.Squeeze()

n = Points.GetNumberOfPoints()
SecondaryStructures.SetNumberOfValues(n)
Residue.SetNumberOfValues(n)

# ids = vtk.vtkIdList()
# ids.SetNumberOfIds(1)
for i in range(Points.GetNumberOfPoints()):
    # ids.SetId(0, i)
    SecondaryStructures.SetValue(i, 99)
    resi = Residue.GetValue(i)
    for j in range(Sheets.GetNumberOfTuples()):
        sheet = []
        Sheets.GetTypedTuple(j, sheet)
        if Chain.GetValue(i) != sheet[0]:
            continue
        if resi < sheet[1]:
            continue
        if resi > sheet[3]:
            continue
        SecondaryStructures.SetValue(i, 's')
        if resi == sheet[1]:
            SecondaryStructuresBegin.SetValue(i, True)
        if resi == sheet[3]:
            SecondaryStructuresEnd.SetValue(i, True)

    for j in range(Helix.GetNumberOfTuples()):
        helix = []
        Helix.GetTypedTuple(j, helix)
        if Helix.GetValue(i) != helix[0]:
            continue
        if resi < helix[1]:
            continue
        if resi > helix[3]:
            continue
        SecondaryStructures.SetValue(i, 's')
        if resi == helix[1]:
            SecondaryStructuresBegin.SetValue(i, True)
        if resi == helix[3]:
            SecondaryStructuresEnd.SetValue(i, True)

output = vtk.vtkPolyData()

AtomType.SetName("atom_type")
output.GetPointData().AddArray(AtomType)

AtomTypeStrings.SetName("atom_types")
output.GetPointData().AddArray(AtomTypeStrings)

Residue.SetName("residue")
output.GetPointData().AddArray(Residue)

Chain.SetName("chain")
output.GetPointData().AddArray(Chain)

IsHetatm.SetName("ishetatm")
output.GetPointData().AddArray(IsHetatm)

SecondaryStructures.SetName("secondary_structures")
output.GetPointData().AddArray(SecondaryStructures)

SecondaryStructuresBegin.SetName("secondary_structures_begin")
output.GetPointData().AddArray(SecondaryStructuresBegin)

SecondaryStructuresEnd.SetName("secondary_structures_end")
output.GetPointData().AddArray(SecondaryStructuresEnd)

rgb = vtk.vtkUnsignedCharArray()
rgb.SetNumberOfComponents(3)
rgb.Allocate(3 * NumberOfAtoms)
rgb.SetName("rgb_colors")

for i in range(NumberOfAtoms):
    rgb.InsertNextTuple(table.GetDefaultRGBTuple(atomic_num[i]))

output.GetPointData().SetScalars(rgb)

Radii = vtk.vtkFloatArray()
Radii.SetNumberOfComponents(3)
Radii.Allocate(3 * NumberOfAtoms)
Radii.SetName("radius")

for i in range(NumberOfAtoms):
    Radii.InsertNextTuple3(table.GetVDWRadius(atomic_num[i]),
                           table.GetVDWRadius(atomic_num[i]),
                           table.GetVDWRadius(atomic_num[i]))

output.GetPointData().SetVectors(Radii)
#print(dir(output))
# arr = output.GetPointData().GetAbstractArray("secondary_structures")

# print(len(vtk_to_numpy(arr)))
# pdata = vtk.vtkCleanPolyData()
# pdata.SetInputData(output)


ribbonFilter = vtk.vtkProteinRibbonFilter()
ribbonFilter.SetInputData(output)

# pointdata = output.GetPointData()
# atomTypes = pointdata.GetAbstractArray("atom_types")
# for i in range(100):
#     print(AtomType.GetValue(i))
# print(atomTypes)

# for i in range(output.GetNumberOfPoints()):
#     Type = vtk.vtkStdString()
#     Type = atomTypes.GetValue(i)

mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(ribbonFilter.GetOutputPort())

actor = vtk.vtkActor()
actor.SetMapper(mapper)

scene = window.Scene()

# # #tb = ui.TextBlock2D(text=protein_name[:4], position=(220, 400), font_size=20, color=(1, 1, 1))
# # #scene.add(tb)
# scene.background((1, 1, 1))
scene.add(actor)
window.show(scene, size=(500, 500))
