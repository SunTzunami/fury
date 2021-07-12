import numpy as np
from fury import window, actor, ui
import vtk
from vtk.util.numpy_support import numpy_to_vtk

NumberOfAtoms = 0

Points = []
AtomType = []
AtomTypeStrings = []
Model = []
Sheets = []
Helix = []
Residue = []
Chain = []
IsHetatm = []
SecondaryStructures = []
current_model_number = 1

table = vtk.vtkPeriodicTable()

# parsing the mmCIF file for information about coordinates and atoms
pdb_code = '1pgb.pdb'
pdbfile = open(pdb_code, 'r')
pdb_lines = pdbfile.readlines()

i = 0
for line in pdb_lines:
        l = line.split()
        try:
            if l[0] == 'ATOM' or l[0] == 'HETATM':
                if l[-1]!='H':
                    i += 1
                    # print(2)
                    coorX, coorY, coorZ = float(l[6]), float(l[7]), float(l[8])
                    resi = l[5]
                    chain = ord(l[4])
                    Points += [[coorX, coorY, coorZ]]
                    Residue += [resi]
                    Chain += [chain]
                    AtomType += [table.GetAtomicNumber(l[-1])]
                    AtomTypeStrings += [l[2]]
                    Model += [current_model_number]
                    NumberOfAtoms += 1
                    if(l[0]=='HETATM'):
                        IsHetatm += [1]
                    else:
                        IsHetatm += [0]
            if l[0] == 'SHEET':
                startChain = ord(l[5])
                startResi = int(l[6])
                endChain = ord(l[8])
                endResi = int(l[9])
                r = [startChain, startResi, endChain, endResi]
                Sheets += [r]
            if l[0] == 'HELIX':
                startChain = ord(l[4])
                startResi = int(l[5])
                endChain = ord(l[7])
                endResi = int(l[8])
                r = [startChain, startResi, endChain, endResi]
                Helix += [r]
        except:
            continue

n = len(Points)
SecondaryStructures = np.ones(n)
Residue = np.array(Residue, dtype='int')

for i in range(n):
    SecondaryStructures[i] = ord('c')
    resi = Residue[i]
    for j in range(len(Sheets)):
        sheet = Sheets[j]
        if Chain[i] != sheet[0] or resi < sheet[1] or resi > sheet[3]:
            continue
        SecondaryStructures[i] = ord('s')

    # print(1)
    for j in range(len(Helix)):
        helix = Helix[j]
        if Chain[i] != helix[0] or resi < helix[1] or resi > helix[3]:
            continue
        SecondaryStructures[i] = ord('h')


output = vtk.vtkPolyData()

# for atom type
AtomType = np.array(AtomType)
atom_type = numpy_to_vtk(num_array=AtomType, deep=True, array_type=vtk.VTK_ID_TYPE)
atom_type.SetName("atom_type")

output.GetPointData().AddArray(atom_type)

# for atom type strings
atom_types = vtk.vtkStringArray()
atom_types.SetName("atom_types")
atom_types.SetNumberOfTuples(len(AtomTypeStrings))
for i in range(len(AtomTypeStrings)):
    atom_types.SetValue(i, AtomTypeStrings[i])

AtomTypeStrings = np.array(AtomTypeStrings)

output.GetPointData().AddArray(atom_types)

# for residue
residue = numpy_to_vtk(num_array=Residue, deep=True, array_type=vtk.VTK_ID_TYPE)
residue.SetName("residue")
output.GetPointData().AddArray(residue)

# for chain
Chain = np.array(Chain)
chain = numpy_to_vtk(num_array=Chain, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR)
chain.SetName("chain")
output.GetPointData().AddArray(chain)

# for secondary structures
s_s = numpy_to_vtk(num_array=SecondaryStructures, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR)
s_s.SetName("secondary_structures")
output.GetPointData().AddArray(s_s)


# for secondary structures begin
newarr = np.ones(n)
s_sb = numpy_to_vtk(num_array=newarr, array_type=vtk.VTK_UNSIGNED_CHAR)
s_sb.SetName("secondary_structures_begin")
output.GetPointData().AddArray(s_sb)

# for secondary structures end
newarr = np.ones(n)
s_se = numpy_to_vtk(num_array=newarr, array_type=vtk.VTK_UNSIGNED_CHAR)
s_se.SetName("secondary_structures_end")
output.GetPointData().AddArray(s_se)


# for ishetatm
IsHetatm = np.array(IsHetatm)
ishetatm = numpy_to_vtk(num_array=IsHetatm, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR)
ishetatm.SetName("ishetatm")
output.GetPointData().AddArray(ishetatm)

# for model
Model = np.array(Model)
model = numpy_to_vtk(num_array=Model, deep=True, array_type=vtk.VTK_UNSIGNED_INT)
model.SetName("model")
output.GetPointData().AddArray(model)


from fury.utils import numpy_to_vtk_points
Points = np.array(Points)
points = numpy_to_vtk_points(Points)
output.SetPoints(points)

ribbonFilter = vtk.vtkProteinRibbonFilter()
ribbonFilter.SetInputData(output)
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(ribbonFilter.GetOutputPort())
actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetDiffuse(1)
actor.GetProperty().SetSpecular(0.9)
actor.GetProperty().SetSpecularPower(100.0)


scene = window.Scene()
scene.background((1, 1, 1))

###############################################################################
# Dimensions of the output screen
screen_x_dim = 600
screen_y_dim = 600
dims = (screen_x_dim, screen_y_dim)

###############################################################################
# creating a ShowManager object
showm = window.ShowManager(scene, size=dims, reset_camera=True,
                           order_transparent=True)

tb = ui.TextBlock2D(text=pdb_code, position=(screen_x_dim/2-60,
                    screen_y_dim/12), font_size=30, color=(0, 0, 0))
scene.add(tb)
scene.add(actor)

interactive = True
if interactive:
    window.show(scene, size=dims, title=pdb_code)
window.record(scene, size=dims, out_path=pdb_code+'.png')
