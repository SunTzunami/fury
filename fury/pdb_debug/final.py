import numpy as np
from fury import window, actor, ui
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk



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
SecondaryStructuresBegin = []
SecondaryStructuresEnd = []
current_model_number = 1

table = vtk.vtkPeriodicTable()

# parsing the mmCIF file for information about coordinates and atoms
pdbx_file_name = '1pgb.pdb'
pdbxfile = open(pdbx_file_name, 'r')
pdbx_lines = pdbxfile.readlines()

i = 0
for line in pdbx_lines:
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
SecondaryStructuresBegin = np.empty(n, dtype=bool)
SecondaryStructuresBegin[:] = np.NaN

SecondaryStructuresEnd = np.empty(n, dtype=bool)
SecondaryStructuresEnd[:] = np.NaN

Residue = np.array(Residue, dtype='int')


for i in range(n):
    SecondaryStructures[i] = ord('c')
    resi = Residue[i]

    for j in range(len(Sheets)):
        sheet = Sheets[j]
        if Chain[i] != sheet[0] or resi < sheet[1] or resi > sheet[3]:
            continue
        # if resi < sheet[1]:
        #     continue
        # if resi > sheet[3]:
        #     continue
        SecondaryStructures[i] = ord('s')
        if resi == sheet[1]:
            SecondaryStructuresBegin[i] = True
        if resi == sheet[3]:
            SecondaryStructuresEnd[i] = True
    # print(1)
    for j in range(len(Helix)):
        helix = Helix[j]
        if Chain[i] != helix[0] or resi < helix[1] or resi > helix[3]:
            continue
        # if resi < helix[1]:
        #     continue
        # if resi > helix[3]:
        #     continue
        SecondaryStructures[i] = ord('h')
        if resi == helix[1]:
            SecondaryStructuresBegin[i] = True
        if resi == helix[3]:
            SecondaryStructuresEnd[i] = True



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

# to check if the arrays obtained after parsing are having same values as the ones supplied by VTK
reader = vtk.vtkPDBReader()
protein_name = '1pgb.pdb'
reader.SetFileName(protein_name)
reader.Update()
polydata = reader.GetOutput()

# secondary structures begin and secondary structures end don't have the required input
arr = polydata.GetPointData().GetArray(5) # secondary structure begin
output.GetPointData().AddArray(arr)

arr = polydata.GetPointData().GetArray(6) # secondary structure end
output.GetPointData().AddArray(arr)

# # for secondary structures begin
# SecondaryStructuresBegin = np.array(SecondaryStructuresBegin)
# s_sb = numpy_to_vtk(num_array=SecondaryStructuresBegin, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR)
# s_sb.SetName("secondary_structures_begin")
# output.GetPointData().AddArray(s_sb)

# # for secondary structures end
# SecondaryStructuresEnd = np.array(SecondaryStructuresEnd)
# s_se = numpy_to_vtk(num_array=SecondaryStructuresEnd, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR)
# s_se.SetName("secondary structures_end")
# output.GetPointData().AddArray(s_se)

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

# to check if the arrays obtained after parsing are the same as generated by vtlPDBParser
for i in range(9):
    if(i==1): # skip if the atomtypestrings is encountered as numpy operations don't work on VTKStrings
        print(i, True)
        continue
    print(i, np.array_equal(vtk_to_numpy(polydata.GetPointData().GetArray(i)), vtk_to_numpy(output.GetPointData().GetArray(i))))

# print(vtk_to_numpy(polydata.GetPointData().GetArray(4)), vtk_to_numpy(output.GetPointData().GetArray(i))))

# print(output.GetPointData())
# print(vtk_to_numpy(output.GetPointData().GetArray(6)))
# print(vtk_to_numpy(polydata.GetPointData().GetArray(6)))


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

scene.add(actor)
window.show(scene, size=(500, 500))
