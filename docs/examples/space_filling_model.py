import numpy as np
from fury import window, actor, ui


# list to store names of the atoms
elem_sym_list = []

# list to store coordinates of the atoms
atom_coords = []

# parsing the mmCIF file for information about coordinates and atoms
pdbx_file_name = '4hhb.cif'
pdbxfile = open(pdbx_file_name, 'r')
pdbx_lines = pdbxfile.readlines()
for line in pdbx_lines:
        l = line.split()
        if l[0] == 'ATOM' or l[0] == 'HETATM':
            # obtaining coordinates of atoms
            coorX, coorY, coorZ = float(l[10]), float(l[11]), float(l[12])
            atom_coords.append([coorX, coorY, coorZ])

            # treating ions (eg: N+, O-) like atoms
            l[2] = ''.join([i for i in l[2] if i.isalpha()])
            # obtaining symbol of the atom (eg: C, N, O etc.)
            elem_sym_list.append(l[2])

elem_sym_list = np.array(elem_sym_list)
atom_coords = np.array(atom_coords)


# function that returns the space filling model
# parameters
    # atom_coords: array of floats having atomic coordinates (in angstroms) of shape (N, 3)
    # elem_sym_list: array of strings having information about the element of the atom
    #                example of elem_sym_list: ['C', 'N', 'C', 'O', 'C']

def space_filling_model(atom_coords, elem_sym_list):
    no_atoms = len(atom_coords)
    colors = np.ones((no_atoms, 3))
    radii = np.ones((no_atoms, 1))
    unique_elem_types = np.unique(elem_sym_list)

    # cpk coloring scheme, kindly note that the fourth entry in each element's list corresponds to it's atomic radii (in angstroms)
    cpkr = {'H': [1, 1, 1, 1.2], 'C': [144/255, 144/255, 144/255, 1.7], 'N': [48/255, 80/255, 248/255, 1.55], 'NA': [171/255, 92/255, 242/255, 2.27],
            'O': [255/255, 13/255, 13/255, 1.52], 'P': [255/255, 128/255, 0, 1.8],
            'S': [255/255, 255/255, 48/255, 1.8], 'K':[143/255, 64/255, 212/255, 2.75], 'MG': [138/255, 1, 0, 1.73], 'CL': [31/255, 240/255, 31/255, 1.75],
            'CA': [61/255, 255/255, 0, 2.31], 'ZN':[125/255, 128/255, 176/255, 1.39], 'FE': [224/255, 102/255, 51/255, 1.16],
            'CD': [255/255, 217/255, 143/255, 1.58], 'CO': [240/255, 144/255, 160/255, 2.4], 'NI': [80/255, 208/255, 80/255, 1.63]
           # No Van der Waals radius data for Iron, using covalent single bond radius data instead
          }

    elements = []
    # assigning the appropriate colors and radii to the atoms
    for i, typ in enumerate(unique_elem_types):
        colors[elem_sym_list == typ] = cpkr[typ][:3]
        radii[elem_sym_list == typ] = cpkr[typ][-1]

        # name and color of elements
        elements.append([typ, cpkr[typ][:3]])
    sf_model = actor.sphere(atom_coords, colors=colors, radii=radii)
    return sf_model, elements

#print(len(atom_coords))    # to check no. of atoms being rendered

axes_actor = actor.axes()
scene = window.Scene()
scene.add(axes_actor)
sf_model, elements = space_filling_model(atom_coords, elem_sym_list)


# Dimensions of the output screen
screen_x_dim = 600
screen_y_dim = 600
# ini variable will help in creating the Key of elements and their colors
ini = 30

# Arranging the elements in alphabetic order
elements = np.array(elements)
elements = elements[elements[:, 0].argsort()[::-1]]

showm = window.ShowManager(scene,
                           size=(screen_x_dim, screen_y_dim), reset_camera=True,
                           order_transparent=True)

scene.add(sf_model)

# Key of elements (with the respective colors)
for i, element in enumerate(elements):
    tb = ui.TextBlock2D(text=' - ' + element[0], position=(screen_x_dim-60, ini), font_size=20, color=(1, 1, 1))
    scene.add(tb)
    tb = ui.TextBlock2D(text='  ', position=(screen_x_dim-70, ini+3), font_size=20, bg_color=element[1], color=(0, 0, 0))
    scene.add(tb)
    ini+=30
tb = ui.TextBlock2D(text='Elements', position=(screen_x_dim-100, ini), font_size=20, color=(1, 1, 1))
scene.add(tb)


scene.set_camera(position=(0, 0, 100), focal_point=(0, 0, 0),
                        view_up=(0, 1, 0))
protein_name = pdbx_file_name[:-4]
window.show(scene, title=protein_name, size=(screen_x_dim, screen_y_dim))
#window.snapshot(scene, fname=protein_name+'.png', size=(500, 500))
