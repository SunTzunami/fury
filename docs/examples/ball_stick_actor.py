import numpy as np
from fury import window, actor, disable_warnings
disable_warnings()

# visualize methane (CH4) via ball and stick diagram


# ball and stick actor
def bs_actor(atom_coordinates, elem_names, bond_coords):
    # bond_elems finds the element which are bonded by comparing bond coordinates and atom's coordinates
    bonds_elems = []
    for bond_coord in bond_coords:
        ename1 = elem_names[np.where(np.all(atom_coords==bond_coord[0],axis=1))]
        ename2 = elem_names[np.where(np.all(atom_coords==bond_coord[1],axis=1))]
        bonds_elems.append([ename1, ename2])
    bonds_elems = np.array(bonds_elems)

    # reshaping bond coordinates and bond_elems
    bonds_elems = bonds_elems.reshape(len(bond_coords)*2)
    bond_coords = bond_coords.reshape(len(bond_coords)*2, 3)

    # cpk coloring scheme
    cpkr = {'H': [1, 1, 1, 1.2], 'C': [144/255, 144/255, 144/255, 1.7], 'N': [48/255, 80/255, 248/255, 1.55], 'NA': [171/255, 92/255, 242/255, 2.27],
                'O': [255/255, 13/255, 13/255, 1.52], 'P': [255/255, 128/255, 0, 1.8],
                'S': [255/255, 255/255, 48/255, 1.8], 'K':[143/255, 64/255, 212/255, 2.75], 'MG': [138/255, 1, 0, 1.73], 'CL': [31/255, 240/255, 31/255, 1.75],
                'CA': [61/255, 255/255, 0, 2.31], 'ZN':[125/255, 128/255, 176/255, 1.39], 'FE': [224/255, 102/255, 51/255, 1.16],
                'CD': [255/255, 217/255, 143/255, 1.58], 'CO': [240/255, 144/255, 160/255, 2.4], 'NI': [80/255, 208/255, 80/255, 1.63]
            # No vanderwaals radius data for Iron, using covalent single bond radius data instead
            }

    # array having names of the unique atoms
    unique_elem_types = np.unique(elem_names)

    # assigning the appropriate colors to bonds and atoms
    colors_bonds = np.ones((len(bonds_elems), 3))
    colors_atoms = np.ones((len(atom_coords), 3))
    for i, typ in enumerate(unique_elem_types):
        colors_atoms[elem_names == typ] = cpkr[typ][:3]
        colors_bonds[bonds_elems == typ] = cpkr[typ][:3]

    # splitting the bond into two to color based on which atoms it connects
    bonds_finals = []
    for i in range(0, len(bond_coords), 2):
        mid = (bond_coords[i]+bond_coords[i+1])/2
        bonds_finals.append([bond_coords[i], mid])
        bonds_finals.append([mid, bond_coords[i+1]])
    bonds_finals = np.array(bonds_finals, dtype=float)

    balls = actor.sphere(atom_coords, colors_atoms, radii=0.25)
    sticks = actor.streamtube(bonds_finals, colors=colors_bonds)

    return balls, sticks



# coordinates of atoms
atom_coords = np.array([[3.875,   0.678,  -8.417],
                [3.800,   1.690,  -8.076],
                [4.907,   0.410,  -8.516],
                [3.406,   0.026,  -7.711],
                [3.389,   0.583,  -9.366]], dtype=float)

# array storing names of elements of the atoms by order
elem_names = np.array(['C', 'H', 'H', 'H', 'H'])

# array storing coordinates of the bonds
bond_coords = np.array([[[3.875,   0.678,  -8.417],
                   [3.800,   1.690,  -8.076]],
                   [[3.875,   0.678,  -8.417],
                   [4.907,   0.410,  -8.516]],
                   [[3.875,   0.678,  -8.417],
                   [3.406,   0.026,  -7.711]],
                   [[3.875,   0.678,  -8.417],
                   [3.389,   0.583,  -9.366]]], dtype=float)

ball, sticks = bs_actor(atom_coords, elem_names, bond_coords)

scene = window.Scene()
scene.add(ball, sticks)
scene.background((1, 1, 1))
scene.set_camera(position=(0, 0, 100), focal_point=(0, 0, 0),
                        view_up=(0, 1, 0))

scene.add()
scene.add()
window.show(scene, size=(500, 500))