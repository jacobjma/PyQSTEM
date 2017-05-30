import numpy as np

def frozen_phonon(atoms, MSD, num_structures=1):

    if (not isinstance(MSD,dict)) and (len(atoms) != len(MSD)):
        raise RuntimeError('Length of ´MSD´ have to match number of atoms')
    elif isinstance(MSD,dict):
        MSD_dict = MSD.copy()
        MSD = []
        for i in range(len(atoms)):
            MSD.append(MSD_dict[atoms[i].symbol])

    atoms_list = []
    for i in range(num_structures):
        atoms_copy=atoms.copy()
        r = np.random.randn(len(atoms))*np.sqrt(MSD)
        theta = np.random.rand(len(atoms))*np.pi
        phi = np.random.rand(len(atoms))*2*np.pi

        displacement = np.array([r*np.sin(theta)*np.cos(phi),
                                r*np.sin(theta)*np.sin(phi),
                                r*np.cos(theta)]).T

        positions = atoms.get_positions() + displacement

        atoms_copy.set_positions(positions)

        atoms_list.append(atoms_copy)

    if num_structures == 1:
        return atoms_list[0]
    else:
        return atoms_list
