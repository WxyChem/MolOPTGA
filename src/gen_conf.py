import os
from rdkit import Chem
from rdkit.Chem.AllChem import EmbedMolecule, MMFFOptimizeMolecule

def genconf(mol, path, seed: int = 42, minimize: bool = False):
    # Add H for the molecule
    mol = Chem.AddHs(mol)
    
    # Generate 3D conformation
    EmbedMolecule(mol, randomSeed=seed)
    
    # Using MMFF forcefield to optimize the conformation
    if minimize:
        MMFFOptimizeMolecule(mol)
    
    # Output the SDF file
    writer = Chem.SDWriter(path)
    writer.write(mol)
    writer.close()


def multi_genconf(mols, folder):
    paths = []
    
    # Create the save folder for the SDF files
    try:
        os.mkdir(folder)
    except FileExistsError:
        pass
    
    for i, mol in enumerate(mols):
        if mol:
            path = f"{folder}/{i}.sdf"
            genconf(mol, path)
            paths.append(path)
        else:
            continue

    return paths
