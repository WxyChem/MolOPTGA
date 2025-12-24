from rdkit import Chem

def get_attachment_points(mol):
    """Search the molecular generation point"""
    attachment_points = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            if atom.GetDegree() == 1:
                attachment_points.append(atom.GetIdx())

    return attachment_points

def connect_fragments(fragA, idxA, fragB, idxB, bond_type=Chem.rdchem.BondType.SINGLE):
    """Connection process"""
    combined = Chem.CombineMols(fragA, fragB)
    rw_combined = Chem.RWMol(combined)
    offset = fragA.GetNumAtoms()

    # The indices of connection point in fragA is not change.
    idxA_combined = idxA
    # The indices of connection point in fragB would be moved
    idxB_combined = idxB + offset

    # Get the neighbor atoms (real atoms) of connection atom （dummy）
    atomA = fragA.GetAtomWithIdx(idxA)

    if atomA.GetDegree() != 1:
        raise ValueError(f"The connection atom of fragment A {idxA} should have only one neighbor")

    neighborA = atomA.GetNeighbors()[0]
    # The indices of neighbor atoms in fragA
    idx_neighborA = neighborA.GetIdx()
    # Not change after combination
    idx_neighborA_combined = idx_neighborA

    atomB = fragB.GetAtomWithIdx(idxB)
    if atomB.GetDegree() != 1:
        raise ValueError(f"The connection atom of fragment A {idxB} should have only one neighbor")

    neighborB = atomB.GetNeighbors()[0]
    # The indices of neighbor atoms in fragB
    idx_neighborB = neighborB.GetIdx()
    # Moved after combination
    idx_neighborB_combined = idx_neighborB + offset

    # Add new bond for neighbor atoms
    rw_combined.AddBond(idx_neighborA_combined, idx_neighborB_combined, bond_type)

    # Remove dummy atoms
    indices_to_remove = sorted([idxA_combined, idxB_combined], reverse=True)
    for idx in indices_to_remove:
        rw_combined.RemoveAtom(idx)

    # Standardize atoms (fix the bond and H)
    new_mol = rw_combined.GetMol()
    Chem.SanitizeMol(new_mol)
    return new_mol


def connections(fragA, fragB):
    """
    Generate two fragments to connect.
    """
    if fragA is None or fragB is None:
        raise ValueError("Invalid SMILES string!")

    # Get the connection position
    connA = get_attachment_points(fragA)
    connB = get_attachment_points(fragB)
    if not connA or not connB:
        raise ValueError("Can not find the connection point!")

    # connection
    molecules = []
    for i in connA:
        for j in connB:
            try:
                new_mol = connect_fragments(fragA, i, fragB, j)
                molecules.append(new_mol)
            except Exception as e:
                print(f"Connection failed: position_A={i}, position_B={j}, fails: {e}")

    return molecules
