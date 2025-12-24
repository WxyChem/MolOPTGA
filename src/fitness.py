import numpy as np
from rdkit import Chem
from .sascores import calculateScore
from .docking import MolecularDocking
from .gen_conf import multi_genconf
from .gen_pdbqt import multi_genpdbqt

from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.Lipinski import NumHDonors, NumHAcceptors, NumRotatableBonds
from rdkit.Chem.Descriptors import MolWt

class FITNESS:
    def __init__(self, program_path, receptor_path):
        self.ring_size = 6
        self.md = MolecularDocking(program_path, receptor_path)
    
    # calculate the number of large ring
    def rings_count(self, mol):
        ring = mol.GetRingInfo()
        ams = ring.AtomRings()
        large_ams = [am for am in ams if len(am) > self.ring_size]
        return len(large_ams)
    
    # calculate the plogp and sascore for each molecule
    def run(self, molecules: list):
        mols = []
        plogps = []
        sas = []
        for mol in molecules:
            mw = MolWt(mol)
            logp = MolLogP(mol)
            hbas = NumHAcceptors(mol)
            hbds = NumHDonors(mol)
            
            sa = calculateScore(mol)
            rings = self.rings_count(mol)
            
            if mw <= 500 and logp <= 5.0 and hbas <= 10 and hbds <= 5:
                plogp = logp - sa - rings
                plogps.append(plogp)
                sas.append(sa)
                mols.append(Chem.MolToSmiles(mol))
            else:
                continue
        
        return mols, plogps, sas
    
    # calculate the vina score for each molecule
    def run_d(self, mols, box_dims, folder_conf, folder_pdbqt, folder_output):
        paths = multi_genconf(mols, folder_conf)
        multi_genpdbqt(paths, folder_pdbqt)
        
        vinas, paths = self.md.run(folder_pdbqt, box_dims, folder_output)
        
        return vinas, paths


class Pareto_Front:
    def __init__(self, objectives: tuple, maximize: tuple, tolerance: float = 1e-12):
        self.objectives = objectives
        self.maximize = maximize
        self.tolerance = tolerance

    def pareto(self, data):
        """
            Compute the Pareto Front（non domination set)
        """
        assert all(col in data.columns for col in self.objectives), "Some objectives not in the DataFrame."
        
        # Build the target martix
        X = data.loc[:, self.objectives].astype(float).to_numpy(copy=True)

        # Convert maximize task to minimize task
        for i, col in enumerate(self.objectives):
            if col in self.maximize:
                X[:, i] = -X[:, i]

        num = X.shape[0]

        # Create the dominated minimize task martix
        dominated = np.zeros(num, dtype=bool)

        for i in range(num):
            if dominated[i]:
                continue

            xi = X[i]

            for j in range(num):
                if i == j or dominated[j]:
                    continue

                xj = X[j]
                # The condition of j dominate i
                le = np.all(xj <= xi + self.tolerance)
                lt = np.any(xj < xi - self.tolerance)
                if le and lt:
                    dominated[i] = True
                    break

        pareto_mask = ~dominated
        return pareto_mask, data.loc[pareto_mask].copy()

def mol_check(mol):
    mc = False
    if mol:
        mw = MolWt(mol)
        logp = MolLogP(mol)
        hbas = NumHAcceptors(mol)
        hbds = NumHDonors(mol)
        if mw <= 500 and logp <= 5.0 and hbas <= 10 and hbds <= 5:
            mc = True
        else:
            pass
    else:
        pass
    
    return mc
    