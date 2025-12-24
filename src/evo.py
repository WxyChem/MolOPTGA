import random
import numpy as np
import pandas as pd
from tqdm import tqdm
from itertools import product
from rdkit import Chem
from .connection import connections
from .utils import *
from .fitness import FITNESS, Pareto_Front, mol_check

class MolecularOPTGA:
    def __init__(self, 
                 molecule: str, 
                 box_center_dim: list, 
                 program_path: str, 
                 receptor_path: str, 
                 convergence: float = 0.1, 
                 tolerance: int = 1,
                 threshold: dict = {'vina': -9.0, 'plogp': 0.0, 'sas': 5}, 
                 min_parents_size: int = 30, 
                 max_children_size: int = 300, 
                 iters: int = 10, 
                 output: str = 'Generation.csv'
        ):
        # define the input molecule
        self.molecule = molecule
        
        # define the threshold of optimization
        self.threshold = threshold
        self.convergence = convergence
        self.tolerance = tolerance
        
        self.records_smiles = []
        self.records_vinas = []
        self.records_plogps = []
        self.records_sascores = []
        
        # define the parents size and children size
        self.min_parents_size = min_parents_size
        self.max_children_size = max_children_size
        
        # define the iters
        self.iters = iters
        
        # load the fitness and pareto front
        self.fitness = FITNESS(program_path, receptor_path)
        self.PF = Pareto_Front(objectives=('vina', 'plogp', 'sascores'), maximize=('plogp'))
        
        # load the fragment library
        df = pd.read_csv('./src/components/fragments-recap.txt')
        self.frags = []
        fragments = df['fragment'].tolist()
        for fragment in fragments:
            try:
                self.frags.append(Chem.MolFromSmiles(fragment))
            except Exception:
                continue
        
        # define the box center dim of molecular docking
        self.box_center_dim = box_center_dim
        
        self.output = output
        
        
    def selection(self, pop: list, folder_conf: str, folder_pdbqt: str, folder_output: str):
        # convert the rdkit mol from smiles
        pop_mol = (Chem.MolFromSmiles(ind) for ind in pop)
        
        # calculate the plogp and sascore from rdkit mols 
        smiles_list, plogps, sascores = self.fitness.run(pop_mol)
        
        # remove the dummy * in smiles and convert the smiles to rdkit mol
        smiles_list_rd = [remove_dummy(smiles) for smiles in smiles_list]
        smiles_list_mol = [Chem.MolFromSmiles(srd) for srd in smiles_list_rd]
        
        # calculate the vina score from rdkit mol
        if len(smiles_list_mol) == len(smiles_list):
            vinas, file_names = self.fitness.run_d(smiles_list_mol, self.box_center_dim, folder_conf, folder_pdbqt, folder_output)
        else:
            raise ValueError("Something Wrong!")
        
        # get the indexes
        file_indexes = [int(extract_number_split(file_name)) for file_name in file_names]
        
        # get the remained smiles and other properties
        remain_smiles = [smiles_list[i] for i in file_indexes]
        remain_plogps = [plogps[i] for i in file_indexes]
        remain_sascores = [sascores[i] for i in file_indexes]
        
        # create the pareto dataframe
        if len(remain_smiles) == len(vinas):
            data = pd.DataFrame({"smiles": remain_smiles, "vina": vinas, "plogp": remain_plogps, "sascores": remain_sascores})
        else:
            raise ValueError("Something wrong happened in selection")
        
        # Pareto_Front
        _, data_pareto = self.PF.pareto(data)
        parents = data_pareto['smiles'].values.tolist()
        vinas = data_pareto['vina'].values.tolist()
        plogps = data_pareto['plogp'].values.tolist()
        sascores = data_pareto['sascores'].values.tolist()
        
        # sample the parents
        if len(parents) < self.min_parents_size:
            ps = parents.copy()
            vs = vinas.copy()
            ls = plogps.copy()
            ss = sascores.copy()
            
            num_sample = self.min_parents_size - len(ps)
            
            data_sv = data.sort_values(by='vina', ascending=True)
            
            p_list1 = data_sv['smiles'].values.tolist()
            v_list1 = data_sv['vina'].values.tolist()
            l_list1 = data_sv['plogp'].values.tolist()
            s_list1 = data_sv['sascores'].values.tolist()
            
            p_list2 = []
            v_list2 = []
            l_list2 = []
            s_list2 = []                                         
            for p, v, l, s in zip(p_list1, v_list1, l_list1, s_list1):
                if p not in ps:
                    p_list2.append(p)
                    v_list2.append(v)
                    l_list2.append(l)
                    s_list2.append(s)
                else:
                    continue
                
                if len(p_list2) == num_sample:
                    break
                else:
                    pass
            
            parents = ps + p_list2
            vinas = vs + v_list2
            plogps = ls + l_list2
            sascores = ss + s_list2
            
            return parents, vinas, plogps, sascores
            
        else:
            
            return parents, vinas, plogps, sascores

    def connecting(self, parents_mol: list):
        mols = set()
        
        for coms in product(parents_mol, self.frags):
            nmols = connections(coms[0], coms[1])
            mols.update(nmols)

        return list(mols)
        
    def crossover(self, parents: list):
        """connecting the parents and predefine fragments."""
        parents_mol = []
        for parent in parents:
            if '*' in parent:
                parents_mol.append(Chem.MolFromSmiles(parent))
            else:
                continue
                
        child_list1 = []
        child_list2 = self.connecting(parents_mol)
        for child in child_list2:
            mc = mol_check(child)
            if mc:
                child_list1.append(Chem.MolToSmiles(child))
        
        if len(child_list1) > self.max_children_size:
            children = random.sample(child_list1, self.max_children_size)
        else:
            children = child_list1
        
        return children
    
    def recording(self, pop: list, vinas: list, plogps: list, sas: list):
        """recording great molecules"""
        for ind, vina, plogp, sa in zip(pop, vinas, plogps, sas):
            if ind in self.records_smiles:
                continue
            else:
                if vina <= self.threshold['vina'] and plogp >= self.threshold['plogp'] and sa <= self.threshold['sas']:
                    ind = remove_dummy(ind)
                    self.records_smiles.append(ind)
                    self.records_vinas.append(vina)
                    self.records_plogps.append(plogp)
                    self.records_sascores.append(sa)
                    
        
    def running(self, folder_conf, folder_pdbqt, folder_output):
        """running evo process"""
        # create init population
        pop = self.crossover([self.molecule])
        
        # define the init parameter
        early_stop = 0
        best_plogp = -10.0
        best_vina = 0.0
        best_sascore = 9.0
        
        # try to create the folders
        try:
            os.mkdir(folder_conf)
        except:
            pass
        
        try:
            os.mkdir(folder_pdbqt)
        except:
            pass
        
        try:
            os.mkdir(folder_output)
        except:
            pass        
        
        # Evo epoch
        for ep in tqdm(range(0, self.iters)):
            fc = f"{folder_conf}/output_{ep}"
            fp = f"{folder_pdbqt}/output_{ep}"
            fo = f"{folder_output}/output_{ep}"
            
            # select the parents
            parents, vinas, plogps, sas = self.selection(pop, folder_conf=fc, folder_pdbqt=fp, folder_output=fo)
            
            # calculate convergence
            convergence1 = np.max(plogps) - best_plogp
            convergence2 = np.min(vinas) - best_vina
            convergence3 = np.min(sas) - best_sascore
            
            if convergence1 > 0:
                best_plogp = np.max(plogps)
            
            if convergence2 < 0:
                best_vina = np.min(vinas)
                
            if convergence3 < 0:
                best_sascore = np.min(sas)
            
            self.recording(parents, vinas, plogps, sas)
            
            # check the convergence
            if convergence1 < self.convergence and -convergence2 < self.convergence and -convergence3 < self.convergence:
                early_stop += 1
                if early_stop > self.tolerance:
                    break
                else:
                    pass
            else:
                pass
            
            # check the molecular completety
            mcc = multi_conn_check(parents)
            if mcc:
                pass
            else:
                break
                
            # create the chidlren from parents
            children = self.crossover(parents)
            pop = parents + children
        
        
        df_final = pd.DataFrame(
            {'smiles': self.records_smiles, 'vinas': self.records_vinas, 'plogps': self.records_plogps, 'sascores': self.records_sascores}
        )

        df_final.to_csv(self.output, index=False)
