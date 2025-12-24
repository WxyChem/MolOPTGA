import argparse
import pandas as pd
from src.evo import MolecularOPTGA
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

def main():
    parser = argparse.ArgumentParser(description='MolOPT with RDKit and OPENBABEL')
    parser.add_argument('-i', '--input', required=True, type=str, help='Input smiles file')
    parser.add_argument('-r', '--receptor', required=True, type=str, help='The protein file')
    
    parser.add_argument('-c', '--conf', required=True, type=str, help='Save PATH for the PDB file')
    parser.add_argument('-p', '--pdbqt', required=True, type=str, help='Save PATH for the PDBQT file')
    parser.add_argument('-o', '--output', required=True, type=str, help='Save PATH for the PDBQT OUT file')
    
    parser.add_argument('-s', '--software', required=True, type=str, help='The docking software to use')
    
    parser.add_argument('-c_x', '--center_x', required=True, type=str, help='The center of box on x dim')
    parser.add_argument('-c_y', '--center_y', required=True, type=str, help='The center of box on y dim')
    parser.add_argument('-c_z', '--center_z', required=True, type=str, help='The center of box on z dim')
    
    parser.add_argument('-ps', '--parents_size', required=True, type=str, help='The docking software to use')
    parser.add_argument('-cs', '--children_size', required=True, type=str, help='The docking software to use')
    parser.add_argument('-t', '--iters', required=True, type=str, help='The docking software to use')
    args = parser.parse_args()
    
    file1 = args.input
    receptor = args.receptor
    
    folder1 = args.conf
    folder2 = args.pdbqt
    folder3 = args.output
    
    program = args.software
    
    center_x = float(args.center_x)
    center_y = float(args.center_y)
    center_z = float(args.center_z)
    
    parents_size = int(args.parents_size)
    children_size = int(args.children_size)
    
    iters = int(args.iters)
    
    df = pd.read_csv(file1)
    molecule = df['smiles'].tolist()[0]
    
    MOGA = MolecularOPTGA(molecule, 
                          box_center_dim=[center_x, center_y, center_z], 
                          program_path=program, 
                          receptor_path=receptor, 
                          min_parents_size=parents_size, 
                          max_children_size=children_size, 
                          iters=iters
                          )
                          
    MOGA.running(folder1, folder2, folder3)
    print("Done!")


if __name__ == "__main__":
    main()

