import os
from openbabel import openbabel as ob

def genpdbqt(file: str, pdbqt: str = None):
    if not os.path.exists(file):
        raise FileNotFoundError(f"Input Molecule file not found: {file}!")

    file_name, file_ex = os.path.splitext(file)
    if pdbqt is None:
        pdbqt = file_name + ".pdbqt"
    else:
        pass

    if file_ex == ".pdb" or file_ex == ".mol2" or file_ex == ".sdf":
        pass
    else:
        raise ValueError("Input file format is not supported, please upload file with pdb, mol2 or sdf format!")

    conv = ob.OBConversion()
    conv.SetInAndOutFormats(file_ex[1:], "pdbqt")

    m = ob.OBMol()

    # read coordinate to mol
    if not conv.ReadFile(m, file):
        print("[!] Failed to convert input file!")

    m.SetAutomaticFormalCharge(False)
    m.SetAutomaticPartialCharge(False)

    # calculate Gasteiger charge
    charge_model = ob.OBChargeModel.FindType("gasteiger")
    if charge_model:
        charge_model.ComputeCharges(m)

    # write PDBQT with charge
    if not conv.WriteFile(m, pdbqt):
        raise RuntimeError("Failed to write PDBQT file")

    conv.CloseOutFile()


def multi_genpdbqt(file_list, folder):
    """ """
    try:
        os.mkdir(folder)
    except FileExistsError:
        pass

    for i, file in enumerate(file_list):
        path = f"{folder}/{i}.pdbqt"
        genpdbqt(file, path)


