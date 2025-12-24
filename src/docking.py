import os
import subprocess

def config_gpu(receptor: str, ligands: str, output: str, box_center: list, box_size: list = None,
               num_modes: int = 10, thread: int = 8192):
    if box_size is None:
        box_size = [20, 20, 20]

    with open("config.txt", "w") as file:
        file.write(f"receptor = {receptor}\n")
        file.write(f"ligand_directory = {ligands}\n")
        file.write(f"output_directory = {output}\n")
        file.write(f"center_x = {box_center[0]}\n")
        file.write(f"center_y = {box_center[1]}\n")
        file.write(f"center_z = {box_center[2]}\n")
        file.write(f"size_x = {box_size[0]}\n")
        file.write(f"size_y = {box_size[1]}\n")
        file.write(f"size_z = {box_size[2]}\n")
        file.write(f"num_modes = {num_modes}\n")
        file.write(f"thread = {thread}")


def run_gpu(program_path: str, config_path: str):
    if not os.path.exists(program_path):
        raise FileNotFoundError(f"Docking software not found: {program_path}")

    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    command = [program_path, "--config", config_path]

    result = subprocess.run(command, capture_output=True, text=False)

    if result.returncode != 0:
        raise RuntimeError(f"The Molecular Docking on Cuda happen error：\n{result.stderr}")


def docking_gpu(program_path: str, receptor: str, ligands: str, box_center: list, outputs: str, box_size=None,
                num_modes: int = 10, thread: int = 8192):
    try:
        os.makedirs(outputs)
    except FileExistsError:
        pass

    config_gpu(
        receptor=receptor,
        ligands=ligands,
        output=outputs,
        box_center=box_center,
        box_size=box_size,
        num_modes=num_modes,
        thread=thread
    )

    run_gpu(
        program_path=program_path,
        config_path='config.txt'
    )


class MolecularDocking:
    def __init__(self, program: str, receptor: str, num_mode: int = 10):
        """
        >>> MD = MolecularDocking('program_path', 'protein.pdbqt', 10)
        >>> MD.run('./ligands', [0.0, 0.0, 0.0], 'output_path')
        """
        self.receptor = receptor
        self.num_mode = num_mode
        self.program = program

    def run(self, ligands: str, box_center: list, outputs: str, box_size=None, num_modes: int = 10, thread: int = 8192):
        """Running Molecular Docking with GPU"""
        docking_gpu(
            program_path=self.program,
            receptor=self.receptor,
            ligands=ligands,
            box_center=box_center,
            box_size=box_size,
            num_modes=num_modes,
            thread=thread,
            outputs=outputs
        )

        files = os.listdir(outputs)
        scores = []
        for file in files:
            f = f"{outputs}/{file}"
            scores.append(get_affinity(f))

        return scores, files

def get_affinity(file: str, pos=3):
    with open(file, "r") as f:
        lines = f.readlines()

    sl = lines[1].strip()

    return float(sl.split()[pos])

