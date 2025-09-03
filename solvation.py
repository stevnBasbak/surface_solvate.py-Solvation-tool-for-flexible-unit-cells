# Create a solvated surface with given lattice parameters and desired solvent molecules (pdb format)
# removes the solvent molecules that are closer then a distance_treshold from the surface
# removes the solvent molecules that have a certain mass percentage outside the unit cell
# Run with: $ python surface_solvate.py

###########################################################################################
#                                  REQUIRED INPUTS
#
# Surface Definition:
surface_pdb = 'diphonix_matrix_unitcell.pdb'  # PDB file for the solid surface.

# Solvent Properties:
solvent_molecule = "water.pdb"    # PDB format file for a molecule of Solvent 1.
density = 0.8                          # Density of Solvent 1 (grams per mL).
molar_mass = 18                      # Molar mass of Solvent 1 (grams per mole).
fraction_solvent_1 = 1                  # Molar fraction of Solvent 1 relative to other solvents (0-1).

# Distance Threshold:
distance_threshold = 1.0  # Minimum distance from the surface to the solvent layer (Angstroms) (0.95 recommended)
mass_threshold = 0  # Percentage of mass allowed outside the unit cell, Boundary condition

###########################################################################################
#                                  OPTIONAL INPUTS
#
#
# Solute Properties (Optional):
solute_molecule = "uranyl.pdb"  # PDB file for the solute molecule.
solute_concentration = 1               # Molar concentration of the solute.
#
solute_cation_molecule = "Li.pdb"       # PDB file for the solute molecule
solute_cation_concentration = 0        # Molar concentration of the solute.

# Additional Solvents:
solvent_molecule_2 = "PC.pdb"  # PDB format file for a molecule of Solvent 2.
density_2 = 1.2                         # Density of Solvent 2 (grams per mL).
molar_mass_2 = 102.09                   # Molar mass of Solvent 2 (grams per mole).
fraction_solvent_2 = 0.0                # Molar fraction of Solvent 2 relative to Solvent 1 (0-1).
#
solvent_molecule_3 = "DMC_molecule.pdb"  # PDB format file for Solvent 3.
density_3 = 1.07                         # Density of Solvent 3 (grams per mL).
molar_mass_3 = 90.08                     # Molar mass of Solvent 3 (grams per mole).
fraction_solvent_3 = 0.0                   # Molar fraction of Solvent 3 relative to Solvent 1 (0-1).
###########################################################################################




## part zero: checking the validity of the input paramters and files

import os
import shutil
import sys

# Define the folder name
folder_name = './temp_workfolder'
folder_path = os.path.join(os.getcwd(), folder_name)

# Check if the folder exists
if os.path.exists(folder_path) and os.path.isdir(folder_path):
    # Delete the folder and all of its contents
    shutil.rmtree(folder_path)
    print(f"Folder '{folder_name}' and all its contents have been deleted.")






## here we check if the surface_pdb file exists

if not os.path.exists(surface_pdb):
    raise FileNotFoundError(f"The file {surface_pdb} was not found.")


## this part checks wether all libraries are installed and shows how to install them if not

def check_and_instruct(library):
    try:
        __import__(library)
    except ImportError:
        print(f"{library} is not installed.")
        if library in custom_install_commands:
            print(f"To install {library}, run: '{custom_install_commands[library]}'")
        else:
            print(f"To install {library}, run: 'pip install {library}'")

# Custom installation commands for libraries that need special treatment or are not directly available via PyPI
custom_install_commands = {
    "openbabel": "conda install -c conda-forge openbabel",  # normally should be system wide installmeent
    "pybel": "pip install pybel"
}

libraries_to_check = [
    "numpy",
    "os",  # This is a built-in module
    "math",  # This is also a built-in module.
    "collections",  # Built-in module
    "networkx",
    "sys",  # Built-in module.
    "openbabel",
    "pybel"
]

for lib in libraries_to_check:
    check_and_instruct(lib)




import shutil
import subprocess

def check_packmol():
    # Check if packmol executable is in PATH
    packmol_path = shutil.which("packmol")
    if packmol_path is None:
        print("Error: 'packmol' is not found in your PATH.")
        print("Please install packmol and make sure it is available in your PATH.")
        print("To install packmol, follow these steps:")
        print("1. Download the source code from: http://m3g.iqm.unicamp.br/packmol/download.shtml")
        print("2. Extract the downloaded tar.gz file.")
        print("3. Navigate to the extracted directory in your terminal.")
        print("4. Run 'make' to compile packmol.")
        print("5. Move the compiled executable to a directory in your PATH, e.g., '/usr/local/bin/'.")
        sys.exit(1)

check_packmol()






# now that all libraries are present we can import them
import numpy as np
import os
from math import radians, sin, cos
import sys
from openbabel import openbabel
import pybel
import re
import numpy as np
from scipy.spatial import cKDTree


# just checking that the total fraction is one :)
def check_solvent_fractions():
    total_fraction = fraction_solvent_1 + fraction_solvent_2 + fraction_solvent_3
    assert total_fraction == 1, "Error: The sum of the solvent fractions must equal to 1."

check_solvent_fractions()



# This is also a check of correct input parameters; It checks that if the solvent fraction >0, then the file must exist
files_and_conditions = [
    (solvent_molecule, fraction_solvent_1),
    (solvent_molecule_2, fraction_solvent_2),
    (solvent_molecule_3, fraction_solvent_3),
    (solute_molecule, solute_concentration),
    (solute_cation_molecule, solute_cation_concentration)
]

# Function to check the existence of files
def check_file_existence(file_path, condition):
    if condition > 0:
        if not os.path.exists(file_path):
            print(f"Error: The required file '{file_path}' does not exist.")
            sys.exit(1)
        else:
            print(f"Check passed: The file '{file_path}' exists as required.")


# Run checks
for file_path, condition in files_and_conditions:
    check_file_existence(file_path, condition)


print(f"part 1/10 completed: All input paramters are correct")
print(f"\n")



### part one: creeer een solvent box met afmetingen die perfect passen op de unit cell, daarvoor moeten we eerste berekenen wat de maximale afstand van de unit cell is
### en ook hoeveel molecules er in de solvent box moeten zijn van elk solvent/solute

## Dus in deze laatste versie van het surface_solvate script maken we een doos die perfect past op de crooked unit cell toch?



# we starten met het kopieren van alle necessary files naar ./temp_workfolder waar alles uitgevoerd zal worden. Clean

# List of files to copy
files_to_copy = [
    surface_pdb, 
    solvent_molecule, 
    solute_molecule, 
    solute_cation_molecule, 
    solvent_molecule_2, 
    solvent_molecule_3
]

# Define the new folder path
new_folder = './temp_workfolder'

# If the folder exists, delete it and its contents
if os.path.exists(new_folder):
    shutil.rmtree(new_folder)

# Create the new folder
os.makedirs(new_folder)

# Copy the files to the new folder
for file_name in files_to_copy:
    if os.path.exists(file_name):
        shutil.copy(file_name, new_folder)

# Change the working directory to the new folder
os.chdir(new_folder)





# hier halen we de unit cell parameters uit de CRYST1 line van de surface_pdb
def read_unit_cell_parameters(surface_pdb):
    with open(surface_pdb, 'r') as file:
        for line in file:
            if line.startswith('CRYST1'):
                params = line.split()
                return {
                    'a': float(params[1]),
                    'b': float(params[2]),
                    'c': float(params[3]),
                    'alpha': float(params[4]),
                    'beta': float(params[5]),
                    'gamma': float(params[6])
                }
    return None

# Set the PDB file path as a variable

unit_cell = read_unit_cell_parameters(surface_pdb)

if unit_cell:
    print("Unit cell parameters:")
    for param, value in unit_cell.items():
        print(f"{param}: {value}")
    
    # You can use these parameters globally in your script
    # For example:
    a, b, c = unit_cell['a'], unit_cell['b'], unit_cell['c']
    alpha, beta, gamma = unit_cell['alpha'], unit_cell['beta'], unit_cell['gamma']
    
else:
    print("CRYST1 line not found in the PDB file.")
    sys.exit(1)







## hier maken we van deze unit cell een perfect passende cuboid box die uiteindelijk gesolvateerd zal worden. op deze manier zal er nooit te veel gesolvateerd
## worden aangezien enkel crooked unit cells te veel solvent zullen hebben die uitiendelijk erwijderd zullen worden

import numpy as np

# Convert angles to radians
alpha_rad = np.radians(alpha)
beta_rad = np.radians(beta)
gamma_rad = np.radians(gamma)

# Compute lattice vectors in Cartesian coordinates
a_vec = np.array([a, 0, 0])
b_vec = np.array([b * np.cos(gamma_rad), b * np.sin(gamma_rad), 0])
c_vec = np.array([
    c * np.cos(beta_rad),
    c * np.cos(alpha_rad) * np.sin(gamma_rad),
    c * np.sqrt(1 - np.cos(beta_rad)**2 - (np.cos(alpha_rad) * np.sin(gamma_rad))**2)
])

# Compute all corners of the unit cell
corners = [
    np.array([0, 0, 0]),
    a_vec,
    b_vec,
    c_vec,
    a_vec + b_vec,
    a_vec + c_vec,
    b_vec + c_vec,
    a_vec + b_vec + c_vec
]

# Convert corners to a numpy array for min/max computations
corners_array = np.array(corners)

# Calculate the extent along each Cartesian axis
global x_min, x_max, y_min, y_max, z_min, z_max
x_min, y_min, z_min = np.min(corners_array, axis=0)
x_max, y_max, z_max = np.max(corners_array, axis=0)

# Global variables for cuboid dimensions
global x_len, y_len, z_len
x_len = x_max - x_min
y_len = y_max - y_min
z_len = z_max - z_min

# Compute cuboid corners based on the extent
cuboid_corners = [
    [x_min, y_min, z_min],
    [x_max, y_min, z_min],
    [x_min, y_max, z_min],
    [x_max, y_max, z_min],
    [x_min, y_min, z_max],
    [x_max, y_min, z_max],
    [x_min, y_max, z_max],
    [x_max, y_max, z_max]
]

# Output results
print("Cuboid corners:")
for i, corner in enumerate(cuboid_corners, start=1):
    print(f"Corner {i}: {corner}")

print(f"x_min: {x_min:.3f}, x_max: {x_max:.3f}")
print(f"y_min: {y_min:.3f}, y_max: {y_max:.3f}")
print(f"z_min: {z_min:.3f}, z_max: {z_max:.3f}")
print(f"x_len: {x_len:.3f} Å")
print(f"y_len: {y_len:.3f} Å")
print(f"z_len: {z_len:.3f} Å")






volume = x_len * y_len * z_len

print(f"Volume cuboid solvent box: {volume} Å³")

print("part 2/10 completed. Crystallographic unit cell parameters extracted and used for creation of the cuboid solvent box")
print(f"\n")





## hieronder wordt nu aangenomen dat dit het volume is dat we gaan nemen voor de berekening van hoeveel moleculen er nodig zijn

# number of solvent molecules in the cubic box
def calculate_number_of_molecules(volume, density, molar_mass, fraction_solvent_1):
    num_molecules = ((volume * density * 0.6022) / molar_mass) * fraction_solvent_1
    return num_molecules

# number of solvent_2 molecules in the cubic box
def calculate_number_of_molecules_solvent_2(volume, density_2, molar_mass_2, fraction_solvent_2):
    num_molecules_2 = ((volume * density_2 * 0.6022) / molar_mass_2) * fraction_solvent_2
    return num_molecules_2

# number of solvent_3 molecules in the cubic box
def calculate_number_of_molecules_solvent_3(volume, density_3, molar_mass_3, fraction_solvent_3):
    num_molecules_3 = ((volume * density_3 * 0.6022) / molar_mass_3) * fraction_solvent_3
    return num_molecules_3


# number of solute molecules based on the molar concentration
def calculate_number_of_solute_molecules(volume, solute_concentration):
    num_molecules_solute = (volume * solute_concentration * 6.022 * 1e-4)
    return int(num_molecules_solute)


# number of solute cation molecules based on the molar concentration
def calculate_number_of_solute_cation_molecules(volume, solute_cation_concentration):
    num_molecules_solute_cation = (volume * solute_cation_concentration * 6.022 * 1e-4)
    return int(num_molecules_solute_cation)







## nu wordt de packmol input.inp file geschreven op basis van de berekende hoeveelheid moleculen en afstanden van de unit cell

# this part defines the actual input.inp file used by packmol
def generate_packmol_input(solvent_molecule, num_molecules, solvent_molecule_2, num_molecules_2, solvent_molecule_3, num_molecules_3, solute_molecule, num_molecules_solute, solute_cation_molecule, num_molecules_solute_cation, x_len, y_len, z_len):
    with open('input.inp', 'w') as f:
        f.write("tolerance 2.0\n")
        f.write(f"filetype pdb\n")
        f.write(f"output solvent_box.pdb\n\n")

        f.write(f"structure {solvent_molecule}\n")
        f.write(f"  number {num_molecules}\n")
        f.write(f"  inside box {x_min} {y_min} {z_min} {x_max} {y_max} {z_max}\n")
        print(f"Solute molecule '{solvent_molecule}' found. Added to Packmol input file.")
        f.write(f"end structure\n")
        
# solvent 2 (optional)
        if os.path.exists(solvent_molecule_2) and num_molecules_2 and fraction_solvent_2 > 0:
            f.write("\n")
            f.write(f"structure {solvent_molecule_2}\n")
            f.write(f"  number {num_molecules_2}\n")
            f.write(f"  inside box {x_min} {y_min} {z_min} {x_max} {y_max} {z_max}\n")
            print(f"Solute molecule '{solvent_molecule_2}' found. Added to Packmol input file.")
            f.write(f"end structure\n")

# solvent 3 (optional)
        if os.path.exists(solvent_molecule_3) and num_molecules_3 and fraction_solvent_3 > 0:
            f.write("\n")
            f.write(f"structure {solvent_molecule_3}\n")
            f.write(f"  number {num_molecules_3}\n")
            f.write(f"  inside box {x_min} {y_min} {z_min} {x_max} {y_max} {z_max}\n")
            print(f"Solute molecule '{solvent_molecule_3}' found. Added to Packmol input file.")
            f.write(f"end structure\n")

# solute (optional)
        if os.path.exists(solute_molecule) and num_molecules_solute > 0:
            f.write("\n")
            f.write(f"structure {solute_molecule}\n")
            f.write(f"  number {num_molecules_solute}\n")
            f.write(f"  inside box {x_min} {y_min} {z_min} {x_max} {y_max} {z_max}\n")
            f.write("end structure\n")
            print(f"Solute molecule '{solute_molecule}' found. Added to Packmol input file.")

# solute cation (optional)
        if os.path.exists(solute_cation_molecule) and num_molecules_solute_cation > 0:
            f.write("\n")
            f.write(f"structure {solute_cation_molecule}\n")
            f.write(f"  number {num_molecules_solute_cation}\n")
            f.write(f"  inside box {x_min} {y_min} {z_min} {x_max} {y_max} {z_max}\n")
            f.write("end structure\n")
            print(f"Solute molecule '{solute_cation_molecule}' found. Added to Packmol input file.")


## hierboven hebben we iegenlijk allemaal defenities gemaakt maar nu gaan we ze uitvoeren

# number of solvent molecules (solvent 1)
num_molecules = calculate_number_of_molecules(volume, density, molar_mass, fraction_solvent_1)

# number of solvent_2 molecules
num_molecules_2 = calculate_number_of_molecules_solvent_2(volume, density_2, molar_mass_2, fraction_solvent_2)

# number of solvent_3 molecules
num_molecules_3 = calculate_number_of_molecules_solvent_3(volume, density_3, molar_mass_3, fraction_solvent_3)

# number of solute molecules
num_molecules_solute = calculate_number_of_solute_molecules(volume, solute_concentration)

# number of solute cation molecules
num_molecules_solute_cation = calculate_number_of_solute_molecules(volume, solute_cation_concentration)

# write packmol file based on above defenitions
generate_packmol_input(solvent_molecule, num_molecules, solvent_molecule_2, num_molecules_2, solvent_molecule_3, num_molecules_3, solute_molecule, num_molecules_solute, solute_cation_molecule, num_molecules_solute_cation, x_len, y_len, z_len)






## now we check in the packmol inp file if every mentioned file has a pdb extention. Packmol cant work without :(
import sys

def check_structure_files(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        sys.stderr.write(f"Error: Input file '{filename}' not found\n")
        sys.exit(1)
    
    errors = []
    for line in lines:
        line = line.strip()
        if line.startswith('structure '):
            parts = line.split()
            if len(parts) < 2:
                continue
            fname = parts[1]
            if not fname.lower().endswith('.pdb'):
                errors.append(fname)
    
    if errors:
        sys.stderr.write("Error: Non-PDB files found in structure declarations:\n")
        for error in errors:
            sys.stderr.write(f"  - {error}\n")
        sys.stderr.write("All structure files must have .pdb extension\n")
        sys.exit(1)
    

if __name__ == "__main__":
    input_file = "input.inp"
    check_structure_files(input_file)


### make the amount of solvent/solute molecules a natural number in the packmol input file
## maybe this step is not necessary

def round_molecules(filename):
    # open and read the file
    with open(filename, 'r') as file:
        content = file.read()

    # find all instances of "number X" and round X to the nearest natural number
    content_modified = re.sub(r'(\bnumber\s+)(\d+(\.\d+)?)', lambda match: f'{match.group(1)}{round(float(match.group(2)))}', content)

    # Write the modified content back to the file
    with open(filename, 'w') as file:
        file.write(content_modified)

if __name__ == "__main__":
    round_molecules('input.inp')

print(f"part 3/10 completed: packmol input file generated")
print('\n')

print(f"now running packmol < input.inp ...")




import subprocess

def run_packmol(input_file, show_output=False):
    command = f"packmol < {input_file}"
    try:
        # Run packmol with redirected output
        result = subprocess.run(
            command,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True  # This returns strings instead of bytes
        )
        
        # Only show output if requested
        if show_output:
            print("Packmol Output:\n", result.stdout)
            if result.stderr:
                print("Stderr:\n", result.stderr)
        else:
            print("Packmol ran successfully")
            
        return True
        
    except subprocess.CalledProcessError as e:
        print("Error running Packmol")
        print("Error message:\n", e.stderr)
        return False

if __name__ == "__main__":
    input_file = "input.inp"
    # Set show_output=True if you want to see packmol's output
    success = run_packmol(input_file, show_output=False)



if os.path.exists(solvent_molecule) and num_molecules and fraction_solvent_1 > 0:
    print(f"Number of {solvent_molecule} molecules: {round(num_molecules)}")

if os.path.exists(solvent_molecule_2) and num_molecules_2 and fraction_solvent_2 > 0:
    print(f"Number of {solvent_molecule_2} molecules: {round(num_molecules_2)}")

if os.path.exists(solvent_molecule_3) and num_molecules_3 and fraction_solvent_3 > 0:
    print(f"Number of {solvent_molecule_3} molecules: {round(num_molecules_3)}")

if os.path.exists(solute_molecule) and num_molecules_solute > 0:
    print(f"Number of {solute_molecule} molecules: {round(num_molecules_solute)}")

if os.path.exists(solute_cation_molecule) and num_molecules_solute_cation > 0:
    print(f"Number of {solute_cation_molecule} molecules: {round(num_molecules_solute_cation)}")

print(f"part 4/10 completed: packmol ran succesfully")
print(f"\n")

#
#
#

#
#
#



print(f"Now erasing overlapping solvent molecules ...")
## part 4: Dit is een groot deel van het script. Nu alle coordinaten van beide files gefixt zijn gaan we nu de merger uitvoeren

import re

def fix_conect_lines(input_pdb_file, output_pdb_file):
    def fix_line(line):
        if line.startswith('CONECT'):
            prefix = line[:6]
            data = line[6:].strip()
            
            # Add a space after every 5 digits if there are no spaces already
            fixed_data = re.sub(r'(\d{5})(?=\d)', r'\1 ', data)
            
            # Combine the prefix and the fixed data part
            fixed_line = f"{prefix} {fixed_data}\n"
            return fixed_line
        else:
            return line
    
    with open(input_pdb_file, 'r') as infile, open(output_pdb_file, 'w') as outfile:
        for line in infile:
            outfile.write(fix_line(line))

# Usage example:
input_pdb_file = './solvent_box.pdb'
output_pdb_file = './solvent_box_conect.pdb'
fix_conect_lines(input_pdb_file, output_pdb_file)


# de eerste stap die we nemen is het verwijderen van solvent moleculen die te dicht in de buurt zitten van de surface_pdb atoms
# de manier waarop dit werkt is door te kijken naar de CONECT lines om te bepalen of de molecules aan elkaar hangen
# dit stuk is gechekt en werkt goed

import numpy as np
from scipy.spatial import cKDTree


def read_pdb(filename):
    atoms = []
    connections = {}
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('HETATM'):
                atom_id = int(line[6:11])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append((atom_id, (x, y, z)))
            elif line.startswith('CONECT'):
                parts = line.split()
                atom_id = int(parts[1])
                connections[atom_id] = [int(p) for p in parts[2:]]
    return atoms, connections

def get_connected_atoms(atom_id, connections):
    connected = set([atom_id])
    to_check = set(connections.get(atom_id, []))
    while to_check:
        current = to_check.pop()
        if current not in connected:
            connected.add(current)
            to_check.update(connections.get(current, []))
    return connected

def minimum_image_distance(vec1, vec2, box_lengths):
    delta = np.array(vec1) - np.array(vec2)
    for i in range(3):
        delta[i] -= np.round(delta[i] / box_lengths[i]) * box_lengths[i]
    return np.linalg.norm(delta)

def remove_close_molecules(solvent_file, surface_file, output_file, distance_threshold):
    # Read atoms and connections from the files
    solvent_atoms, solvent_connections = read_pdb(solvent_file)
    surface_atoms, _ = read_pdb(surface_file)

    # Extract coordinates
    solvent_coords = np.array([atom[1] for atom in solvent_atoms])
    surface_coords = np.array([atom[1] for atom in surface_atoms])

    # Define the box lengths for periodic boundaries
    box_lengths = np.array([a, b, c])

    # Identify solvent atoms to remove
    atoms_to_remove = set()
    for i, solvent_atom in enumerate(solvent_coords):
        for surface_atom in surface_coords:
            dist = minimum_image_distance(solvent_atom, surface_atom, box_lengths)
            if dist < distance_threshold:
                atom_id = solvent_atoms[i][0]
                atoms_to_remove.update(get_connected_atoms(atom_id, solvent_connections))
                break

    # Write the output file, excluding atoms to be removed
    with open(solvent_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('HETATM'):
                atom_id = int(line[6:11])
                if atom_id not in atoms_to_remove:
                    outfile.write(line)
            elif line.startswith('CONECT'):
                parts = line.split()
                atom_id = int(parts[1])
                if atom_id not in atoms_to_remove:
                    connected = [int(p) for p in parts[2:] if int(p) not in atoms_to_remove]
                    if connected:
                        outfile.write(f"CONECT {atom_id:5d}" + "".join(f"{c:5d}" for c in connected) + "\n")
            else:
                outfile.write(line)

# Usage
solvent_file = './solvent_box_conect.pdb'
surface_file = surface_pdb
output_file = './hollowed_solvent_box.pdb'

remove_close_molecules(solvent_file, surface_file, output_file, distance_threshold)


print(f"part 5/10 completed: Solvent molecules with an atom closer then {distance_threshold} Å have been deleted ")
print(f"\n")



print(f"Processing PDB file with {mass_threshold}% mass threshold ...")

def clean_conect_lines(file_path):
    # Open the file and read all lines
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Process each line
    cleaned_lines = []
    for line in lines:
        if line.startswith("CONECT"):
            # Split the line into elements and join them with a single space
            cleaned_line = " ".join(line.split())
            cleaned_lines.append(cleaned_line + "\n")
        else:
            cleaned_lines.append(line)

    # Write the cleaned lines back to the file
    with open(file_path, 'w') as file:
        file.writelines(cleaned_lines)

# Specify the path to the PDB file
file_path = './hollowed_solvent_box.pdb'
clean_conect_lines(file_path)




import re

def fix_conect_lines(input_pdb_file, output_pdb_file):
    def fix_line(line):
        if line.startswith('CONECT'):
            prefix = line[:6]
            data = line[6:].strip()
            
            # Add a space after every 5 digits if there are no spaces already
            fixed_data = re.sub(r'(\d{5})(?=\d)', r'\1 ', data)
            
            # Combine the prefix and the fixed data part
            fixed_line = f"{prefix} {fixed_data}\n"
            return fixed_line
        else:
            return line
    
    with open(input_pdb_file, 'r') as infile, open(output_pdb_file, 'w') as outfile:
        for line in infile:
            outfile.write(fix_line(line))

# Usage example:
input_pdb_file = './hollowed_solvent_box.pdb'
output_pdb_file = './hollowed_solvent_box_conect.pdb'
fix_conect_lines(input_pdb_file, output_pdb_file)



## nu halen we de solvent molecules eruit die uit de unit cell liggen, met een massa/percentage cutoff variable gedefineerd door de user
# oiutput is gecheckt en werkt

import numpy as np
from math import cos, sin, radians

# Extended dictionary of atomic masses (in atomic mass units)
ATOMIC_MASSES = {
    'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.811, 'C': 12.011,
    'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305,
    'Al': 26.982, 'Si': 28.086, 'P': 30.974, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948,
    'K': 39.098, 'Ca': 40.078, 'Fe': 55.845, 'Cu': 63.546, 'Zn': 65.38, 'Br': 79.904,
    'I': 126.904, 'Ba': 137.327, 'Au': 196.967
}

def cart_to_frac(cart_coords, cell_params):
    a, b, c, alpha, beta, gamma = cell_params
    alpha, beta, gamma = map(radians, (alpha, beta, gamma))
    
    v = a * b * c * np.sqrt(1 - cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2 + 2*cos(alpha)*cos(beta)*cos(gamma))
    
    transform_matrix = np.array([
        [1/a, -cos(gamma)/(a*sin(gamma)), b*c*(cos(alpha)*cos(gamma) - cos(beta))/(v*sin(gamma))],
        [0, 1/(b*sin(gamma)), a*c*(cos(beta)*cos(gamma) - cos(alpha))/(v*sin(gamma))],
        [0, 0, a*b*sin(gamma)/v]
    ])
    
    return np.dot(cart_coords, transform_matrix.T)

def is_inside_unit_cell(frac_coords):
    return all(0 <= coord < 1 for coord in frac_coords)

def parse_pdb(filename):
    atoms = []
    connections = {}
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('HETATM'):
                atom_id = int(line[6:11])
                element = line[76:78].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append((atom_id, element, np.array([x, y, z])))
            elif line.startswith('CONECT'):
                parts = line.split()
                atom_id = int(parts[1])
                connections[atom_id] = [int(p) for p in parts[2:]]
    
    return atoms, connections

def find_molecules(atoms, connections):
    molecules = []
    visited = set()
    
    def dfs(atom_id):
        molecule = set()
        stack = [atom_id]
        while stack:
            current = stack.pop()
            if current not in visited:
                visited.add(current)
                molecule.add(current)
                for neighbor in connections.get(current, []):
                    if neighbor not in visited:
                        stack.append(neighbor)
        return molecule
    
    for atom_id, _, _ in atoms:
        if atom_id not in visited:
            molecule = dfs(atom_id)
            molecules.append(molecule)
    
    return molecules

def calculate_mass_inside(molecule, atoms, inside_atoms):
    total_mass = 0
    inside_mass = 0
    for atom_id in molecule:
        atom = next((a for a in atoms if a[0] == atom_id), None)
        if atom is None:
            continue  # Skip if atom_id is not found in atoms
        mass = ATOMIC_MASSES.get(atom[1], 0)
        total_mass += mass
        if atom_id in inside_atoms:
            inside_mass += mass
    return inside_mass, total_mass


def remove_outside_molecules(pdb_file, output_file, cell_params, mass_threshold):
    atoms, connections = parse_pdb(pdb_file)
    
    inside_atoms = set()
    for atom_id, _, coords in atoms:
        frac_coords = cart_to_frac(coords, cell_params)
        if is_inside_unit_cell(frac_coords):
            inside_atoms.add(atom_id)
    
    molecules = find_molecules(atoms, connections)
    
    kept_molecules = []
    for molecule in molecules:
        inside_mass, total_mass = calculate_mass_inside(molecule, atoms, inside_atoms)
        if inside_mass / total_mass >= (1 - mass_threshold / 100):
            kept_molecules.append(molecule)
    
    with open(pdb_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('HETATM'):
                atom_id = int(line[6:11])
                if any(atom_id in mol for mol in kept_molecules):
                    outfile.write(line)
            elif line.startswith('CONECT'):
                parts = line.split()
                atom_id = int(parts[1])
                if any(atom_id in mol for mol in kept_molecules):
                    outfile.write(line)
            else:
                outfile.write(line)

# Example usage
if __name__ == "__main__":
    cell_params = (a, b, c, alpha, beta, gamma)
    remove_outside_molecules('./hollowed_solvent_box_conect.pdb', './filtered_solvent_inside.pdb', cell_params, mass_threshold)


print(f"part 6/10 completed: Solvent molecules with a mass percentage higher then {mass_threshold} outside the unit cell have been removed")
print(f"\n")
# de eerste stap die we nemen is het verwijderen van solvent moleculen die te dicht in de buurt zitten van de surface_pdb atoms
# de manier waarop dit werkt is door te kijken naar de CONECT lines om te bepalen of de molecules aan elkaar hangen









## nu hebben we simpelweg twee pdb files die klaar zijn om gemerged te worden
## is het niet het simpelste om te mergen met openbabel? zo zal alles wel behouden blijven
from openbabel import openbabel, pybel

# Load the first PDB file
mol1 = next(pybel.readfile("pdb", "./filtered_solvent_inside.pdb"))

# Load the second PDB file
mol2 = next(pybel.readfile("pdb", surface_pdb))

# Combine the two molecules
merged_mol = pybel.Molecule(openbabel.OBMol())
merged_mol.OBMol += mol1.OBMol
merged_mol.OBMol += mol2.OBMol

# Write the merged molecule to a new PDB file
merged_mol.write("pdb", "./merged_output.pdb", overwrite=True)

print(f"part 7/10 completed: {surface_pdb} and solvent_box have been merged")
print('\n')
# ziet er precies wel deftig uit




## nu voegen we de CRYST1 line toe aan de derde rij indien deze al niet aanwezig was

def check_and_add_cryst1_line(output_merged_path, surface_pdb):
    # Check if the output_merged.pdb file exists
    if not os.path.isfile(output_merged_path):
        print(f"Error: The file {output_merged_path} does not exist.")
        return
    
    # Check if the surface_pdb file exists
    if not os.path.isfile(surface_pdb):
        print(f"Error: The file {surface_pdb} does not exist.")
        return
    
    # Read the contents of output_merged.pdb
    with open(output_merged_path, 'r') as f:
        lines = f.readlines()

    # Check if there's a line starting with 'CRYST1'
    has_cryst1 = any(line.startswith('CRYST1') for line in lines)

    if not has_cryst1:
        # Read the surface_pdb file to find the CRYST1 line
        with open(surface_pdb, 'r') as f:
            surface_lines = f.readlines()

        cryst1_line = next((line for line in surface_lines if line.startswith('CRYST1')), None)

        if cryst1_line:
            # Insert the CRYST1 line at the third position
            if len(lines) >= 2:
                lines.insert(2, cryst1_line)
            else:
                lines.append(cryst1_line)

            # Write the modified contents back to output_merged.pdb
            with open(output_merged_path, 'w') as f:
                f.writelines(lines)
            
            print(f"CRYST1 line added to {output_merged_path} at the third position.")
        else:
            print(f"Error: No CRYST1 line found in {surface_pdb}.")
    else:
        print(f"The file {output_merged_path} already contains a CRYST1 line.")

# Example usage
output_merged_path = './merged_output.pdb'
check_and_add_cryst1_line(output_merged_path, surface_pdb)






# nu gaan we de output file naar de parent directory brengen en verwijderen we de temp folder
def rename_and_move_file(surface_pdb):
    base_name = os.path.splitext(os.path.basename(surface_pdb))[0]
    
    new_file_name = f"{base_name}_solvated.pdb"

    old_file_name = "./merged_output.pdb"
    
    # Check if the old file exists
    if os.path.exists(old_file_name):
        # Rename the file
        os.rename(old_file_name, new_file_name)
        print(f"Renamed '{old_file_name}' to '{new_file_name}'")
        
        # Define the parent directory
        parent_directory = os.path.dirname(os.getcwd())
        
        # Move the new file to the parent directory
        shutil.move(new_file_name, os.path.join(parent_directory, new_file_name))
        print(f"Moved '{new_file_name}' to '{parent_directory}'")
        
        # Get the current working directory
        current_directory = os.getcwd()
        
        # Remove the current working directory
        os.chdir(parent_directory)  # Change to parent directory to avoid being in the directory to be removed
        shutil.rmtree(current_directory)
        print(f"Removed the directory '{current_directory}'")
    else:
        print(f"File '{old_file_name}' does not exist")

rename_and_move_file(surface_pdb)





print(f"part 10/10 completed: output file moved to the current directory and removed the temp_workfolder")

# Print the result
print('\n')
print('script ended succesfully!')

# that's all folks
# Created by Stijn De Vos
# For questions you can always email me x Stijn.De.Vos@vub.be