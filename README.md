# vaspToQE
Python script to generate QE-input files from vasp POSCAR files on a system with numpy installed.

## Advantages
- This script can generate Quantum Espresso input files directly from VASP POSCAR files. You just need to specify a psedupotential directory inside "ppPath.txt" file. 
- Can be highly useful to generate QE input files on the fly without a need to edit it for the atomic-weight and pseudopotential file description
- Create an alias on your working system and start generating QE input files from vasp POSCAR files and vice-versa.

## Usage
- python genQE.py -inputs
- Use, python genQE.py -h , for help regarding input parameters
