#!/usr/bin/env python3

import os
import sys
import numpy as np
import argparse

### Check the input arguments
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--VASP_file", help="The name of the VASP File", required=False, default="POSCAR")
parser.add_argument("-o", "--QE_file", help="The name of the QE output File", required=False, default="DEFAULT")
parser.add_argument("-t", "--QE_calType", help="Type of QE file (e.g, scf, nscf, dos, etc)", required=False, default="scf")
parser.add_argument("-p", "--PressureKb", help="Pressure for vc-relax in Kbar", required=False, default=0)
parser.add_argument("-c", "--CoordType", help="Coordinate unit (Cartesian/Fractional)", required=False, default="Cartesian")
args = parser.parse_args()

iFile, cal, cType, press, oFile = args.VASP_file, args.QE_calType, args.CoordType, args.PressureKb, args.QE_file

#Get pp path
try:
    with open("ppPath.txt", 'r') as f0:
        lines = f0.readlines()
        for l in lines:
                ls = l.split()
                if "=" not in ls[-1]:
                    ppPath = ls[-1]
                    break
                else:
                    ppPath = ls[-1].split("=")[-1]
                    break
except:
    ppPath = "./"

#class to extract the relevant information from vasp POSCAR
class cellInfo:

    def __init__(self, File):
        self.File = File

    # Getting all the lines of the POSCAR file
    def Data(self):
        with open(self.File, 'r') as f:
            lines = f.readlines()
            return lines

    def pref(self):
        line0 = "Unknown"
        lines = self.Data()
        if self.File.endswith(".in") or self.File.endswith(".pwscf"):
            for l in lines:
                if "pref" in l:
                    line0 = l.split()[2]
        elif self.File.endswith(".o"):
            pass
        else:
            line0 = lines[0].split()[0]
        return line0.strip("'")


    def coordType(self):
        lines = self.Data()
        cType = "angstrom"
        if self.File.endswith(".in") or self.File.endswith(".pwscf"):
            for i, l in enumerate(lines):
                if "ATOMIC_POSITIONS" in l:
                    if "crystal" in l:
                        cType = "crystal"
                    else:
                        cType = "angstrom"
                    break

        elif self.File.endswith(".o"):
            s0 = "Begin final coordinates"
            s1 = "End final coordinates"
            index0 = 0
            for i, l in enumerate(lines):
                if s0 in l:
                    index0 = i
                    break
            if "angstrom" in lines[index0+9]:
                cType = "angstrom"
            elif "bohr" in lines[index0+9]:
                cType = "bohr"
            else:
                cType = "crystal"

        else:
            Type = lines[7].split()[0]
            if Type.lower().startswith("s"):
                Type = lines[8].split()[0]
                if Type.lower().startswith("d"):
                    cType = "crystal"
                else:
                    cType = "angstrom"
            else:
                if Type.lower().startswith("d"):
                    cType = "crystal"
                else:
                    cType = "angstrom"
        return cType



    def cell(self):
        lines = self.Data()
        Cell = []
        if self.File.endswith(".in") or self.File.endswith(".pwscf"):
            cell = []
            for i, l in enumerate(lines):
                if "CELL_PARAMETERS" in l:
                    cell = lines[i+1:i+4]
            for v in cell:
                Cell.append([float(a) for a in v.split()])

        elif self.File.endswith(".o"):
            s0 = "Begin final coordinates"
            s1 = "End final coordinates"
            index0 = 0
            for i, l in enumerate(lines):
                if s0 in l:
                    index0 = i
                    break
            for l1 in lines[index0+5:index0+8]:
                if self.coordType() == "bohr":
                    Cell.append([float(a)*0.529177 for a in l1.split()])
                else:
                    Cell.append([float(a) for a in l1.split()])

        else:
            scale = float(lines[1].split()[0])
            for l in lines[2:5]:
                ls = [scale * float(a) for a in l.split()]
                Cell.append(ls)
        return np.array(Cell)


    def nType(self):
        lines = self.Data()
        ntyp = 0
        if self.File.endswith(".in") or self.File.endswith(".pwscf"):
            for i, l in enumerate(lines):
                if "ntyp" in l:
                    ntyp = int(l.split()[-1])
                    break
        elif self.File.endswith(".o"):
            for i, l in enumerate(lines):
                if "number of atomic types" in l:
                    ntyp = int(l.split()[-1])
        else:
            ntyp = len(lines[5].split())
        return ntyp



    def Atoms(self):
        lines = self.Data()
        newLines = []
        Atoms = []
        if self.File.endswith(".in") or self.File.endswith(".pwscf"):
            for i, l in enumerate(lines):
                if "ATOMIC_SPECIES" in l:
                    Atoms = [a.split()[0] for a in lines[i+1:i+self.nType()+1]]
                    break
        elif self.File.endswith(".o"):
            for i, l in enumerate(lines):
                if "atomic species   valence    mass     pseudopotential" in l:
                    if len(newLines) > 0:
                        continue
                    else:
                        newLines = lines[i+1:i+1+self.nType()]
                    break

            for nl in newLines:
                ls = nl.split()
                #Wts.append(ls[2])
                Atoms.append(ls[0])
        else:
            Atoms = lines[5].split()

        properAtoms = []
        for a in Atoms:
            if a[-1].isdigit():
                properAtoms.append(a[:-1])
            else:
                properAtoms.append(a)
        return Atoms


    def nTotal(self):
        lines = self.Data()
        nat = 0
        if self.File.endswith(".in") or self.File.endswith(".pwscf"):
            for i, l in enumerate(lines):
                if "nat" in l:
                    nat = int(l.split()[-1])
                    break
        elif self.File.endswith(".o"):
            for i, l in enumerate(lines):
                if "number of atoms/cell" in l:
                    nat = int(l.split()[-1])
                    break
        else:
            nat = sum([int(n) for n in lines[6].split()])
        return nat


    def nAtoms(self):
        lines = self.Data()
        nAtoms = []
        nDict = {}
        if self.File.endswith(".in") or self.File.endswith(".pwscf"):
            for i, l in enumerate(lines):
                if "ATOMIC_POSITIONS" in l:
                    atomicPosLines = lines[i+1:i+self.nTotal()+1]
                    break
            atomicLabels = [a.split()[0] for a in atomicPosLines]
            unique, counts = np.unique(np.array(atomicLabels), return_counts=True)
            for u, c in zip(unique, counts):
                nDict[u] = c
            for atom in self.Atoms():
                nAtoms.append(nDict[atom])
        elif self.File.endswith(".o"):
            s0 = "Begin final coordinates"
            s1 = "End final coordinates"
            for i, l in enumerate(lines):
                if s0 in l:
                    index0 = i
                elif s1 in l:
                    index1 = i
                    break
            posLines = lines[index0+10:index0+10+self.nTotal()]
            atomicLabels = [a.split()[0] for a in posLines]
            unique, counts = np.unique(np.array(atomicLabels), return_counts=True)
            for u, c in zip(unique, counts):
                nDict[u] = c
            for atom in self.Atoms():
                nAtoms.append(nDict[atom])
        else:
            nAtoms = [int(n) for n in lines[6].split()]
        return nAtoms




    def coords(self):
        lines = self.Data()
        coords = []
        if self.File.endswith(".in") or self.File.endswith(".pwscf"):
            C0 = []
            for i, l in enumerate(lines):
                if "ATOMIC_POSITIONS" in l:
                    C0 = lines[i+1:i+1+self.nTotal()]
                    break
            for l in C0:
                coords.append([float(a) for a in l.split()[1:4]])

        elif self.File.endswith(".o"):
            s0 = "Begin final coordinates"
            s1 = "End final coordinates"
            index0 = 0
            for i, l in enumerate(lines):
                if s0 in l:
                    index0 = i + 10
                    break
            C0 = lines[index0:index0+self.nTotal()]
            for l in C0:
                if self.coordType() == "bohr":
                    coords.append([float(a)*0.529177 for a in l.split()[1:4]])
                else:
                    coords.append([float(a) for a in l.split()[1:4]])

        else:
            base = 8
            Type = lines[7].split()[0]
            if Type.lower().startswith("s"):
                base = 9
            totalN = self.nTotal()
            for l in lines[base:base + totalN]:
                ls = [float(a) for a in l.split()[:3]]
                coords.append(ls)

        return coords


    def latPara(self):
        return np.array([np.linalg.norm(v) for v in self.cell()])

    def genKpts(self, Type='scf', kspacing=0.04):
        LV = np.reciprocal(self.latPara())
        if Type == "nscf":
            kspacing = kspacing / 2
        elif "fit" in Type.lower():
            kspacing = kspacing / 3

        kx, ky, kz = int(max(1, np.ceil(LV[0]/kspacing))), int(max(1, np.ceil(LV[1]/kspacing))),\
        int(max(1, np.ceil(LV[2]/kspacing)))
        minK = min([kx, ky, kz])
        kx, ky, kz = minK * (kx // minK), minK * (ky // minK), minK * (kz // minK)
        return [kx, ky, kz]


    def toCart(self):
        cartCoords = []
        if self.coordType() == "angstrom":
            cartCoords = self.coords()
        elif self.coordType() == "bohr":
            cartCoords = self.coords()
        else:
            lines = self.coords()
            tMatrix = self.cell().T # transpose of the cell
            for l in lines:
                cartCoords.append(np.matmul(tMatrix, np.array(l)))
        return cartCoords


    def toFrac(self):
        fracCoords = []
        if self.coordType() == "crystal":
            fracCoords = self.coords()
        else:
            lines = self.coords()
            tMatrix = np.linalg.inv(self.cell().T)
            for l in lines:
                fracCoords.append(np.matmul(tMatrix, np.array(l)))
        return fracCoords


    def volume(self):
        return np.linalg.det(self.cell())

    def atomicMass(self, atom="NA"):
        Elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds ', 'Rg ', 'Cn ', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

        Weights = [1.0070000000000001, 4.002, 6.941, 9.012, 10.811, 12.011, 14.007, 15.999, 18.998, 20.18, 22.99, 24.305, 26.982, 28.086, 30.974, 32.065, 35.453, 39.948, 39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693000000000005, 63.54600000000001, 65.38, 69.723, 72.64, 74.922, 78.96, 79.904, 83.79799999999999, 85.46799999999999, 87.62, 88.906, 91.22399999999999, 92.906, 95.96, 98.0, 101.07, 102.906, 106.42, 107.868, 112.411, 114.818, 118.71, 121.76, 127.6, 126.904, 131.293, 132.905, 137.327, 138.905, 140.116, 140.908, 144.24200000000002, 145.0, 150.36, 151.964, 157.25, 158.925, 162.5, 164.93, 167.25900000000001, 168.93400000000003, 173.054, 174.967, 178.49, 180.94799999999998, 183.84, 186.207, 190.23, 192.217, 195.084, 196.967, 200.59, 204.38299999999998, 207.2, 208.98, 210.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.03799999999998, 231.03599999999997, 238.02900000000002, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0, 258.0, 259.0, 262.0, 261.0, 262.0, 266.0, 264.0, 267.0, 268.0, 271.0, 272.0, 285.0, 284.0, 289.0, 288.0, 292.0, 295.0, 294.0]
        weights = {}
        for atom in self.Atoms():
            if atom in Elements:
                ind = Elements.index(atom)
                weights[atom] = Weights[ind]
            else:
                weights[atom] = "NA"
        return weights

####END####

#Function to generate kpoints

def genKpts(Cell, Type='scf', kspacing=0.04):
    latticePara = np.array([np.linalg.norm(v) for v in Cell])
    LV = np.reciprocal(latticePara)
    if Type == "nscf":
        kspacing = kspacing / 2
    elif "fit" in Type.lower():
        kspacing = kspacing / 3
    kx, ky, kz = int(max(1, np.ceil(LV[0]/kspacing))), int(max(1, np.ceil(LV[1]/kspacing))), int(max(1, np.ceil(LV[2]/kspacing)))
    minK = min([kx, ky, kz])
    kx, ky, kz = minK * (kx // minK), minK * (ky // minK), minK * (kz // minK)
    return [kx, ky, kz]

###BEGIN-qeToQE###
class qeToQE:

    def __init__(self, iFile="out.o"):
        self.iFile = iFile


    def latticeInfo(self):
        F1 = self.iFile
        nat = 0
        ntyp = 0
        pref = None
        AtomicSpeciesLines = []
        kPointsLines = []
        cellParaLines = []
        atomicPosLines = []
        cellAtomicLines = []
        if F1.endswith(".in"):
            with open(F1, 'r') as f1:
                lines = f1.readlines()
                for i, l in enumerate(lines):
                    if "pref" in l:
                        pref = l.split()[-1]
                    elif "nat" in l:
                        nat = int(l.split()[-1])
                    elif "ntyp" in l:
                        ntyp = int(l.split()[-1])
                    elif "ATOMIC_SPECIES" in l:
                        AtomicSpeciesLines = lines[i:i + ntyp + 1] + ["\n"]
                    elif "K_POINTS" in l:
                        kPointsLines = lines[i:i+2] + ["\n"]
                    elif "CELL_PARAMETERS" in l:
                        cellParaLines = lines[i:i+4]
                    elif "ATOMIC_POSITIONS" in l:
                        atomicPosLines = lines[i:i+nat+1]

        elif F1 == "out.o":
            s0 = "Begin final coordinates"
            s1 = "End final coordinates"
            Atoms, Wts = [], []
            pseudos = []
            index0, index1 = 0, 1
            newLines = []
            Cell = []
            with open(F1, 'r') as f1:
                lines = f1.readlines()
                for i, l in enumerate(lines):
                    if "number of atoms/cell" in l:
                        nat = int(l.split()[-1])
                    elif "number of atomic types" in l:
                        ntyp = int(l.split()[-1])
                    elif "atomic species   valence    mass     pseudopotential" in l:
                        if len(newLines) > 0:
                            continue
                        else:
                            newLines = lines[i+1:i+1+ntyp]
                        for nl in newLines:
                            ls = nl.split()
                            Wts.append(ls[2])
                            Atoms.append(ls[0])

                    elif "Pseudo file" in l:
                        pseudos.append(l.split()[2])
                    elif s0 in l:
                        index0 = i
                        for l1 in lines[i+5:i+8]:
                            Cell.append([float(a) for a in l1.split()])
                    elif s1 in l:
                        index1 = i
                cellAtomicLines = lines[index0+4:index1]
                #print(len(cellAtomicLines))
                #print(np.array(Cell))
                kp = genKpts(np.array(Cell), Type="fit")
                kPointsLines = ["\n", "K_POINTS automatic\n", f"{kp[0]} {kp[1]} {kp[2]}  0 0 0\n"]

            AtomicSpeciesLines.append("ATOMIC_SPECIES\n")
            for Atm, Wt, PP in zip(Atoms, Wts, pseudos):
                AtomicSpeciesLines.append(f"{Atm}  {float(Wt):>3f}   {PP}\n")


        cellLines = AtomicSpeciesLines + kPointsLines + cellParaLines + ["\n"] + atomicPosLines + cellAtomicLines
        return nat, ntyp, cellLines, Atoms
###END-qeToQE###



class vaspToQE(cellInfo):

    def __init__(self, File, Type, cType="default"):
        self.File = File
        self.Type = Type
        self.cType = cType

    def cellStr(self):
        Atoms = self.Atoms()
        Wts = self.atomicMass()
        atomicSpecies = ["\n", f"ATOMIC_SPECIES\n"]
        print(f'\n\tVASP STR {self.File} has \"{" ".join(Atoms)}\" Elements.\n')
        for a in Atoms:
            pp = qePseudo(a) # qePseudo function
            wt = f"{Wts[a]:.6f}"
            atomicSpecies.append(f"{a:<5}{wt:>8}   {pp:>8}\n")

        kp = self.genKpts(self.Type)
        kpts = ["\n", "K_POINTS automatic\n", f"{kp[0]} {kp[1]} {kp[2]}  0 0 0\n"]

        cellPara = ["\n", f"CELL_PARAMETERS {{angstrom}}\n"]
        for v in self.cell():
            v = [f"{a:.15f}" for a in v]
            cellPara.append(f"{v[0]:>20}{v[1]:>20}{v[2]:>20}\n")

        if self.cType.lower().startswith("c") and not self.coordType().startswith("a"):
            atomicPos = ["\n", f"ATOMIC_POSITIONS {{angstrom}}\n"]
            lines = self.toCart()
        elif self.cType.lower().startswith("f") and not self.coordType().startswith("c"):
            atomicPos = ["\n", f"ATOMIC_POSITIONS {{crystal}}\n"]
            lines = self.toFrac()
        else:
            atomicPos = ["\n", f"ATOMIC_POSITIONS {{{self.coordType()}}}\n"]
            lines = self.coords()

        Pos = []
        n0 = 0
        for n in self.nAtoms():
            Pos.append(lines[n0:n0 + n])
            n0 += n

        for i, coord in enumerate(Pos):
            for l in coord:
                l = [f"{a:.12f}" for a in l]
                atomicPos.append(f"{Atoms[i]:<4} {l[0]:>18}{l[1]:>18}{l[2]:>18}\n")

        cellLines = atomicSpecies + kpts + cellPara + atomicPos
        return cellLines

    def writeCell(self):
        lines = self.cellStr()
        with open("cellFile.txt", 'w') as f0:
            for l in lines:
                f0.write(l)


    def control_tag(self):
        Type = self.Type.lower()
        if "vc" in Type or "relax" in Type:
            dataDict = {"calculation":f"\'{self.Type}\'", "restart_mode":"\'from_scratch\'",
                    "prefix":f"\'{self.pref()}\'", "pseudo_dir":"\'./\'", "outdir":"\'./tmp/\'",
                    "forc_conv_thr":"1d-5"}

        elif "tc" in Type:
            dataDict = {"calculation":f"\'scf\'", "max_seconds":"8.64000e+04",
                    "prefix":f"\'Tc_calculation\'", "pseudo_dir":"\'./\'", "tprnfor":".TRUE.", "tstress":".TRUE."}

        else:
            if "fm" in Type:
                Type = "scf"
            dataDict = {"calculation":f"\'{Type}\'", "prefix":f"\'{self.pref()}\'",
                    "pseudo_dir":"\'./\'", "outdir":"\'./tmp/\'", "etot_conv_thr":"1d-5"}

        lines = ["&CONTROL\n"]
        for k, v in dataDict.items():
            lines.append(f"   {k:<18}{'='}  {v:<4}\n")
        return lines + ["/\n\n"]


    def system_tag(self):
        Type = self.Type.lower()

        if "fit" in Type:
            dataDict = {"constrained_magnetization":"\"none\"", "degauss":1.0000e-02, "ecutrho":700,
                    "ecutwfc":90, "ibrav":0, "nat":self.nTotal(), "ntyp":self.nType(),
                    "occupations":"\"smearing\"", "smearing":"\"mv\"", "la2F":".TRUE."}
        elif self.Type.lower().endswith("fm"):
            dataDict = {"ibrav":0, "degauss":0.01, "ecutrho":600,
                    "ecutwfc":80, "nat":self.nTotal(), "ntyp":self.nType(), "occupations":"\"smearing\"",
                    "smearing":"\"mv\"", "nspin":2}
            for l in range(self.nType()):
                dataDict[f"starting_magnetization({l+1})"] = 0
        else:
            dataDict = {"constrained_magnetization":"\"none\"", "degauss":1.0000e-02, "ecutrho":700,
                    "ecutwfc":90, "ibrav":0, "nat":self.nTotal(), "ntyp":self.nType(),
                    "occupations":"\"smearing\"", "smearing":"\"mv\""}

        lines = ["&SYSTEM\n"]
        for k, v in dataDict.items():
            lines.append(f"   {k:<26}{'='}  {v:<4}\n")
        return lines + ["/\n\n"]


    def magCellStr(self, mag="fm"):
        systemCard, atomicSpcies_pos = [], []
        nTotal, nType, cellLines, Atoms = self.nTotal(), self.nType(), self.cellStr(), self.Atoms()
        PosHead = len(cellLines)-nTotal - 1
        AtomicSpeciesLines = []
        AtomicPositions = cellLines[PosHead:]
        KPTS = []
        cellPara = []
        electronsCard = self.electrons_tag()
        for i, l in enumerate(cellLines):
            if "K_POINTS" in l:
                KPTS = cellLines[i-1:i+2]
            elif "CELL_PARAMETERS" in l:
                cellPara = cellLines[i-1:i+4]
            elif "ATOMIC_SPECIES" in l:
                AtomicSpeciesLines = cellLines[i-1:i+self.nType()+1]
        #_atomicSpecies = AtomicSpeciesLines[2:-1].copy()
        if self.Type.lower().endswith("afm"):
            AtomicPositions = [cellLines[PosHead]]
            AtomicSpeciesLines.clear()
            #nTotal, nType, cellLines, Atoms = qeToQE(self.File).latticeInfo()
            magMoments = {}
            print("\n\t" + "="*80)
            print("\n\tPlease Carefully ENTER Magnetic Parameters!\n")
            magAtom = input(f"\n\tPlease Choose Atom for AFM from {Atoms}: ")
            #nType += 1
            #Atoms.append(magAtom)
            for i in range(nType):
                if magAtom == Atoms[i]:
                    try:
                        magMoments[f"{Atoms[i]}1"] = input(f"\n\tInitial MagMoment for {Atoms[i]}: ")
                    except:
                        print(f"\tInvalid Input! DEFAULT value (0.00) is set.")
                        magMoments[f"{Atoms[i]}1"] = 0.0
                    magMoments[f"{Atoms[i]}2"] = "-" + magMoments[f"{Atoms[i]}1"]
                else:
                    try:
                        magMoments[f"{Atoms[i]}"] = input(f"\n\tInitial MagMoment for {Atoms[i]}: ")
                    except:
                        print(f"\tInvalid Input! DEFAULT value (0.00) is set.")
                        magMoments[f"{Atoms[i]}"] = 0.0

            dataDict = {"ibrav":0, "degauss":0.01, "ecutrho":600,
                    "ecutwfc":80, "nat":nTotal, "ntyp":nType+1, "occupations":"\"smearing\"",
                    "smearing":"\"mv\"", "nspin":2}
            for l in range(nType+1):
                dataDict[f"starting_magnetization({l+1})"] = magMoments[list(magMoments.keys())[l]]

            #Setting atomicPositions and atomicCoordinates
            #magAtom = input("\n\tPlease Choose Atom for AFM: ")
            for l in cellLines[:nType+2]:
                if magAtom in l and len(l.split()) == 3:
                    ls = l.split()
                    AtomicSpeciesLines.append(f"{magAtom}1  {ls[1]}   {ls[2]}\n")
                    AtomicSpeciesLines.append(f"{magAtom}2  {ls[1]}   {ls[2]}\n")
                else:
                    AtomicSpeciesLines.append(l)
            sign = "plus"
            for ind, l in enumerate(cellLines[PosHead + 1:]):
                ls = l.split()
                if magAtom in ls[0] and ind % 2 == 0:
                    atm = f"{ls[0]}1"
                    AtomicPositions.append(f"{atm:<4} {ls[1]:>18}{ls[2]:>18}{ls[3]:>18}\n")
                elif magAtom in ls[0] and ind % 2 != 0:
                    atm = f"{ls[0]}2"
                    AtomicPositions.append(f"{atm:<4} {ls[1]:>18}{ls[2]:>18}{ls[3]:>18}\n")
                else:
                    AtomicPositions.append(f"{ls[0]:<4} {ls[1]:>18}{ls[2]:>18}{ls[3]:>18}\n")


        elif self.Type == "fm":
            #nTotal, nType, cellLines, Atoms = qeToQE(self.File).latticeInfo()
            print("\n\t" + "="*80)
            print("\n\tPlease Carefully ENTER Magnetic Parameters!\n")
            magMoments = {}
            for i in range(nType):
                magMoments[Atoms[i]] = input(f"\n\tInitial MagMoment for {Atoms[i]}: ")
            dataDict = {"ibrav":0, "degauss":0.01, "ecutrho":600,
                    "ecutwfc":80, "nat":nTotal, "ntyp":nType, "occupations":"\"smearing\"",
                    "smearing":"\"mv\"", "nspin":2}
            for l in range(nType):
                dataDict[f"starting_magnetization({l+1})"] = magMoments[list(magMoments.keys())[l]]

        lines = ["&SYSTEM\n"]
        for k, v in dataDict.items():
            lines.append(f"   {k:<26}{'='}  {v:<4}\n")
        systemCard = lines + ["/\n\n"]
        return systemCard + electronsCard + KPTS + AtomicSpeciesLines + ["\n"] + cellPara + ["\n\n"] + AtomicPositions


    def MagCellLines(self):
        if self.Type.lower().endswith("afm"):
            nTotal, nType, cellLines, Atoms = self.nTotal(), self.nType(), self.cellStr(), self.Atoms()
            AtomicSpeciesLines = []
            PosHead = len(cellLines)-nTotal - 1
            AtomicPositions = [cellLines[PosHead]]
            magAtom = input("\n\tPlease Choose Atom for AFM: ")
            for l in cellLines[:nType+2]:
                if magAtom in l:
                    ls = l.split()
                    AtomicSpeciesLines.append(f"{magAtom}1  {ls[1]}   {ls[2]}\n")
                    AtomicSpeciesLines.append(f"{magAtom}2  {ls[1]}   {ls[2]}\n")
                else:
                    AtomicSpeciesLines.append(l)
            sign = "plus"
            for ind, l in enumerate(cellLines[PosHead + 1:]):
                ls = l.split()
                if magAtom in ls[0] and ind % 2 == 0:
                    atm = f"{ls[0]}1"
                    AtomicPositions.append(f"{atm:<4} {ls[1]:>18}{ls[2]:>18}{ls[3]:>18}\n")
                elif magAtom in ls[0] and ind % 2 != 0:
                    atm = f"{ls[0]}2"
                    AtomicPositions.append(f"{atm:<4} {ls[1]:>18}{ls[2]:>18}{ls[3]:>18}\n")
                else:
                    AtomicPositions.append(f"{ls[0]:<4} {ls[1]:>18}{ls[2]:>18}{ls[3]:>18}\n")
            return AtomicSpeciesLines + ["/\n"] + AtomicPositions
        else:
            cellLines = self.cellStr()
            return cellLines


    def electrons_tag(self):
        dataDict = {"conv_thr":1.000e-08, "electron_maxstep":400, "mixing_beta":0.7,
                "startingpot":"\"atomic\"", "startingwfc":"\"atomic+random\""}
        lines = ["&ELECTRONS\n"]
        for k, v in dataDict.items():
            lines.append(f"   {k:<18}{'='}  {v:<4}\n")
        return lines + ["/\n\n"]


    def ions_tag(self):
        dataDict = {"ion_dynamics":"\"bfgs\""}
        lines = ["&IONS\n"]
        for k, v in dataDict.items():
            lines.append(f"   {k:<18}{'='}  {v:<4}\n")
        return lines + ["/\n\n"]


    def cell_tag(self, press=0):
        dataDict = {"cell_dynamics":"\"bfgs\"", "press":press, "press_conv_thr":0.5}
        lines = ["&CELL\n"]
        for k, v in dataDict.items():
            lines.append(f"   {k:<18}{'='}  {v:<4}\n")
        return lines + ["/\n\n"]


    def qeToVasp(self):
        mainLines = [f"{self.pref()}\n", f"1.0\n"]
        coordType = "Cartesian"
        if self.coordType() == "angstrom":
            pass
        else:
            coordType = "Direct"
        Coords = []
        for l in self.cell():
            l = [f"{a:.15f}" for a in l]
            mainLines.append(f"{       l[0]:>22}{l[1]:>22}{l[2]:>22}\n")
        mainLines += [f'   {"  ".join(self.Atoms())}\n', f'   {"  ".join([str(a) for a in self.nAtoms()])}\n']
        if self.cType.lower().startswith("c") and not self.coordType().startswith("a"):
            coordType = f"Cartesian"
            Coords = self.toCart()
        elif self.cType.lower().startswith("f") and not self.coordType().startswith("c"):
            coordType = f"Direct"
            Coords = self.toFrac()
        else:
            coordType = f"{coordType}"
            Coords = self.coords()
        mainLines.append(f"{coordType}\n")
        for l in Coords:
            l = [f"{a:.12f}" for a in l]
            mainLines.append(f"   {l[0]:>18}{l[1]:>18}{l[2]:>18}\n")
        return mainLines


    def get_pseudo(self):
        pseudos = []
        for e in self.Atoms():
            pseudos.append(qePseudo(e))
        return pseudos


    def dos_tags(self):
        dataDict = {"prefix":f"\"{self.pref()}\"", "outdir":"\"./tmp/\"", "fildos":f"\"{self.pref()}.dos\"", "emin":-10.0, "emax":15.0}
        lines = ["&DOS\n"]
        for k, v in dataDict.items():
            lines.append(f"   {k:<18}{'='}  {v:<4}\n")
        return lines + ["/"]

    def bands_tags(self):
        pass

    def elph_tags(self):
        pass

    def writeInputFile(self, press=0, oFile=f"NA"):
        if oFile == "NA" or oFile == "DEFAULT":
            oFile = f"{self.File}.in"
        if self.Type == "vc-relax":
            lines = self.control_tag() + self.system_tag() + self.electrons_tag() + self.ions_tag()\
                    + self.cell_tag(press) + self.cellStr()
        elif "fm" in self.Type:
            lines = self.control_tag() + self.magCellStr() #self.MagSystem_tag() + self.MagCellLines()
        elif self.Type == "elph":
            lines = self.elph()
        elif "vasp" in self.Type.lower():
            lines = self.qeToVasp()
            if oFile == "NA":
                oFile = f"{self.File}.vasp"
            else:
                pass

        elif self.Type == "dos":
            lines = self.dos_tags()
        else:
            lines = self.control_tag() + self.system_tag() + self.electrons_tag() + self.cellStr()

        with open(oFile, 'w') as f0:
            for l in lines:
                f0.write(l)
        print(f"\n{'='*80}")
        print(f"\t{oFile} pwscf input file generated!")
        print(f"\n{'='*80}")



    def relax(self):
        pass

    def scf(self):
        pass

    def nscf(self):
        pass


    def elph(self):
        try:
            kpts = self.genKpts("scf", 0.08)
        except:
            kpts = (1,1,1)
        dataDict = {"tr2_ph":"1.0d-14", "prefix":"\'Tc_calculation\'", "fildvscf":"\'Tc_calculationdv\'", "outdir":"\'./tmp/\'",\
                "fildyn":"\'Tc_calculation.dyn\'", "electron_phonon":"\'interpolated\'",  "el_ph_sigma":0.005,\
                "el_ph_nsigma":10, "trans":".true.", "ldisp":".true.", "nq1":kpts[0], "nq2":kpts[1], "nq3":kpts[2]}
        lines = ["&inputph\n"]
        for k, v in dataDict.items():
            lines.append(f"   {k:<18}{'='}  {v:<4}\n")
        return ["Electron-phonon coefficients\n"] + lines + [ "/"]

    def dos(self):
        pass

    def bands(self):
        pass

    def lambdax(self):
        pass

#Getting QE Pseudopotentials

def qePseudo(element):
    """ function to get appropriate qe-pps from a specied directory location """
    #src = "/mnt/lustre3/toBeSync/Abdul/scripts/qe_pseudos/"
    src = ppPath


    files = next(os.walk(src))[2]
    #files = [f for f in os.listdir(src) if os.path.isfile(f)]
    def getInfo(f):
        wav, chg = "NA", "NA"
        with open(f, 'r') as f0:
            lines = [next(f0) for x in range(50)]
            #lines = f0.readlines()
            for l in lines:
                if "Suggested minimum cutoff for wavefunctions" in l:
                    wav = l.split()[-2]
                    chg = lines[lines.index(l)+1].split()[-2]
        return  wav, chg


    current = []
    for f in files:
        if element == f.split(".")[0]:
            current.append(f)

    print(f"\n\n\tAvailable PP for {element} (SELECT):\n")
    for i, f in enumerate(current):
        wav, chg = getInfo(src + f)
        print(f"\t({i+1}). {f} INFO -> cutoff WAV: {wav} | cutoff CHG: {chg}")
    if len(current) == 0:
        print(f"\tNo valid element specified!")
    print("\n\n")
    e = input('\tInput: ')
    if e.isdigit():
        if int(e) in range(1, len(current)+1):
            pseudo = current[int(e)-1]

            os.system(f"cp {src}{pseudo} .")

            print(f"\n\tQE Pseudopotential \"{pseudo}\" has been copied in the PWD!\n")
            return pseudo
        else:
            print("\n\tInvalid choice! (SET MANUALLY!)\n")
            return "set_pseudoPot_manually.UPF"
    else:
        print("\n\tInvalid choice! (SET MANUALLY!)\n")
        return "set_pseudoPot_manually.UPF"



cell = vaspToQE(iFile, cal, cType)
cell.writeInputFile(press, oFile)
