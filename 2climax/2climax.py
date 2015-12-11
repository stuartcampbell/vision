#! /usr/bin/env python

import sys, os, re


def unix2dos(file_aclimax):
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# replace \n with \r\n in the *.aclimax file
#
    f = open(file_aclimax, 'r')
    full_text = f.read()
    f.close()
    f = open(file_aclimax, 'w')
    f.write(full_text.replace("\n", "\r\n"))
    f.close()
#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


def get_elements(file_elements, elem):
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Input: file name of the database, list of symbols/types (Z) to look up
# Output: lists of type/symbol, mass, and sigma
#
    f = open(file_elements, 'r')
    
    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    for line in f:
        aux = line.split()
        col1.append(int(aux[0]))        # column 1: type (Z number)
        col2.append(aux[1])             # column 2: symbol
        col3.append(float(aux[2]))      # column 3: mass
        col4.append(float(aux[3]))      # column 4: sigma-total
        col5.append(float(aux[4]))      # column 5: sigma-inc
    f.close()
    
    mass = []
    sigma = []

    if str(elem[0]).isdigit():          # input list is "type"
        table = dict(zip(col1,range(len(col1)))) 
        symbol = []
        for tp in elem:
            symbol.append(col2[table[tp]])
            mass.append(col3[table[tp]])
            sigma.append(col4[table[tp]])
        return symbol, mass, sigma
    else:                               # input list is symbol
        table = dict(zip(col2,range(len(col2))))
        type = []
        for syb in elem:
            type.append(col1[table[syb]])
            mass.append(col3[table[syb]])
            sigma.append(col4[table[syb]])
        return type, mass, sigma
#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


def get_poscar(file_poscar):
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Input: VASP POSCAR of the unit cell
# Output: lists of symbols and coordinates
#
    f = open(file_poscar, 'r')
    
    f.readline()
    f.readline()     # discard the first two lines
    
    lattice = [[0 for _ in range(3)] for _ in range(3)]  # unit cell lattice
    line = f.readline()
    lattice[0] = map(float, line.split()[:3])
    line = f.readline()
    lattice[1] = map(float, line.split()[:3])
    line = f.readline()
    lattice[2] = map(float, line.split()[:3])
    
    line = f.readline()
    elem = line.split()              # elements in unit cell
    
    line = f.readline()
    nelem = map(int, line.split())   # number of atoms for each element
    f.readline()
    
    na = sum(nelem)
    xyz = [[0 for x in range(3)] for x in range(na)] 
    for i in range(na):
    	line = f.readline()
    	xyz[i] = map(float, line.split()[:3])     # fractional coordinates
    f.close()

    symbol = []
    for i in range(len(elem)):
    	for j in range(nelem[i]):
    		symbol.append(elem[i])        # symbol list

    for i in range(na):
        sx = xyz[i][0]*lattice[0][0] + \
             xyz[i][1]*lattice[1][0] + \
             xyz[i][2]*lattice[2][0]
        sy = xyz[i][0]*lattice[0][1] + \
             xyz[i][1]*lattice[1][1] + \
             xyz[i][2]*lattice[2][1]
        sz = xyz[i][0]*lattice[0][2] + \
             xyz[i][1]*lattice[1][2] + \
             xyz[i][2]*lattice[2][2]
        xyz[i] = [sx, sy, sz]                 # absolute xyz

    return symbol, xyz    
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


def get_phonopy(file_phonopy):
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Input: mesh.yaml from phonopy (VASP force constants) 
# Output: phonon frequencies and modes
#
    f = open(file_phonopy, 'r')
    
    line = f.readline()
    qmesh = map(int, re.findall(r"[-+]?\d*\.\d+|\d+", line))
    qtot = 1
    for i in qmesh:
        qtot = qtot*i             # number of total q points in BZ
    
    line = f.readline()                               
    nq = int(line.split()[1])     # number of irreducible q points

    line = f.readline()
    na = int(line.split()[1])     # number of atoms
    nb = na*3                     # number of bands at each q point
    f.readline()
    
    wq = [0 for i in range(nq)]   # weight of each irreducible q point
    freq = [[0 for i in range(nb)] for i in range(nq)]
    mode = [[[[[0 for i in range(2)] for i in range(3)]
           for i in range(na)] for i in range(nb)] for i in range(nq)]
    dxyz = [[[[0 for i in range(3)]
           for i in range(na)] for i in range(nb)] for i in range(nq)]
    for iq in range(nq):
        f.readline()
        line = f.readline()
        wq[iq] = float(line.split()[1])
        f.readline()
        for ib in range(nb):
            f.readline()
            line = f.readline()
            freq[iq][ib] = float(line.split()[1])
            freq[iq][ib] *= 33.3564  # convert from THz to cm-1
            f.readline()
            for ia in range(na):
                f.readline()
                for i in range(3):
                    line=f.readline()
                    mode[iq][ib][ia][i] = map(float, 
                        re.findall(r"[-+]?\d*\.\d+|\d+", line))
                xr = mode[iq][ib][ia][0][0]  # normalized sqrt(mass)*u
                xi = mode[iq][ib][ia][0][1]  
                yr = mode[iq][ib][ia][1][0]
                yi = mode[iq][ib][ia][1][1]
                zr = mode[iq][ib][ia][2][0]
                zi = mode[iq][ib][ia][2][1]
                dxyz[iq][ib][ia][0] = cmp(xr,0)* \
                    ((xr**2+xi**2)*wq[iq]/qtot)**(0.5)
                dxyz[iq][ib][ia][1] = cmp(yr,0)* \
                    ((yr**2+yi**2)*wq[iq]/qtot)**(0.5)
                dxyz[iq][ib][ia][2] = cmp(zr,0)* \
                    ((zr**2+zi**2)*wq[iq]/qtot)**(0.5)
        f.readline()
    
#    if(sum(wq)!=qtot): print "error in yaml file, sum(wq) != qtot"
    
    f.close()

    return nq, nb, freq, dxyz
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


def get_phout(file_phout):
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Input: *.ph.out from Quantum Espresso (ph.x)
# Output: unit cell information
#
    f = open(file_phout, 'r')
    celldm = [0 for _ in range(6)]
    lattice = [[0 for _ in range(3)] for _ in range(3)]
    xyz = []
    symbol = []
    finishing = False
    for line in f:
        aux = line.split()
        if len(aux) == 0 and finishing:
            break
        if len(aux) == 0:
            continue
        if aux[0] == r"celldm(1)=":
            celldm[0] = float(aux[1])
            celldm[1] = float(aux[3])
            celldm[2] = float(aux[5])
        if aux[0] == r"celldm(4)=":
            celldm[3] = float(aux[1])
            celldm[4] = float(aux[3])
            celldm[5] = float(aux[5])
        if aux[0] == r"a(1)":
            lattice[0][0] = float(aux[3])
            lattice[0][1] = float(aux[4])
            lattice[0][2] = float(aux[5])
        if aux[0] == r"a(2)":
            lattice[1][0] = float(aux[3])
            lattice[1][1] = float(aux[4])
            lattice[1][2] = float(aux[5])
        if aux[0] == r"a(3)":
            lattice[2][0] = float(aux[3])
            lattice[2][1] = float(aux[4])
            lattice[2][2] = float(aux[5])
        if len(aux) > 4 and re.match(r"tau.*", aux[2]) != None:
            xyz.append(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", line))[2:5])
            symbol.append(aux[1])
            finishing = True

    if celldm[1] == 0.0:
       celldm[1] = celldm[0]
    else:
       celldm[1] = celldm[1]*celldm[0]
    if celldm[2] == 0.0:
       celldm[2] = celldm[2]*celldm[0]
    for i in range(len(xyz)):    # only orthorhombic cell for now
        x = [0 for _ in range(3)]
        for j in range(3):
            for k in range(3):
                x[k] = x[k] + celldm[k] * lattice[k][j] * xyz[i][k] 
        for k in range(3):
            xyz[i][k] = x[k] 

    na = len(symbol)
    nb = 3*na

    return symbol, xyz, nb, na
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


def get_mesh(file_mesh):
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Input: mesh_k file generated by Quantum Espresso (kpoints.x) 
# Output: irreducible q points and degeneration (weight)
#
    f = open(file_mesh, 'r')
    wq = []
    line = f.readline()
    nq = int(line.split()[0])
    for i in range(nq):
        line = f.readline()
        wq.append(float(line.split()[4]))
    f.close()
    
    return nq, wq
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


def get_matdyn(file_matdyn, nq, nb, na, wq):
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Input: matdyn.modes from Quantum Espresso 
# Output: phonon frequencies and modes
#
    f = open(file_matdyn, 'r')
    qtot = sum(wq)
    freq = [[0 for i in range(nb)] for i in range(nq)]
    mode = [[[[0 for i in range(6)]
           for i in range(na)] for i in range(nb)] for i in range(nq)]
    dxyz = [[[[0 for i in range(3)]
           for i in range(na)] for i in range(nb)] for i in range(nq)]
    for iq in range(nq):
        for _ in range(4):
            f.readline()     # discard the first four lines
        for ib in range(nb):
            line = f.readline()
            freq[iq][ib] = map(float, 
                         re.findall(r"[-+]?\d*\.\d+|\d+", line))[2]
            for ia in range(na):
                line = f.readline()
                mode[iq][ib][ia] = map(float, 
                        re.findall(r"[-+]?\d*\.\d+|\d+", line))
                xr = mode[iq][ib][ia][0]  # normalized sqrt(mass)*u
                xi = mode[iq][ib][ia][1] 
                yr = mode[iq][ib][ia][2]
                yi = mode[iq][ib][ia][3]
                zr = mode[iq][ib][ia][4]
                zi = mode[iq][ib][ia][5]
                dxyz[iq][ib][ia][0] = cmp(xr,0)* \
                    ((xr**2+xi**2)*wq[iq]/qtot)**(0.5)
                dxyz[iq][ib][ia][1] = cmp(yr,0)* \
                    ((yr**2+yi**2)*wq[iq]/qtot)**(0.5)
                dxyz[iq][ib][ia][2] = cmp(zr,0)* \
                    ((zr**2+zi**2)*wq[iq]/qtot)**(0.5)
        f.readline()

    return freq, dxyz
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


def get_gaussian(file_gaussian, file_elements):
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Input: Gaussian output file 
# Output: atomic coordinates and phonon information
# Note: copied from Timmy's code, minor changes for consistency in notations
#
    f = open(file_gaussian, 'r')
    full_text = f.read()
    f.close() 

# Selects the text between 2 markers for positions (last appearance)
#
    aux = re.compile('(?=Standard orientation)(.*?)(?=Rotational constants)',
          re.DOTALL).findall(full_text)  
    chopped = aux[len(aux)-1].split('\n')

# Reads the atomic positions
#
    atoms = []
    type = []
    xyz = []
    for i in range(5,len(chopped)-2):
        atoms.append(chopped[i].split())
        na = len(atoms)
    for i in range(na):
        atom = atoms[i]
        type.append(int(atom[1]))
        xyz.append(map(float, atom[3:6]))
    symbol, mass, sigma = get_elements(file_elements, type)

# Selects the text between 2 markers for positions Frequencies
#
    aux = re.compile('(?=Frequencies)(.*?)Thermochemistry', 
          re.DOTALL).search(full_text)
    chopped = aux.group(1)
    aux = re.compile('(?=Frequencies)(.*?)-------------------', 
          re.DOTALL).findall(chopped)
    A = aux[0].split('\n')

# Reads the frequencies and stores them in multidimensional array
#
    freq = []
    dxyz = []
    for i in range(len(A)):
        aux = A[i].split()
        if re.search('Frequencies', A[i]):
            disp1 = []
            disp2 = []
            disp3 = []
            reading = True
            for j in range(2,5):
                freq.append(float(aux[j]))
        if len(aux) > 6:
            if (aux[0].isdigit()):
                disp1.append([float(aux[2]),float(aux[3]),float(aux[4])])
                disp2.append([float(aux[5]),float(aux[6]),float(aux[7])])
                disp3.append([float(aux[8]),float(aux[9]),float(aux[10])])
        if (len(aux) < 6) and reading:
            dxyz.append(disp1)
            dxyz.append(disp2)
            dxyz.append(disp3)
            reading = False

# Produces the mass weighted normal modes eigenvectors (normalised)
#
  
    for i in range(len(dxyz)):
    # iteration through frequency
        aux = 0.0
        for j in range(len(dxyz[0])):
        # Iteration through atoms
            for k in range(0,len(dxyz[0][0])):
                aux = aux + (dxyz[i][j][k]**2) * mass[j]
        for j in range(len(dxyz[0])):
        # Iteration through atoms
            for k in range(len(dxyz[0][0])):
                dxyz[i][j][k] = dxyz[i][j][k] * (mass[j]/aux)**(0.5)

    nq = 1            # single molecule, gamma point only
    nb = len(freq)
    freq = [freq]
    dxyz = [dxyz]

    return type, symbol, mass, sigma, xyz, nq, nb, freq, dxyz  
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


def get_castep_line(file_castep, file_elements):
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Input: *.phonon from CASTEP 
# Output: atomic coordinates and phonon information
# Note: file is read line-by-line with all information retained
#
    f = open(file_castep, 'r')
    
    f.readline()
    line = f.readline()
    na = int(re.findall(r"\d+", line)[0])
    line = f.readline()
    nb = int(re.findall(r"\d+", line)[0])
    line = f.readline()
    nq = int(re.findall(r"\d+", line)[0])

    for _ in range(4):
        f.readline()
    lattice = [[0 for _ in range(3)] for _ in range(3)]
    for i in range(3):
        line = f.readline()
        lattice[i] = map(float, line.split())
    f.readline()
    xyz = [[0 for _ in range(3)] for _ in range(na)]
    symbol = [0 for _ in range(na)]
    type = []
    mass = []
    sigma = []
    for i in range(na):
        line = f.readline()
        xyz[i] = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", line)[1:4])
        symbol[i] = re.split(r"\s+", line)[5]    
    type, mass, sigma = get_elements(file_elements, symbol)

    for i in range(len(xyz)):
        x = [0 for _ in range(3)]
        for j in range(3):
            for k in range(3):
                x[k] = x[k] + lattice[k][j] * xyz[i][k] 
        for k in range(3):
            xyz[i][k] = x[k] 
    f.readline()
    wq = [0 for i in range(nq)]     # weight of each irreducible q point
    freq = [[0 for i in range(nb)] for i in range(nq)]
    dxyz = [[[[0 for i in range(3)]
           for i in range(na)] for i in range(nb)] for i in range(nq)]
    for iq in range(nq):
        line = f.readline()
        wq[iq] = float(line.split()[5])
        for ib in range(nb):
            line = f.readline()
            freq[iq][ib] = float(line.split()[1])
        f.readline()
        f.readline()
        for ib in range(nb):
            for ia in range(na):
                line = f.readline()
                aux = map(float, line.split()[2:8])
                dxyz[iq][ib][ia][0] = cmp(aux[0],0)* \
                    ((aux[0]**2+aux[1]**2)*wq[iq])**(0.5)
                dxyz[iq][ib][ia][1] = cmp(aux[2],0)* \
                    ((aux[2]**2+aux[3]**2)*wq[iq])**(0.5)
                dxyz[iq][ib][ia][2] = cmp(aux[4],0)* \
                    ((aux[4]**2+aux[5]**2)*wq[iq])**(0.5)
    f.close()

    return type, symbol, mass, sigma, xyz, nq, nb, freq, dxyz
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

def get_castep_all(file_castep,file_elements):
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Input: *.phonon from CASTEP 
# Output: atomic coordinates and phonon information
# Note: copied from Timmy's code, file read all at once
#

#  Read the file in memory (may be a problem with CASTEP)
#
    f = open(file_castep,'r')
    full_text = f.read()
    f.close() 
    full_text = full_text + 'q-pt='
# Last line added to improve the search of eigenvectors

# Selects the text between 2 markers for positions
# Reads Unit Cell Vectors
#
    aux = re.compile('Unit cell vectors(.*?)Fractional Co-ordinates',
                    re.DOTALL).search(full_text)
    chopped = aux.group(1).split('\n')
    lattice = []
    for i in range(len(chopped)):
        extra = chopped[i].split()
        if len(extra) == 3 :
            lattice.append([float(extra[0]),float(extra[1]),float(extra[2])])

# Selects the text between 2 markers for positions
# Fractional Coordinates and atomic detail
#
    aux = re.compile('Fractional Co-ordinates(.*?)END header', 
                     re.DOTALL).search(full_text)
    chopped = aux.group(1).split('\n')
    xyz = []
    symbol = []
    mass = []
    type = []
    sigma = []
    for i in range(0,len(chopped)):
        extra = chopped[i].split()
        if len(extra) == 6 :
            xyz.append([float(extra[1]),float(extra[2]),float(extra[3])])
            symbol.append(extra[4])
            mass.append(float(extra[5]))
    type, mass, sigma = get_elements(file_elements, symbol)
    for i in range(len(xyz)):
        x = [0 for _ in range(3)]
        for j in range(3):
            for k in range(3):
                x[k] = x[k] + lattice[k][j] * xyz[i][k] 
        for k in range(3):
            xyz[i][k] = x[k] 

# Selects the text between 2 markers for positions
# stores in variable aux
#
    aux = re.compile('(?=q-pt=)(.*?)Phonon Eigenvectors', 
          re.DOTALL).findall(full_text)
    freq = []
    weights = []
    for i in range(len(aux)):
        extra = aux[i].split('\n')
        Peso = extra[0].split()
        for j in range(len(extra)):
            InQpt = extra[j].split()
            if len(InQpt) > 0 :
                if InQpt[0].isdigit():
                    freq.append(float(InQpt[1]))
                    weights.append(float(Peso[5]))

# Selects the text between 2 markers for positions
# stores in variable aux
#
    aux = re.compile('(?=Mode Ion)(.*?)(?=q-pt=)', 
          re.DOTALL).findall(full_text)
    FrqNum = 0
    dxyz = []

    for i in range(len(aux)):
        extra = aux[i].split('\n')
        disp1 = []
        for j in range(len(extra)-1):
            InQpt = extra[j].split()
            InQpt2 = extra[j+1].split()
            if len(InQpt) > 7 :
                MassAtom = (mass[int(InQpt[1])-1])**(0.5)
      
                DX = (float(InQpt[2])**2+float(InQpt[3])**2)**(0.5)
                DY = (float(InQpt[4])**2+float(InQpt[5])**2)**(0.5)
                DZ = (float(InQpt[6])**2+float(InQpt[7])**2)**(0.5)
        
                DX = DX * cmp(float(InQpt[2]), 0)
                DY = DY * cmp(float(InQpt[4]), 0)
                DZ = DZ * cmp(float(InQpt[6]), 0)
        
                disp1.insert(int(InQpt[1])-1,[DX,DY,DZ])
    
                if len(InQpt2) == 0:
                    InQpt2.append('Carajo')
                if InQpt[0] != InQpt2[0]:
                    dxyz.append(disp1)
                    disp1 = []

    for i in range(len(dxyz)):
        acum = 0
        for j in range(0,len(dxyz[0])):
            for k in range(0,3):
                acum += dxyz[i][j][k]**2
        for j in range(0,len(dxyz[0])):
            for k in range(0,3):
                dxyz[i][j][k] = dxyz[i][j][k] * (weights[i]/acum)**(0.5)
  
    nq = 1           # get_castep_all does not differentiate q points and
    nb = len(freq)   # branches. Equivalent to gamma point only.
    freq = [freq]    # Use get_castep_line to extract the q-mesh info.
    dxyz = [dxyz] 
    return type, symbol, mass, sigma, xyz, nq, nb, freq, dxyz
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


def write_aclimax(file_aclimax, type, symbol, mass, sigma, xyz, 
                  nq, nb, freq, dxyz):
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Write aClimax file
#
    f = open(file_aclimax, 'w')
    
    print >> f, "BEGIN ATOMS"
    for i in range(len(type)):
        print >> f, '%5d' % (i+1), '%2s' % symbol[i],'%4d' % type[i],\
                  '%12.6f' % mass[i],'%12.6f' % xyz[i][0],\
                  '%12.6f' % xyz[i][1],'%12.6f' % xyz[i][2],\
                  '%12.6f' % sigma[i]
    
    print >> f, "END ATOMS"     
    print >> f
    print >> f
    print >> f, "BEGIN FREQUENCIES"
    for iq in range(nq):
        for ib in range(nb):
            print >> f, "!MODE# FREQUENCY"
            print >> f, '%5d' % (iq*nb+ib+1),' %12.4f' % freq[iq][ib]
            print >> f
            for ia in range(len(type)):
                print >> f, '%4d'% (ia+1), '%2s' % symbol[ia], \
                             '%4d'% type[ia], \
                             '%12.6f' % dxyz[iq][ib][ia][0], \
                             '%12.6f' % dxyz[iq][ib][ia][1], \
                             '%12.6f' % dxyz[iq][ib][ia][2]
    print >> f, "END FREQUENCIES"        
    
    f.close()

    unix2dos(file_aclimax)

#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


#=============================================================================
# Main program starts here
#

file_elements = os.path.dirname(os.path.realpath(__file__)) + "/elements.txt"
file_aclimax = "out.aclimax" 

if sys.argv[1] in ['-v', '-V', '-vasp', '-VASP']:  # from VASP
    if len(sys.argv) < 4: 
        print "VASP: two input files (POSCAR and mesh.yaml)"
        quit()
    if len(sys.argv) == 5: 
        file_aclimax = sys.argv[4]
    if len(sys.argv) > 5: 
        print "Too many arguments"
        quit()
    symbol, xyz = get_poscar(sys.argv[2])
    type, mass, sigma = get_elements(file_elements, symbol)
    nq, nb, freq, dxyz = get_phonopy(sys.argv[3])
    write_aclimax(file_aclimax, type, symbol, mass, sigma, xyz, 
                  nq, nb, freq, dxyz)

elif sys.argv[1] in ['-q', '-Q', '-qe', '-QE']:   # from Quantum Espreso
    if len(sys.argv) < 5: 
        print "QE: three input files (*.ph.out, mesh_k, and matdyn.modes)"
        quit()
    if len(sys.argv) == 6: 
        file_aclimax = sys.argv[5]
    if len(sys.argv) > 6: 
        print "Too many arguments"
        quit()
    symbol, xyz, nb, na = get_phout(sys.argv[2])
    type, mass, sigma = get_elements(file_elements, symbol)
    nq, wq = get_mesh(sys.argv[3])
    freq, dxyz = get_matdyn(sys.argv[4], nq, nb, na, wq)
    write_aclimax(file_aclimax, type, symbol, mass, sigma, xyz, 
                  nq, nb, freq, dxyz)

elif sys.argv[1] in ['-c', '-C', '-castep', '-CASTEP']: # from CASTEP
    if len(sys.argv) < 3: 
        print "CASTEP: *.phonon file required"
        quit()
    if len(sys.argv) == 4: 
        file_aclimax = sys.argv[3]
    if len(sys.argv) > 4: 
        print "Too many arguments"
        quit()
    type, symbol, mass, sigma, xyz, nq, nb, freq, dxyz = \
                 get_castep_line(sys.argv[2], file_elements)
    write_aclimax(file_aclimax, type, symbol, mass, sigma, xyz, 
                  nq, nb, freq, dxyz)

elif sys.argv[1] in ['-g', '-G', '-gaussian', '-GAUSSIAN']: # from Gaussian
    if len(sys.argv) < 3: 
        print "Gaussian: *.out file required"
        quit()
    if len(sys.argv) == 4: 
        file_aclimax = sys.argv[3]
    if len(sys.argv) > 4: 
        print "Too many arguments"
        quit()
    type, symbol, mass, sigma, xyz, nq, nb, freq, dxyz = \
                 get_gaussian(sys.argv[2], file_elements)
    write_aclimax(file_aclimax, type, symbol, mass, sigma, xyz, 
                  nq, nb, freq, dxyz)

else:
    print "Command option not specified/recognized"
    quit()

