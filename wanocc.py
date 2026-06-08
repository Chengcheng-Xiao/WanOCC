#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function
import argparse
import xml.etree.ElementTree as ET
import numpy as np


# Command line praser
#----------------------------
parser = argparse.ArgumentParser(description='A script to calculate Wannier Occupation.')
parser.add_argument('-seed_name', action="store", default=str('wannier90'), dest="seed_name",
                    help='seed_name')
# parser.add_argument('-R', action="store", default=[0,0,0], dest="R_latt",
#                     help='Which R_latt to use? Default: 0,0,0')
parser.add_argument('-v', action="store_true", default=False, dest="verbose",
                    help='Verbose output, print the diagonal part in the home cell. Default: True')
parser.add_argument('-spin', action="store", default="unpolarized", dest="spin",
                    help='Which spin channel? (Default=unpolarized, available "up"/"down"/"unpolarized")')
parser.add_argument('-bnd_exc', action="store", default=None, dest="bnd_exclude",
                    help='Which bands to exclude? (Default=None)')
parser.add_argument('-dis', action="store_true", default=False, dest="disentangle",
                    help='Do we want to use disentanglement procedure?')
parser.add_argument('-dft', action="store", default="qe", dest="dft",
                    help='What DFT code? (Default=qe, available "qe"/"vasp")')

prm = parser.parse_args()

#%% parseIntSet
def parseIntSet(nputstr=""):
    selection = []
    invalid = []
    # tokens are comma seperated values
    tokens = [x.strip() for x in nputstr.split(',')]
    for i in tokens:
        if len(i) > 0:
            if i[:1] == "<":
                i = "1-%s"%(i[1:])
        try:
            # typically tokens are plain old integers
            selection.append(int(i))
        except:
            # if not, then it might be a range
            try:
                token = [int(k.strip()) for k in i.split('-')]
                if len(token) > 1:
                    token.sort()
                    # we have items seperated by a dash
                    # try to build a valid range
                    first = token[0]
                    last = token[len(token)-1]
                    for x in range(first, last+1):
                        selection.append(x)
            except:
                # not an int and not a range...
                invalid.append(i)
    # Report invalid tokens before returning valid selection
    if len(invalid) > 0:
        print("Invalid set: " + str(invalid))
    return np.array(selection, dtype=int)

def extract_occupations_array_vasp(vasprun_file):
    tree = ET.parse(vasprun_file)
    root = tree.getroot()

    # Find eigenvalues section
    eigenvalues = root.find(".//eigenvalues")
    if eigenvalues is None:
        raise ValueError("No <eigenvalues> section found.")

    # Top-level <set>
    top_set = eigenvalues.find(".//set")
    spin_sets = top_set.findall("set")

    all_occ = []

    for spin_set in spin_sets:
        spin_occ = []

        # k-points inside this spin
        kpoint_sets = spin_set.findall("set")

        for k_set in kpoint_sets:
            band_occ = []

            # each <r> contains: eigenvalue occupation
            for r in k_set.findall("r"):
                vals = r.text.split()

                occ = float(vals[1])   # second column = occupation
                band_occ.append(occ)

            spin_occ.append(band_occ)

        all_occ.append(spin_occ)

    # Convert to numpy array
    occ_array = np.array(all_occ)

    return occ_array


def read_wannier90_kpoints(filename="wannier90.win"):
    """
    Read k-points from the block:
    
        begin kpoints
        ...
        end kpoints

    in a wannier90.win file.

    Returns
    -------
    kpoints : list of tuple
        List of (kx, ky, kz)
    """

    kpoints = []
    inside_block = False

    with open(filename, "r") as f:
        for line in f:
            line_stripped = line.strip()
            # Skip empty lines and comments
            if not line_stripped or line_stripped.startswith(("#", "!")):
                continue
            lower = line_stripped.lower()

            if lower == "begin kpoints":
                inside_block = True
                continue
            if lower == "end kpoints":
                inside_block = False
                break

            if inside_block:
                parts = line_stripped.split()
                if len(parts) >= 3:
                    kx, ky, kz = map(float, parts[:3])
                    kpoints.append((kx, ky, kz))
    return kpoints

def read_wannier90_hr(filename="wannier90_hr.dat"):
    """
    Read Wannier90 *_hr.dat file.

    Returns
    -------
    nrpts : int
        Number of R-points.

    degeneracies : ndarray of shape (nrpts,)
        Degeneracy for each R-point.

    R_vectors : ndarray of shape (nrpts, 3)
        Integer lattice vectors R = (Rx, Ry, Rz).
    """

    with open(filename, "r") as f:
        lines = f.readlines()

    # Header
    comment = lines[0].strip()

    num_wann = int(lines[1])
    nrpts = int(lines[2])

    # ------------------------------------------------------------
    # Read degeneracies
    # Degeneracies are written 15 integers per line
    # ------------------------------------------------------------
    degeneracies = []

    idx = 3
    while len(degeneracies) < nrpts:
        degeneracies.extend(map(int, lines[idx].split()))
        idx += 1

    degeneracies = np.array(degeneracies[:nrpts], dtype=int)

    # ------------------------------------------------------------
    # Read Hamiltonian matrix elements
    #
    # Each line format:
    #   Rx Ry Rz   m   n   Re(H)   Im(H)
    #
    # There are nrpts * num_wann^2 lines
    # ------------------------------------------------------------
    R_vectors = []
    total_lines = nrpts * num_wann * num_wann
    current_R = None

    for i in range(total_lines):
        parts = lines[idx + i].split()
        Rx, Ry, Rz = map(int, parts[:3])
        R = (Rx, Ry, Rz)

        # Each R appears num_wann^2 times;
        # keep only unique consecutive R's
        if R != current_R:
            R_vectors.append(R)
            current_R = R

    R_vectors = np.array(R_vectors, dtype=int)

    return nrpts, degeneracies, R_vectors

#%% input
if prm.bnd_exclude != None:
    bnd_exclude = parseIntSet(prm.bnd_exclude)
    tot_bnd_exc = int(bnd_exclude.shape[0])
else:
    bnd_exclude = []
    tot_bnd_exc = 0
seed_name = prm.seed_name

#-------------------------------------------------------------------------------
if prm.dft == "qe":
    # get me starting index of ks_energies in root
    tree = ET.parse(seed_name+'.xml')
    root = tree.getroot()

    starting_index = [i.tag for i in root[3][-1]].index('ks_energies')

    # get me occupation of KS states.
    occupation_mat = []
    for i in root[3][-1][starting_index:]:
        occupation_mat.append(i[3].text.split())

elif prm.dft == "vasp":
    occ = extract_occupations_array_vasp("vasprun.xml")
    # occ Shape = (nspin, nkpts, nbands)
    if occ.shape[0] == 1:
        # unpolarized case, duplicate the spin channel
        occupation_mat = occ[0]
    elif occ.shape[0] == 2:
        # spin-polarized case, concatenate the two spin channels
        occupation_mat = np.column_stack([
            occ[0, :, :],
            occ[1, :, :]
        ])
    else:
        raise ValueError("Unexpected number of spin channels in vasprun.xml: {}".format(occ.shape[0]))

else:
    raise ValueError("Unsupported DFT code: {}. Available options: qe/vasp".format(prm.dft))

#-------------------------------------------------------------------------------
# get me k point list from seedname.win
k_points = read_wannier90_kpoints(seed_name+'.win')

# get me R_vectors and degeneracies from seedname_hr.dat
nrpts, degen, R_vectors = read_wannier90_hr(seed_name+'_hr.dat')


#-------------------------------------------------------------------------------
if prm.spin == "up":
    occupation_mat_up=np.zeros([len(occupation_mat),int(len(occupation_mat[0])/2-tot_bnd_exc)],dtype=float)
    for i in range(len(occupation_mat)):
        i_bnd_tot = 0
        for i_bnd in range(int(len(occupation_mat[0])/2)):
            if i_bnd+1 in bnd_exclude:
                continue
            occupation_mat_up[i,i_bnd_tot]=occupation_mat[i][i_bnd]
            i_bnd_tot += 1
elif prm.spin == "down":
    occupation_mat_up=np.zeros([len(occupation_mat),int(len(occupation_mat[0])/2-tot_bnd_exc)],dtype=float)
    for i in range(len(occupation_mat)):
        i_bnd_tot = 0
        for i_bnd in range(int(len(occupation_mat[0])/2)):
            if i_bnd+1 in bnd_exclude:
                continue
            occupation_mat_up[i,i_bnd_tot]=occupation_mat[i][int(len(occupation_mat[0])/2+i_bnd)]
            i_bnd_tot += 1
elif prm.spin == "unpolarized":
    occupation_mat_up=np.zeros([len(occupation_mat),len(occupation_mat[0])-tot_bnd_exc],dtype=float)
    for i in range(len(occupation_mat)):
        i_bnd_tot = 0
        for i_bnd in range(len(occupation_mat[0])):
            if i_bnd+1 in bnd_exclude:
                continue
            occupation_mat_up[i,i_bnd_tot]=occupation_mat[i][i_bnd]
            i_bnd_tot += 1
else:
    raise ValueError("spin must be up/down/unpolarized")


if prm.disentangle:
    # get me k point list and u_matrix_opt(disentanglement)
    k_list=[]
    with open(seed_name+"_u_dis.mat") as f:
        header = f.readline()
        nkpt, num_wann, num_bnd = [int(x) for x in f.readline().split()]
        f.readline()

        u_matrix_opt = np.zeros([nkpt, num_bnd, num_wann], dtype=complex)

        for ikpt in range(nkpt):
            k_list.append(f.readline().split())
            for iwann in range(num_wann):
                for ibnd in range(num_bnd):
                    u_real, u_img = [float(x) for x in f.readline().split()]
                    u_matrix_opt[ikpt, ibnd, iwann] = complex(u_real, u_img)
            f.readline()


# get me k point list and u_matrix(MLWF)
k_list=[]
with open(seed_name+"_u.mat") as f:
    header = f.readline()
    nkpt, num_wann, num_bnd = [int(x) for x in f.readline().split()]
    f.readline()

    u_matrix = np.zeros([nkpt, num_wann, num_bnd], dtype=complex)

    for ikpt in range(nkpt):
        k_list.append(f.readline().split())
        for iwann in range(num_wann):
            for ibnd in range(num_bnd):
                u_real, u_img = [float(x) for x in f.readline().split()]

                u_matrix[ikpt, ibnd, iwann] = complex(u_real, u_img)
        f.readline()


# get me occ_mat under wannier basis
#-----------------------------------

#------------------------------------------------------------------------------
# original verison, kept for reference, but it's very slow due to the triple
# loop over k-points and R-vectors
#------------------------------------------------------------------------------
#  WF_occ = np.zeros([num_wann,num_wann,len(R_vectors)],dtype=complex)
#  counter = 0
#  print("Calculating Wannier occupation in real space...")
#  for R_latt in R_vectors:
#      print("  R = ({:d}, {:d}, {:d})".format(R_latt[0], R_latt[1], R_latt[2]))
#      occ_mat=[]
#      occ_mat_pre_f = []
#      for ikpt in range(nkpt):
#          a = np.zeros([occupation_mat_up[0].size, occupation_mat_up[0].size])
#          np.fill_diagonal(a,occupation_mat_up[ikpt])
#          if prm.disentangle:
#              # pre-rotate occ_matrix with u_matrix_opt
#              occ_mat_pre = np.dot(np.dot(u_matrix_opt[ikpt],a),np.matrix(u_matrix_opt[ikpt]).H)
#              occ_mat_pre_f.append(occ_mat_pre)
#              # calculate the k factor
#              phase_factor = np.exp(complex(0,np.dot(k_points[ikpt],R_latt)*2*np.pi))
#              # rotate occ_matrix with u_matrix
#              occ_mat.append(phase_factor*np.dot(np.dot(u_matrix[ikpt],occ_mat_pre),np.matrix(u_matrix[ikpt]).H))
#          else:
#              # calculate the k factor
#              phase_factor = np.exp(complex(0,np.dot(k_points[ikpt],R_latt)*2*np.pi))
#              # rotate occ_matrix with u_matrix
#              occ_mat.append(phase_factor*np.dot(np.dot(u_matrix[ikpt],a),np.matrix(u_matrix[ikpt]).H))
#
#      # summing up kpoints weighted by k-weight (1/nkpt for wannier90 uses
#      # uniform-kpoint grid)
#      for ikpt in range(nkpt):
#          for i in range(num_wann):
#              for j in range(num_wann):
#                  WF_occ[i,j,counter] += occ_mat[ikpt][i,j]/nkpt
#
#      counter += 1

#------------------------------------------------------------------------------
# optmized version of the above code using numpy broadcasting and vectorization
#------------------------------------------------------------------------------
print("Calculating Wannier occupation in real space...")
# Precompute phase factors
phase_factors = np.exp(1j * 2 * np.pi * np.dot(k_points, R_vectors.T))

WF_occ = np.zeros((num_wann, num_wann, nrpts), dtype=complex)

for ikpt in range(nkpt):
    if prm.disentangle:
        occ_mat_pre = (u_matrix_opt[ikpt].conj().T * occupation_mat_up[ikpt]) @ u_matrix_opt[ikpt]
        #  occ_mat_pre = np.diag(occ_mat_pre)  # Extract diagonal elements after pre-rotation

        occ = occ_mat_pre # 1D diagonal elements
        U = u_matrix[ikpt]  # Shape: (num_wann, nbands)

        # Optimized rotation: U^\dagger * diag(occ) * U
        # Using broadcasting to multiply U columns by occ
        rotated = (U.conj().T @ occ) @ U

        # Vectorized accumulation across all R vectors
        WF_occ += rotated[:, :, np.newaxis] * phase_factors[ikpt, np.newaxis, np.newaxis, :]
    else:
        occ_mat_pre = occupation_mat_up[ikpt]

        occ = occ_mat_pre # 1D diagonal elements
        U = u_matrix[ikpt]  # Shape: (num_wann, nbands)

        # Optimized rotation: U^\dagger * diag(occ) * U
        # Using broadcasting to multiply U columns by occ
        rotated = (U.conj().T * occ) @ U

        # Vectorized accumulation across all R vectors
        WF_occ += rotated[:, :, np.newaxis] * phase_factors[ikpt, np.newaxis, np.newaxis, :]

# Average over k-points,assuming we have uniform k-point grid, so each k-point
# has equal weight of 1/nkpt
WF_occ /= nkpt


#-------------------------------------------------------------------------------
from datetime import datetime
now = datetime.now() # current date and time
date = now.strftime("%m%b%Y")
time = now.strftime("%H:%M:%S")

#  if (False):
if (True):
    with open(seed_name+"_occ.dat",'w') as f:
        f.write(" written on {} at {} \n".format(date,time))
        f.write("        {:d}\n".format(num_wann))
        f.write("        {:d}\n".format(nrpts))
        separator='    '
        for start in range(0, nrpts, 15):
            #  count = min(15, nrpts - start)
            count = min(start+15, nrpts)
            line = separator + separator.join(degen[start:count].astype(str))
            f.write(line)
            f.write('\n')
        counter = 0
        for R_latt in R_vectors:
            for j in range(num_wann):
                for i in range(num_wann):
                    f.write("   {: d}   {: d}   {: d}   {: d}   {: d}   {: 9.6f} {:9.6f} \n".format(R_latt[0],R_latt[1],R_latt[2],i+1,j+1,WF_occ[j,i,counter].real,WF_occ[j,i,counter].imag))
            counter += 1


#-------------------------------------------------------------------------------
if prm.verbose == True:
    idx = np.where(np.all(R_vectors == (0, 0, 0), axis=1))[0]
    print(''' ''')
    print('''         WanOCC v1.1.0''')
    print(''' ''')
    print('''#------------------------------''')
    print(''' Orbital index      Occupation''')
    print('''#------------------------------''')
    total_charge = 0
    for i in range(num_wann):
        total_charge += np.real(WF_occ[i,i,idx])[0]
        print("   {:6d}       |   {:8.4f}".format(i+1,np.real(WF_occ[i,i,idx])[0]))
    print("#------------------------------")
    print(" Charge from diagonal: {:6.4f}".format(total_charge))
    print("#------------------------------")



# sum over K point for occ_mat_pre_F to get back to real space.
# this is to show the origina rotation matrix doesn't have the preserved.
# WF_occ = np.zeros(3,dtype=complex)
# for ikpt in range(nkpt):
#     for i in range(3):
#         WF_occ[i] += occ_mat_pre_f[ikpt][i,i]
# print np.real(WF_occ[:]/64.)
