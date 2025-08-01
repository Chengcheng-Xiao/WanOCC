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


#%% input
if prm.bnd_exclude != None:
    bnd_exclude = parseIntSet(prm.bnd_exclude)
    tot_bnd_exc = int(bnd_exclude.shape[0])
else:
    bnd_exclude = []
    tot_bnd_exc = 0
seed_name = prm.seed_name

#-------------------------------------------------------------------------------
# get me starting index of ks_energies in root
tree = ET.parse(seed_name+'.xml')
root = tree.getroot()

starting_index = [i.tag for i in root[3][-1]].index('ks_energies')

# get me occupation matrix
occupation_mat = []
for i in root[3][-1][starting_index:]:
    occupation_mat.append(i[3].text.split())

#-------------------------------------------------------------------------------
# get me starting index of kpoints in root
tree = ET.parse(seed_name+'.xml')
root = tree.getroot()

# not sure if this index can change in QE but here it is
starting_index = [i.tag for i in root[2][8]].index('k_point')

# get me kpoints
k_points = []
for i in root[2][8][starting_index:]:
    k_points.append(np.array(i.text.split(),dtype=float))

# convert data format
k_points = np.array(k_points)

#%% Calculate R limit with k points
# NOTE: assuming k points are evenly spaced
R_lim = np.zeros(3,dtype=int)
for dir in range(3):
    if len(np.sort(list(set(np.abs(k_points[:,dir])))))>1: # multi-kpt
       smallest_k = np.sort(list(set(np.abs(k_points[:,dir]))))[1]
       R_lim[dir]=np.round(1/smallest_k)-1
    else: # gamma only
       R_lim[dir]=1

# sanity check for R latt
# try:
#     R_latt = prm.R_latt.split()
#     R_latt = np.array(R_latt, dtype=int)
# except ValueError:
#     raise ValueError("**ERROR: R lattice must be integers")
#     # print("**ERROR: R lattice must be integers")
#
# if all(R_latt < R_lim):
#     pass
# else:
#     raise ValueError("R lattice no commensurate with K lattice. (R needs to be smaller than {})".format(R_lim))

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

        u_matrix_opt = np.zeros([nkpt, num_wann, num_bnd], dtype=complex)

        for ikpt in range(nkpt):
            k_list.append(f.readline().split())
            for iwann in range(num_wann):
                for ibnd in range(num_bnd):
                    u_real, u_img = [float(x) for x in f.readline().split()]
                    u_matrix_opt[ikpt, iwann, ibnd] = complex(u_real, u_img)
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

                u_matrix[ikpt, iwann, ibnd] = complex(u_real, u_img)
        f.readline()


# get me occ_mat under wannier basis
#-----------------------------------
# R_latt = [0,0,0]
WF_occ = np.zeros([num_wann,num_wann,2*R_lim[0]+1,2*R_lim[1]+1,2*R_lim[2]+1],dtype=complex)
for R1 in range(-R_lim[0],R_lim[0]+1):
    for R2 in range(-R_lim[1],R_lim[1]+1):
        for R3 in range(-R_lim[2],R_lim[2]+1):
            R_latt = [R1,R2,R3]
            occ_mat=[]
            occ_mat_pre_f = []
            for ikpt in range(nkpt):
                a = np.zeros([occupation_mat_up[0].size, occupation_mat_up[0].size])
                np.fill_diagonal(a,occupation_mat_up[ikpt])
                if prm.disentangle:
                    # pre-rotate occ_matrix with u_matrix_opt
                    occ_mat_pre = np.dot(np.dot(u_matrix_opt[ikpt],a),np.matrix(u_matrix_opt[ikpt]).H)
                    occ_mat_pre_f.append(occ_mat_pre)
                    # calculate the k factor
                    phase_factor = np.exp(complex(0,np.dot(k_points[ikpt],R_latt)*2*np.pi))
                    # rotate occ_matrix with u_matrix
                    occ_mat.append(phase_factor*np.dot(np.dot(u_matrix[ikpt],occ_mat_pre),np.matrix(u_matrix[ikpt]).H))
                else:
                    # calculate the k factor
                    phase_factor = np.exp(complex(0,np.dot(k_points[ikpt],R_latt)*2*np.pi))
                    # rotate occ_matrix with u_matrix
                    occ_mat.append(phase_factor*np.dot(np.dot(u_matrix[ikpt],a),np.matrix(u_matrix[ikpt]).H))

            # summing up kpoints
            for ikpt in range(nkpt):
                for i in range(num_wann):
                    for j in range(num_wann):
                        WF_occ[i,j,R1,R2,R3] += occ_mat[ikpt][i,j]/nkpt

#-------------------------------------------------------------------------------
from datetime import datetime
now = datetime.now() # current date and time
date = now.strftime("%m%b%Y")
time = now.strftime("%H:%M:%S")

with open(seed_name+"_occ.mat",'w') as f:
    f.write("written on {} at {} \n".format(date,time))
    f.write("        {:d}\n".format(num_wann))
    f.write("        {:d}\n".format((2*R_lim[0]+1)*(2*R_lim[1]+1)*(2*R_lim[2]+1)))
    for i in range(num_wann):
        for j in range(num_wann):
            for R1 in range(-R_lim[0],R_lim[0]+1):
                for R2 in range(-R_lim[1],R_lim[1]+1):
                    for R3 in range(-R_lim[2],R_lim[2]+1):
                        f.write("   {: d}   {: d}   {: d}   {: d}   {: d}  {: 9.6f}  {: 9.6f} \n".format(i,j,R1,R2,R3, WF_occ[i,j,R1,R2,R3].real,WF_occ[i,j,R1,R2,R3].imag))


#-------------------------------------------------------------------------------
if prm.verbose == True:
    print(''' ''')
    print('''         WanOCC v1.0.0''')
    print(''' ''')
    print('''#------------------------------''')
    print(''' Orbital index      Occupation''')
    print('''#------------------------------''')
    total_charge = 0
    for i in range(num_wann):
        total_charge += np.real(WF_occ[i,i,0,0,0])
        print("   {:6d}       |   {:8.4f}".format(i+1,np.real(WF_occ[i,i,0,0,0])))
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
