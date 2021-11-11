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
parser.add_argument('-R', action="store", default=[0,0,0], dest="R_latt",
                    help='Which R_latt to use? Default: 0,0,0')
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
    return np.array(selection, dtype=np.int)


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
tree = ET.parse('si.xml')
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
# assuming k points are evenly spaced
R_lim = np.zeros(3,dtype=int)
for dir in range(3):
    smallest_k = np.sort(list(set(np.abs(k_points[:,dir]))))[1]
    R_lim[dir]=np.floor(1/smallest_k)

# sanity check for R latt
try:
    R_latt = prm.R_latt.split()
    R_latt = np.array(R_latt, dtype=int)
except ValueError:
    raise ValueError("**ERROR: R lattice must be integers")
    # print("**ERROR: R lattice must be integers")

if all(R_latt < R_lim):
    pass
else:
    raise ValueError("R lattice no commensurate with K lattice. (R needs to be smaller than {})".format(R_lim))

#-------------------------------------------------------------------------------
if prm.spin == "up":
    occupation_mat_up=np.zeros([len(occupation_mat),len(occupation_mat[0])/2-tot_bnd_exc],dtype=np.float)
    for i in range(len(occupation_mat)):
        i_bnd_tot = 0
        for i_bnd in range(len(occupation_mat[0])/2):
            if i_bnd+1 in bnd_exclude:
                continue
            occupation_mat_up[i,i_bnd_tot]=occupation_mat[i][i_bnd]
            i_bnd_tot += 1
elif prm.spin == "down":
    occupation_mat_up=np.zeros([len(occupation_mat),len(occupation_mat[0])/2-tot_bnd_exc],dtype=np.float)
    for i in range(len(occupation_mat)):
        i_bnd_tot = 0
        for i_bnd in range(len(occupation_mat[0])/2):
            if i_bnd+1 in bnd_exclude:
                continue
            occupation_mat_up[i,i_bnd_tot]=occupation_mat[i][len(occupation_mat[0])/2+i_bnd]
            i_bnd_tot += 1
else:
    occupation_mat_up=np.zeros([len(occupation_mat),len(occupation_mat[0])-tot_bnd_exc],dtype=np.float)
    for i in range(len(occupation_mat)):
        i_bnd_tot = 0
        for i_bnd in range(len(occupation_mat[0])):
            if i_bnd+1 in bnd_exclude:
                continue
            occupation_mat_up[i,i_bnd_tot]=occupation_mat[i][i_bnd]
            i_bnd_tot += 1


if prm.disentangle:
    # get me k point list and u_matrix_opt(disentanglement)
    k_list=[]
    with open(seed_name+"_u_dis.mat") as f:
        header = f.readline()
        nkpt, num_wann, num_bnd = [int(x) for x in f.readline().split()]
        f.readline()

        u_matrix_opt = np.zeros([nkpt, num_wann, num_bnd], dtype=np.complex)

        for ikpt in range(nkpt):
            k_list.append(f.readline().split())
            for iwann in range(num_wann):
                for ibnd in range(num_bnd):
                    u_real, u_img = [np.float(x) for x in f.readline().split()]
                    u_matrix_opt[ikpt, iwann, ibnd] = complex(u_real, u_img)
            f.readline()


# get me k point list and u_matrix(MLWF)
k_list=[]
with open(seed_name+"_u.mat") as f:
    header = f.readline()
    nkpt, num_wann, num_bnd = [int(x) for x in f.readline().split()]
    f.readline()

    u_matrix = np.zeros([nkpt, num_wann, num_bnd], dtype=np.complex)

    for ikpt in range(nkpt):
        k_list.append(f.readline().split())
        for iwann in range(num_wann):
            for ibnd in range(num_bnd):
                u_real, u_img = [np.float(x) for x in f.readline().split()]

                u_matrix[ikpt, iwann, ibnd] = complex(u_real, u_img)
        f.readline()


# get me occ_mat under wannier basis
#-----------------------------------
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

WF_occ = np.zeros(num_wann,dtype=np.complex)
for ikpt in range(nkpt):
    for i in range(num_wann):
        WF_occ[i] += occ_mat[ikpt][i,i]


#-------------------------------------------------------------------------------
print('''WanOCC v1.0.0
#------------------------------''')
print("R lattice: {}".format(R_latt))
print('''#------------------------------
Orbital index       Occupation
#------------------------------''')



for i in range(len(WF_occ)):
    print("{:6d}          |   {:8.4f}".format(i+1,np.real(WF_occ[i]/nkpt)))

print("#------------------------------")



# sum over K point for occ_mat_pre_F to get back to real space.
# this is to show the origina rotation matrix doesn't have the preserved.
# WF_occ = np.zeros(3,dtype=np.complex)
# for ikpt in range(nkpt):
#     for i in range(3):
#         WF_occ[i] += occ_mat_pre_f[ikpt][i,i]
# print np.real(WF_occ[:]/64.)
