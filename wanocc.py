#!/usr/bin/env python
# coding: utf-8

import argparse
import xml.etree.ElementTree as ET
import numpy as np

# Command line praser
#----------------------------
parser = argparse.ArgumentParser(description='A script to calculate Wannier Occupation.')
parser.add_argument('-seed_name', action="store", default=str('wannier90'), dest="seed_name",
                    help='seed_name')
parser.add_argument('-spin', action="store_true", default=False, dest="spin",
                    help='spin polarized? (Default=False)')
parser.add_argument('-bnd_exc', action="store", default=0, dest="bnd_exclude",
                    help='Get distance? (Default=False)')
parser.add_argument('-dis', action="store_true", default=False, dest="disentangle",
                    help='Do we want to use disentanglement procedure?')

prm = parser.parse_args()

#%% input
bnd_exclude = int(prm.bnd_exclude)
seed_name = prm.seed_name


if prm.spin:
    # get me starting index of ks_energies in root
    tree = ET.parse(seed_name+'.xml')
    root = tree.getroot()

    starting_index = [i.tag for i in root[3][-1]].index('ks_energies')

    # get me occupation matrix
    occupation_mat = []
    for i in root[3][-1][starting_index:]:
    #     print "kpoints: ", i[0].text
        occupation_mat.append(i[3].text.split())

    occupation_mat_up=np.zeros([len(occupation_mat),len(occupation_mat[0])/2-bnd_exclude],dtype=np.float)
    occupation_mat_dn=np.zeros([len(occupation_mat),len(occupation_mat[0])/2-bnd_exclude],dtype=np.float)
    for i in range(len(occupation_mat)):
        occupation_mat_up[i]=occupation_mat[i][0+bnd_exclude:len(occupation_mat[0])/2]
        occupation_mat_dn[i]=occupation_mat[i][len(occupation_mat[0])/2+bnd_exclude:]
else:
    # get me starting index of ks_energies in root
    tree = ET.parse(seed_name+'.xml')
    root = tree.getroot()

    starting_index = [i.tag for i in root[3][-1]].index('ks_energies')

    # get me occupation matrix
    occupation_mat = []
    for i in root[3][-1][starting_index:]:
    #     print "kpoints: ", i[0].text
        occupation_mat.append(i[3].text.split())

    occupation_mat_up=np.zeros([len(occupation_mat),len(occupation_mat[0])-bnd_exclude],dtype=np.float)
    for i in range(len(occupation_mat)):
        occupation_mat_up[i]=occupation_mat[i][0+bnd_exclude:len(occupation_mat[0])]

if prm.disentangle:
    # get me k point list and u_matrix_opt(disentanglement)
    k_list=[]
    with open(seed_name+"_up_u_dis.mat") as f:
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
with open(seed_name+"_up_u.mat") as f:
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
    # pre-rotate occ_matrix with u_matrix_opt
    occ_mat_pre = np.dot(np.dot(u_matrix_opt[ikpt],a),np.matrix(u_matrix_opt[ikpt]).H)
    occ_mat_pre_f.append(occ_mat_pre)
    # rotate occ_matrix with u_matrix
    occ_mat.append(np.dot(np.dot(u_matrix[ikpt],occ_mat_pre),np.matrix(u_matrix[ikpt]).H))

WF_occ = np.zeros(3,dtype=np.complex)
for ikpt in range(nkpt):
    for i in range(3):
        WF_occ[i] += occ_mat[ikpt][i,i]
if prm.spin:
    print "spin[up]:", np.real(WF_occ[:]/64.)
else:
    print "WanOCC: ", np.real(WF_occ[:]/64.)


# sum over K point for occ_mat_pre_F to get back to real space.
# this is to show the origina rotation matrix doesn't have the preserved.
# WF_occ = np.zeros(3,dtype=np.complex)
# for ikpt in range(nkpt):
#     for i in range(3):
#         WF_occ[i] += occ_mat_pre_f[ikpt][i,i]
# print np.real(WF_occ[:]/64.)

if prm.spin:
    if prm.disentangle:
        # get me k point list and u_matrix_opt(disentanglement)
        k_list=[]
        with open(seed_name+"_dn_u_dis.mat") as f:
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
    with open(seed_name+"_dn_u.mat") as f:
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
        a = np.zeros([occupation_mat_dn[0].size, occupation_mat_dn[0].size])
        np.fill_diagonal(a,occupation_mat_dn[ikpt])
        # pre-rotate occ_matrix with u_matrix_opt
        occ_mat_pre = np.dot(np.dot(u_matrix_opt[ikpt],a),np.matrix(u_matrix_opt[ikpt]).H)
        occ_mat_pre_f.append(occ_mat_pre)
        # rotate occ_matrix with u_matrix
        occ_mat.append(np.dot(np.dot(u_matrix[ikpt],occ_mat_pre),np.matrix(u_matrix[ikpt]).H))

    WF_occ = np.zeros(3,dtype=np.complex)
    for ikpt in range(nkpt):
        for i in range(3):
            WF_occ[i] += occ_mat[ikpt][i,i]
    if prm.spin:
        print "spin[dn]:", np.real(WF_occ[:]/64.)

    # sum over K point for occ_mat_pre_F to get back to real space.
    # this is to show the origina rotation matrix doesn't have the preserved.
    # WF_occ = np.zeros(3,dtype=np.complex)
    # for ikpt in range(nkpt):
    #     for i in range(3):
    #         WF_occ[i] += occ_mat_pre_f[ikpt][i,i]
    # print np.real(WF_occ[:]/64.)
