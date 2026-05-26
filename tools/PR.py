import numpy as np
import argparse

"""
This script reads the Wannier90 seedname_r.dat file, and the occ.dat file,
and calculates the band energy using the formula:
$E = \\sum_R 1/(d_R)^2 Tr[P(R) * r(R)^\\dagger]$
where r(R) is the position matrix in real space, P(R) is the occupation matrix
in real space, and d_R is the degeneracy of the R point.
"""

# Command line praser
#----------------------------
parser = argparse.ArgumentParser(description='A script to calculate Wannier Occupation.')
parser.add_argument('-seed_name', action="store", default=str('wannier90'), dest="seed_name",
                    help='seed_name')

prm = parser.parse_args()

def read_wannier90_r(filename):

    with open(filename, 'r') as f:
        lines = f.readlines()

    num_wann = int(lines[1])
    nrpts = int(lines[2])

    ndegen_lines = 0
    degeneracies = [1] * nrpts

    start = 3 + ndegen_lines

    # Read hopping data
    hoppings_dict = {}

    line_index = start

    for irpt in range(nrpts):

        # Each R has num_wann^2 matrix elements
        for _ in range(num_wann * num_wann):

            data = lines[line_index].split()

            Rx, Ry, Rz = map(int, data[:3])

            m = int(data[3]) - 1
            n = int(data[4]) - 1

            X_re = float(data[5])
            X_im = float(data[6])
            Y_re = float(data[7])
            Y_im = float(data[8])
            Z_re = float(data[9])
            Z_im = float(data[10])

            R = (Rx, Ry, Rz)

            if R not in hoppings_dict:
                hoppings_dict[R] = np.zeros(
                    (num_wann, num_wann, 3),
                    dtype=complex
                )

            hoppings_dict[R][m, n, 0] = X_re + 1j * X_im
            hoppings_dict[R][m, n, 1] = Y_re + 1j * Y_im
            hoppings_dict[R][m, n, 2] = Z_re + 1j * Z_im

            line_index += 1

    R_vectors = list(hoppings_dict.keys())
    with open("R_vectors.txt", "w") as f:
        for R in R_vectors:
            f.write(f"{R[0]} {R[1]} {R[2]} \n")
    hoppings = [hoppings_dict[R] for R in R_vectors]

    return num_wann, nrpts, degeneracies, R_vectors, hoppings

def read_wannier90_hr(filename):

    with open(filename, 'r') as f:
        lines = f.readlines()

    num_wann = int(lines[1])
    nrpts = int(lines[2])

    # Read degeneracies
    ndegen_lines = (nrpts + 14) // 15

    degeneracies = []
    for i in range(ndegen_lines):
        degeneracies.extend(
            [int(x) for x in lines[3 + i].split()]
        )

    start = 3 + ndegen_lines

    # Read hopping data
    hoppings_dict = {}

    line_index = start

    for irpt in range(nrpts):

        # Each R has num_wann^2 matrix elements
        for _ in range(num_wann * num_wann):

            data = lines[line_index].split()

            Rx, Ry, Rz = map(int, data[:3])

            m = int(data[3]) - 1
            n = int(data[4]) - 1

            re = float(data[5])
            im = float(data[6])

            R = (Rx, Ry, Rz)

            if R not in hoppings_dict:
                hoppings_dict[R] = np.zeros(
                    (num_wann, num_wann),
                    dtype=complex
                )

            hoppings_dict[R][m, n] = re + 1j * im

            line_index += 1

    R_vectors = list(hoppings_dict.keys())
    with open("R_vectors.txt", "w") as f:
        for R in R_vectors:
            f.write(f"{R[0]} {R[1]} {R[2]} \n")
    hoppings = [hoppings_dict[R] for R in R_vectors]

    return num_wann, nrpts, degeneracies, R_vectors, hoppings


if __name__ == "__main__":

    r_file = prm.seed_name + "_r.dat"
    (
        num_wann,
        nrpts,
        degeneracies,
        R_vectors,
        hoppings
    ) = read_wannier90_r(r_file)

    occ_file = prm.seed_name + "_occ.dat"
    (
        num_wann,
        nrpts,
        degeneracies,
        R_vectors,
        occupations
    ) = read_wannier90_hr(occ_file)

    # calculate \sum_R 1/d^2_R Tr[H(R) * P(R)]
    for dir in range(3):
        result = 0.0
        for R_idx in range(nrpts):
            result += (1.0 / (degeneracies[R_idx] ** 2)) * np.trace(
                    np.dot(np.matrix(occupations[R_idx]),np.matrix(hoppings[R_idx][:,:,dir]).H)
            )

        print(f"dipole moment dir {dir}: {np.real(result):+0.5f}")
