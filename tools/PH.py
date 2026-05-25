import numpy as np
import argparse

"""
This script reads the Wannier90 hr.dat file, and the corresponding occ.dat file,
and calculates the band energy using the formula:
$E = \\sum_R 1/(d_R)^2 Tr[P(R) * H(R)^\\dagger]$
where H(R) is the Hamiltonian in real space, P(R) is the occupation matrix in
real space, and d_R is the degeneracy of the R point. 
"""

# Command line praser
#----------------------------
parser = argparse.ArgumentParser(description='A script to calculate Wannier Occupation.')
parser.add_argument('-seed_name', action="store", default=str('wannier90'), dest="seed_name",
                    help='seed_name')

prm = parser.parse_args()

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

    hr_file = prm.seed_name + "_hr.dat"
    (
        num_wann,
        nrpts,
        degeneracies,
        R_vectors,
        hoppings
    ) = read_wannier90_hr(hr_file)

    occ_file = prm.seed_name + "_occ.dat"
    (
        num_wann,
        nrpts,
        degeneracies,
        R_vectors,
        occupations
    ) = read_wannier90_hr(occ_file)

    # calculate \sum_R 1/d^2_R Tr[H(R) * P(R)]
    result = 0.0
    for R_idx in range(nrpts):
        result += (1.0 / (degeneracies[R_idx] ** 2)) * np.trace(
            np.dot(np.matrix(occupations[R_idx]),np.matrix(hoppings[R_idx]).H)
        )

    print("band energy:", np.real(result))
