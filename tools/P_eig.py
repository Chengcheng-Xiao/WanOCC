import numpy as np
import argparse

'''
This script reads the occupation file seedname_occ.mat, constructs the density
matrix at a secified k-point, and computes its eigenvalues.
'''

# Command line praser
#----------------------------
parser = argparse.ArgumentParser(description='A script to calculate Wannier Occupation.')
parser.add_argument('-seed_name', action="store", default=str('wannier90'), dest="seed_name",
                    help='seed_name')
parser.add_argument('-k', action="store", default=[0,0,0], dest="kpt",
                   help='Which kpt to use? Default: 0,0,0')

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
    hoppings = [hoppings_dict[R] for R in R_vectors]

    return num_wann, nrpts, degeneracies, R_vectors, hoppings


def construct_hk(kvec, R_vectors, hoppings, degeneracies):

    nwann = hoppings[0].shape[0]

    Hk = np.zeros((nwann, nwann), dtype=complex)

    for R, HR, deg in zip(R_vectors, hoppings, degeneracies):

        phase = np.exp(2j * np.pi * np.dot(kvec, R))

        Hk += (HR / deg) * phase

    return Hk


def eigenvalues_at_k(hr_file, kvec):

    (
        num_wann,
        nrpts,
        degeneracies,
        R_vectors,
        hoppings
    ) = read_wannier90_hr(hr_file)

    Hk = construct_hk(
        kvec,
        R_vectors,
        hoppings,
        degeneracies
    )

    eigvals = np.linalg.eigvalsh(Hk)

    return np.real(eigvals)


if __name__ == "__main__":

    hr_file = prm.seed_name + "_occ.dat"

    # Example k-point (fractional coordinates)
    #  kvec = [0.0, 0.0, 0.0]
    kvec = [float(x) for x in prm.kpt.split(",")]

    eigvals = eigenvalues_at_k(hr_file, kvec)

    print("Eigenvalues at k =", kvec)
    print(np.sort(eigvals)[::-1])
