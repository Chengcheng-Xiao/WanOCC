# WANOCC
Python script to calculate Wannier occupation matrix from Maximally localized Wannier functions.
<!-- [`Wannier90`](http://www.wannier.org/). -->

## Theory
Rotating the DFT occupation matrix using the unitary matrix calculated by Wannier90 to obtain Wannier occupation matrix:

$$
\mathscr{F}_{n n'}^{\mathbf{R} \mathbf{R'}} = \sum_{m \mathbf{k}} e^{i\mathbf{k}(\mathbf{R}-\mathbf{R'})} U_{m n'}^{(\mathbf{k})} {U_{m n}^{(\mathbf{k})}}^* f_{m \mathbf{k}}
$$

<!-- ![image](WanOcc_eq.svg) -->

## Usage
```
wanocc.py -seed_name test -spin unpolarized -bnd_exc 5-12 -dis
```

## Command line option
See all option by `wanocc.py -h`

- `seed_name`, seed name for all files. Default = "wannier90".
<!-- - `R`, which R_latt to use. Default: 0,0,0 -->
- `dis`, enable subspace selection. Default = False
- `spin`, which spin channel to calculate? Default = unpolarized. Choose from: up/down/unpolarized
- `bnd_exc`, which bands to exclude. Default = 'empty'. Format "10-12 15"
- `v`, verbose mode, print out charge on WFs in home cell. Default = False

## Input files
1. U matrix:`seed_name_u.mat`

   [Optional] subspace selection: `seed_name_u_dis.mat`

   [Optional] spin polarization: `seed_name_up[dn]_u_dis.mat`

2. DFT occupation matrix is read from `seed_name.xml` [from Quantum espresso]

## Output file
The Wannier occupation matrix will be print in `seed_name_occ.mat` with the following format:

- First line: comment
- Second line: number of wannier functions
- Third line: number of R lattices
- Forth line to EOF: 1st WF index | 2nd WF index | R1 | R2 | R3 | Occupaiton (real) | Occupation (imag).

