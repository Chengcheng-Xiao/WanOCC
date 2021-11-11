# Si example
A simple example of calculating Wannier occupation

- Si-all_bands: calculates the occupation number for all 8 sp3 orbitals.
- Si-valence_bands: calculates the occupation number for 4 valence sp3 orbitals.

## Run
### Si-all_bands
1. Run `pw.x < Si.scf > Si.scf.out` to generate charge density.
2. Run `pw.x < Si.nscf > Si.nscf.out` to generate wavefunctions on a uniform k grid.
3. Run `wannier90.x -pp Si` to generate k-points neighbors.
4. Run `pw2wanier90.x` to generate `.amn` and `.mmn` files.
5. Run `wanier90.x Si` to generate `Si_u.mat` and `Si_u_dis.mat` files.
6. Run `wanocc.py -seed_name Si -spin unpolarized -dis -R "0 0 0"` to get Wannier occupations


### Si-valence_bands
1. Run `pw.x < Si.scf > Si.scf.out` to generate charge density.
2. Run `pw.x < Si.nscf > Si.nscf.out` to generate wavefunctions on a uniform k grid.
3. Run `wannier90.x -pp Si` to generate k-points neighbors.
4. Run `pw2wanier90.x` to generate `.amn` and `.mmn` files.
5. Run `wanier90.x Si` to generate `Si_u.mat` and `Si_u_dis.mat` files.
6. Run `./wanocc.py -seed_name Si -spin unpolarized -bnd_exc 5-12 -R "0 0 0"` to get Wannier occupations
