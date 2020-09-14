# WANOCC
Python script to calculate Wannier occupation matrix from Maximally localized Wannier functions.
<!-- [`Wannier90`](http://www.wannier.org/). -->

## Input files
1. u matrix:`seed_name_u.mat`

   [optional] `seed_name_u_dis.mat`

   [optional] `seed_name_up[dn]_u_dis.mat`

2. DFT occupation matrix is read from `seed_name.xml`

## Command line option
See all option by `wanocc.py -h`

- seed_name, seed name for all files. Default = "wannier90".
- dis, disentanglement used? Default = False. Choose from: True/False
- spin, which spin channel to calculate? Default = unpolarized. Choose from: up/down/unpolarized
- bnd_exc, which bands to exclude. Default = 'empty'


## Usage
```
occtest.py -seed_name test -spin -bnd_exc 26 -dis
```
