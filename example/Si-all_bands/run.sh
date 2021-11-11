pw.x < Si.scf > Si.scf.out
pw.x < Si.nscf > Si.nscf.out
wannier90.x -pp Si
pw2wannier90.x < Si.pw2wan > Si.pw2wan.out
wannier90.x Si
wanocc.py -seed_name Si -spin unpolarized -dis -R "0 0 0"
