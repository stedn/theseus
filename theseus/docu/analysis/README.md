# activnet analysis scripts

These files are used to generate analysis figures

In these files we will assume the following convention for the parameters in allp


allp(:,1)  = zeta - medium viscosity (nN ds /um^2)
allp(:,2)  = L    - filament length (um)
allp(:,3)  = mu_c - filament compressional modulus (nN)
allp(:,4)  = kap  - filament bending modulus (ignored)
allp(:,5)  = lc   - average distance between cross-links (um)
allp(:,6)  = xi   - inter-cross-link slip coeffficient (nN ds /um)
allp(:,7)  = ups  - active cross-link force (nN)
allp(:,8)  = phi  - active cross-link decorating fraction
allp(:,9)  = psi  - active cross-link spatial distribution (see simulation docs)
allp(:,10) = r    - recycling rate (1/ds)
allp(:,11) = sig  - applied 1d stress (nN/um)
allp(:,12) = Dx   - domain size in x dimension
allp(:,13) = Dy   - domain size in y dimension
allp(:,14) = Df   - location of applied force or edge of domain or both
allp(:,15) = Dw   - width of region for forces and constraints (see simulation docs)
allp(:,16) = nonl - multiplication factor between mu_c and mu_e (if positive ignore)

where nN means nanoNewtons, um means micrometers, ds means deciseconds (10^-1 s)


In addition
all times are recorded in ds,
all forces are recorded in nN
all stresses in nN/um
all distances in um
all velocities in um/ds
