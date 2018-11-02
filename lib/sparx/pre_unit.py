##################
# unit converter #
##################


### MASS ###
# Mass of sun to kg
Msun2kg = 0.19889E+31 # kg
kg2Msun = 1./ Msun2kg # Msun

### LENGTH ###
# AU to m
AU2MKS = 0.14960E+12 / (100*365.25*24*60*60)
AU2m = 149597870700.0
m2cm = 100.0
# AU to parsec
AU2pc = 1. / 206260.
# AU to cm
AU2cm = AU2m *m2cm
# km to m
km2m = 1e3
# parsec to m
pc2m = 30.857e15
volume_pc2m =  pc2m**3
# m to parsec
m2pc = 1. / pc2m

# parsec to km
pc2km = 30.857e12

### TIME ###
yr2sec = 31536000.



######################
# physical constants #
######################
# mass of molecular hydrogen
mH2 = 2.*1.660538921e-27 # kg
# meam molecular mass, mass / mass_H2 ~ 1.67
MeanMolecularMass = 2. * 1.67 * 1.6726219e-27 # kg
