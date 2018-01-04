# Grid Type of Zeus data: SPH, CYL, or REC
GridType =              'SPH'
# data type : BINARY or ASCII
DataType =              'BINARY'
# The path of the directory of Zeus Data
#ZeusDataDir =           '/tiara/ara/data/outflows/zeustw2sparx/ly/n=6_case/'
ZeusDataDir =           '/tiara/ara/data/ithsieh/datab335'
# the time stamp of the files (INTEGER)
#TimeStamp =             10
TimeStamp =             35
# number of merged cell, to reduce resolution of one dimension,
# the value 0 or 1 will do nothing on griding
Nmerged =               8
# use the post-processed density profile (in the file  'o_d__'+time_stamp+'_vThr')
# for the outflow case
PostProcessedDensity =  False
# the constantant temperature in unit of Kelvin, 
# use the file 'o_T__'+'%5d' % TimeStamp if the value is zero or without the attribute
ConstantTemperature =   10.0
# maximum extension in the outflow axis (AU)
Rmax_AU =               12500.0 # AU
# turbulent velocity (m/s)
TurbulentVelocity =     500.0 
# the name of the molecule to be simulated
MolecularSpecie =       'co@xpol'
# the molecular abundance
MolecularAbundance =    1e-3
# The absorption profile of the dust
DustKappa =             'table,jena_thin_e5'
# the ratio between the dust and the gas
DustToGas =             0.01
# the temperature of CMB emission
T_cmb =                 2.73
# Dust polarized efficiency
DustAlpha = 0.15
