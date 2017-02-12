import numpy as np

PSE = {}
# PSE['Element'] = [atomic number, mass, LJ-sigma, LJ-Eps]
# Property --- Unit
# Atmic Number --- N.A.
# Mass --- amu
# LJ-Sigma --- m
# LJ-Eps  --- kcal/mol
PSE['H'] = np.array([1, 1.00794, 'placeholder_sigma', 'placeholer_eps'])
PSE['He']= np.array([2, 4.002602, 'placeholder_sigma', 'placeholer_eps'])
PSE['Li']= np.array([3, 6.941, 1.715, 0.05766])
PSE['Na']= np.array([11, 22.98977, 2.497, 0.078226 ])
PSE['K'] = np.array([19, 39.0983, 3.184, 0.1183])
PSE['Rb']= np.array([37, 85.4678, 3.302, 0.2405])
PSE['Cs']= np.array([55, 132.90545, 3.440, 0.5013])
PSE['Fr']= np.array([87, 2223.0197, 'placeholder_sigma', 'placeholer_eps'])
PSE['Be']= np.array([4, 9.012182, 'placeholder_sigma', 'placeholer_eps'])
PSE['Mg']= np.array([12, 24.3050, 'placeholder_sigma', 'placeholer_eps'])
PSE['Ca']= np.array([20, 40.078, 'placeholder_sigma', 'placeholer_eps'])
PSE['Sr']= np.array([38, 87.62, 'placeholder_sigma', 'placeholer_eps'])
PSE['Ba']= np.array([56, 137.327, 'placeholder_sigma', 'placeholer_eps'])
PSE['Ra']= np.array([88, 226.08, 'placeholder_sigma', 'placeholer_eps'])
PSE['O'] = np.array([8, 15.999, 'placeholder_sigma', 'placeholer_eps'])
PSE['S'] = np.array([16, 32.06, 'placeholder_sigma', 'placeholer_eps'])
PSE['Se']= np.array([34, 78.96, 'placeholder_sigma', 'placeholer_eps'])
PSE['Te']= np.array([52, 127.60, 'placeholder_sigma', 'placeholer_eps'])
PSE['Po']= np.array([84, 209.98, 'placeholder_sigma', 'placeholer_eps'])
PSE['F'] = np.array([9, 18.984032, 3.954, 0.006465])
PSE['Cl']= np.array([17, 35.45, 4.612, 0.02502])
PSE['Br']= np.array([35, 79.904, 4.812, 0.03596])
PSE['I'] = np.array([53, 126.90447, 5.197,0.04220])
PSE['As']= np.array([85, 209.9871, 'placeholder_sigma', 'placeholer_eps'])
PSE['Sc']= np.array([21, 44.95591, 'placeholder_sigma', 'placeholer_eps'])
PSE['Ti']= np.array([22, 47.867, 'placeholder_sigma', 'placeholer_eps'])
PSE['V'] = np.array([23, 50.9415, 'placeholder_sigma', 'placeholer_eps'])
PSE['Cr']= np.array([24, 51.9981, 'placeholder_sigma', 'placeholer_eps'])
PSE['Mn']= np.array([25, 54.938049, 'placeholder_sigma', 'placeholer_eps'])
PSE['Fe']= np.array([26, 55.845, 'placeholder_sigma', 'placeholer_eps'])
PSE['Co']= np.array([27, 58.93320, 'placeholder_sigma', 'placeholer_eps'])
PSE['Ni']= np.array([28, 58.6934, 'placeholder_sigma', 'placeholer_eps'])
PSE['Cu']= np.array([29, 63.546, 'placeholder_sigma', 'placeholer_eps'])
PSE['Zn']= np.array([30, 65.409, 'placeholder_sigma', 'placeholer_eps'])

# Source for Group(I) and Group(XVII):
#A.H. Mao, R.V. Pappu (2012) Crystal Lattice Properties Fully Determine Short-Range Interaction Parameters for Alkali and Halide Ions. Journal of Chemical Physics, 137: 064104(1-9) 
#doi: http://dx.doi.org/10.1063/1.4742068
