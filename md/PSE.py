import numpy as np

PSE = {}
# PSE['Element'] = [atomic number, mass, LJ-sigma, LJ-Eps]
PSE['H'] = np.array([1, 1.00794, 'placeholder_sigma', 'placeholer_eps'])
PSE['He']= np.array([2, 4.002602, 'placeholder_sigma', 'placeholer_eps'])
PSE['Li']= np.array([3, 6.941, 'placeholder_sigma', 'placeholer_eps'])
PSE['Na']= np.array([11, 22.98977, 0.385, 0.0022 ])
PSE['K'] = np.array([19, 39.0983, 'placeholder_sigma', 'placeholer_eps'])
PSE['Rb']= np.array([37, 85.4678, 'placeholder_sigma', 'placeholer_eps'])
PSE['Cs']= np.array([55, 132.90545, 'placeholder_sigma', 'placeholer_eps'])
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
PSE['F'] = np.array([9, 18.984032, 'placeholder_sigma', 'placeholer_eps'])
PSE['Cl']= np.array([17, 35.45, 0.385, 2.972])
PSE['Br']= np.array([35, 79.904, 'placeholder_sigma', 'placeholer_eps'])
PSE['I'] = np.array([53, 126.90447, 'placeholder_sigma', 'placeholer_eps'])
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

# Source Masses: 
#https://upload.wikimedia.org/wikipedia/commons/thumb/e/e6/Periodensystem_der_Elemente.svg/2000px-Periodensystem_der_Elemente.svg.png

### Source for Values of sigma and epsilon (Na and Cl)
#I. Gladich et al. / Chemical Physics Letters 489 (2010) 113–117

