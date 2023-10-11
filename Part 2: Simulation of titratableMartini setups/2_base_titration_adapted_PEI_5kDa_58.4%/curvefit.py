 
import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy.optimize import curve_fit
import pandas as pd 

curve=pd.read_csv("protdata.csv")
print(curve)

ph=curve["pH"].to_numpy()
degreeprot=curve["Degree of prot"].to_numpy()

print(ph,degreeprot)
def calc_degree_of_deprot(ph, pka, q):
    return 1/(10**(q*(pka-ph))+1)

popt, _ = curve_fit(calc_degree_of_deprot,ph, degreeprot, p0=(6.0, 1.0),maxfev=1000)

pka,q=popt

plt.scatter(curve["pH"],curve["Degree of prot"])

x_line = np.arange(min(curve["pH"]), max(curve["pH"]), 1)
y_line = calc_degree_of_deprot(x_line, pka, q)
plt.plot(x_line, y_line, '--', color='red')
print('y = %.5f * x + %.5f' % (pka, q))
print( popt)
plt.show()
