import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy.optimize import curve_fit
import pandas as pd 


def calc_degree_of_deprot(ph, pka, q):
    return 1/(10**(q*(pka-ph))+1)


parser = argparse.ArgumentParser(description='Plot the results')
parser.add_argument('-f', dest='file', type=str, default='results.txt', help='results file (.txt)')
parser.add_argument('-o', dest='out', type=str, default='results.pdf', help='output file (.pdf)')
args = parser.parse_args()

ph, degree_of_deprot, std = np.loadtxt(args.file, unpack=True)
degree_of_prot=1-degree_of_deprot

fit = curve_fit(calc_degree_of_deprot, ph, degree_of_prot, p0=(5, 3.0),maxfev=2000)[0]
ph_range = np.arange(3, 8.1, 0.1)

fig, ax1 = plt.subplots()

ax2 = ax1.twinx()

degree_of_prot=1-degree_of_deprot

ax1.plot(ph_range, calc_degree_of_deprot(ph_range, *fit),marker="o",color="green")
ph_ex=[10.198,9.68,9.244,8.772,8.054,7.094,5.93,4.464,3.052,2.305,2.11]
hcl_ex=[1,2,3,4,5,6,7,8,9,10,11]
ax2.plot(ph_ex,hcl_ex,marker="^")

ax1.set_xlabel("pH")
ax1.set_ylabel('Degree of Protonation')
ax2.set_ylabel('V HCl')
ax1.set_ylim(0,0.55)

fig.tight_layout()
fig.savefig(args.out, bbox_inches='tight')

