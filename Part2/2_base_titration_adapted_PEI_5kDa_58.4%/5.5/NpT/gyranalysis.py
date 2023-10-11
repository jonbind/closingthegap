import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt

gyr=pd.read_csv("gyr.xvg",delim_whitespace=True,skiprows=26)
print(gyr)
gyr.plot(y="@",use_index=True)
mean=gyr["@"].mean()
f = open("../../gyr.txt", "a")
f.write("{},".format(mean))
f.close()
