import numpy as np
import pandas as pd
import matplotlib.pyplot as mplt

data = pd.read_table("dat.txt",names=["x","p","z"], sep=' ')
print(data)
t = np.arange(0.1,10.1,0.1)
x = np.sin(t)
p = np.cos(t)
mplt.plot(data["x"].values,data["p"].values,'--',label='integration')
mplt.plot(x,p,label ='real')
mplt.legend()
mplt.show()
mplt.plot(t,data["x"].values,'--',label ='integration')
mplt.plot(t,x,label = "real")
mplt.legend()
mplt.show()
