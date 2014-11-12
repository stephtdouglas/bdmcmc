import datetime
import logging

## Third-party
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np

sigma_0 = 1e-16
lns = -35.8
s = np.exp(lns)

def gaussian(sigma):
    x = np.arange(-5e-15,5e-15,1e-19)
    y = np.exp(-0.5*x**2/sigma**2)/(sigma*np.sqrt(2*np.pi))
    return x,y

a,b = gaussian(sigma_0)
d,f = gaussian(np.sqrt(s**2 + sigma_0**2))

plt.plot(a,b/sum(b)/1e-19,label='orig')
plt.plot(d,f/sum(f)/1e-19,label='+tolerance')
plt.legend()
