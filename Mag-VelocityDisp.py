# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 10:04:20 2019

@author: Angelo Hollett angelokh26@gmail.com

A progrom to analyze the magnitude-velocity relationship.

"""

from astropy.io import fits
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy import optimize
from scipy.optimize import curve_fit
import numpy as np

#-----------------Open the fits file--------------------#
fits = fits.open('subset_fits.fits')
data = Table.read(fits[1], format='fits')
#-------------------------------------------------------#

#data.info()

#----------------Read in data to arrays-----------------#
velDisp = data['velDisp']
raw_velDispErr = data['velDispErr']
fits_mag = data['ABSMAG']
#-------------------------------------------------------#

#-------Select only R magnitudes from the matrix--------#
raw_Rmag = fits_mag[:,2].data
#-------------------------------------------------------#

cor_velDisp = []
Rmag = []
velDispErr = []

for i,j,k in zip(velDisp, raw_Rmag, raw_velDispErr):
    if all ( [(i > 0), (i < 840)] ):
        cor_velDisp.append(i)
        Rmag.append(j)
        velDispErr.append(k)


print (len(cor_velDisp))
print (len(Rmag))
print (len(velDispErr))

log_velDisp = np.log10(cor_velDisp)
        
#------------Fit the data with a linear model-----------#
def model(x,m,b):
    return m*x+b

# LINEAR fit
init_guess = [-3,-16]
fit = curve_fit(model, log_velDisp, Rmag, sigma=velDispErr, p0=init_guess, absolute_sigma=True)

# unpack the results LINEAR
ans,cov = fit
fit_m,fit_b = ans
fit_sm,fit_sb = np.sqrt(np.diag(np.absolute(cov)))

# print the LINEAR fit results:
print("The fited values for the linear fit model are:")
print("Linear slope m: %.2f +/- %.2f"%(fit_m,fit_sm))
print("Y-axis intercept / vertical translation b: %.2f +/- %.2f"%(fit_b,fit_sb))

#---------------------------------plot the figure etc...---------------------------------------#
plt.figure(figsize=(10,9))
ax1 = plt.subplot()

#ax1.plot(log_velDisp, Rmag, 'o', color='darkcyan', markersize=2)
ax1.errorbar(log_velDisp, Rmag, velDispErr,fmt='k.', label="data")
ax1.set_ylim((-30,0))
ax1.set_title('Determination of Corellation Between Magnitudes and Velocity Dispersion')
ax1.set_ylabel('Absolute Magnitude in R band ($M_{r}$)', fontsize=12)
ax1.set_xlabel('Velocity Dispersion (log)', fontsize=12)
#ax1.set_xscale('log')
ax1.xaxis.set_minor_formatter(mticker.ScalarFormatter())

# # initial guess testing
# curve = (-1)*(log_velDisp) - 16
# ax1.plot(log_velDisp, curve, '-r', label='Model')

# plot the LINEAR data and fit results
print("covariance:")
print(cov)

curve = (fit_m)*(log_velDisp) + fit_b
ax1.plot(log_velDisp, curve, '-r', label='Model')

#ax1.plot(t,model(t,fit_m,fit_b), label="model")
ax1.legend()

# compute LINEAR chi-square
# chisq = sum((Rmag - model(velDisp,fit_m,fit_b))**2/velDispErr**2)
# ax1.text(0.5,0.2,"chi-square: %.2f"%chisq,fontweight="bold",
#      horizontalalignment='center',
#      verticalalignment='center',
#      transform = ax1.transAxes)


#plt.savefig("Rmag_velDisp.png", bbox_inches="tight")