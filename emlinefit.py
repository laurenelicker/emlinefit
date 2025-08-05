import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class emlinefit(object):
    """
    This package will fit a gaussian to an emission line.
    """
    def __init__(self, wavelength, flux, line_l, line_u):
        """
        emlinefit requires a wavelength array ('wavelength') and flux array ('flux')
        containing the emission line to be fit
        """
        self.wavelength=wavelength
        self.flux=flux
        self.line_l=line_l
        self.line_u=line_u
        
    def gaussian(self, amp, mu, sigma):
        exponential = -1 * (self.wavelength - mu)**2 / (2 * sigma**2)
        return amp * np.exp(exponential)

    def asym_gaussian(l, A, a_asym, d):
        delta_vel=l-lpeak
        return A * np.exp((-delta_vel**2 / (2 * (a_asym*delta_vel + d)**2))

    def cal_fwhm(asym,d):
        return 2*np.sqrt(2*np.log(2))*d/(1-2*np.log(2)*asym**2)

    def fitting(self):
        ind = (self.wavelength>self.line_l) & (self.wavelength<self.line_u)
        popt,pcov=curve_fit(gaussian, self.wavelength[ind], self.flux[ind],
                    bounds=([0,self.line_l,0],[np.max(self.flux[ind]),self.line_u,self.line_u-self.line_l]))
        return popt,pcov
