import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class emlinefit(object):
    """
    This package will fit a gaussian to an emission line.
    """
    def __init__(self, wavelength, flux, line_l, line_u, fit_type='gaussian'):
        """
        emlinefit requires a wavelength array ('wavelength') and flux array ('flux')
        containing the emission line to be fit
        The line_l and line_u parameters define the lower and upper bounds of the fitting range.
        Fit_type: 'gaussian', 'asymmetric gaussian'
        """
        self.wavelength=wavelength
        self.flux=flux
        self.line_l=line_l
        self.line_u=line_u
        self.fit_type=fit_type
        
    def gaussian(self, wavelength, amp, mu, sigma):
        exponential = -1 * (wavelength - mu)**2 / (2 * sigma**2)
        return amp * np.exp(exponential)

    def asym_gaussian(self, wavelength, A, a_asym, d):
        ind = (wavelength>self.line_l) & (wavelength<self.line_u)
        ind_full = (self.wavelength>self.line_l) & (self.wavelength<self.line_u)
        delta_vel=wavelength[ind]-self.wavelength[ind_full][np.argmax(self.flux[ind_full])]
        return A * np.exp((-delta_vel**2) / (2 * (a_asym*delta_vel + d)**2))

    def cal_fwhm(self):
        return 2*np.sqrt(2*np.log(2))*self.d/(1-2*np.log(2)*self.asym**2)

    def gaussfitting(self):
        ind = (self.wavelength>self.line_l) & (self.wavelength<self.line_u)
        popt,pcov=curve_fit(self.gaussian, self.wavelength[ind], self.flux[ind],
                    bounds=([0,self.line_l,0],[np.max(self.flux[ind]),self.line_u,self.line_u-self.line_l]))
        return popt,pcov
    
    def asymfitting(self):
        ind = (self.wavelength>self.line_l) & (self.wavelength<self.line_u)
        popt,pcov=curve_fit(self.asym_gaussian, self.wavelength[ind], self.flux[ind],
                    p0=[np.max(self.flux[ind]),1,0.4])
        return popt,pcov
    
    def asym_width(self):
        popt, _ = self.asymfitting()
        self.A, self.asym, self.d = popt
        width = self.cal_fwhm()
        return width
    
    def return_result(self):
        if self.fit_type == 'gaussian':
            popt, pcov = self.gaussfitting()
            return {'amplitude': popt[0], 'mean': popt[1], 'stddev': popt[2], 'covariance': pcov}
        elif self.fit_type == 'asymmetric':
            popt, pcov = self.asymfitting()
            self.A, self.asym, self.d = popt
            width = self.cal_fwhm()
            return {'amplitude': self.A, 'asymmetry': self.asym, 'd': self.d, 'width': width, 'covariance': pcov}
        else:
            raise ValueError("fit_type must be either 'gaussian' or 'asymmetric'")
        
    def plot_fit(self):
        ind = (self.wavelength>self.line_l) & (self.wavelength<self.line_u)
        plt.figure(figsize=(8,5))
        plt.plot(self.wavelength[ind], self.flux[ind], label='Data', color='black')
        if self.fit_type == 'gaussian':
            popt, _ = self.gaussfitting()
            plt.plot(self.wavelength[ind], self.gaussian(self.wavelength[ind], *popt), label='Gaussian Fit', color='red')
        elif self.fit_type == 'asymmetric':
            popt, _ = self.asymfitting()
            plt.plot(self.wavelength[ind], self.asym_gaussian(self.wavelength[ind], *popt), label='Asymmetric Gaussian Fit', color='blue')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.title('Emission Line Fit')
        plt.legend()
        plt.show()
        
wavelength=np.genfromtxt('data/wavelength_sample.txt')
flux=np.genfromtxt('data/flux_sample.txt')
emline=emlinefit(wavelength, flux, 1620, 1640, fit_type='asymmetric')
result=emline.return_result()
print(result)
emline.plot_fit()
