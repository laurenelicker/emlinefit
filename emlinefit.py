class emlinefit(object):
    """
    This package will fit a gaussian to an emission line
    """
    def __init__(self, wavelength, flux):
        """
        emlinefit requires a wavelength array ('wavelength') and flux array ('flux')
        containing the emission line to be fit
        """
        
        