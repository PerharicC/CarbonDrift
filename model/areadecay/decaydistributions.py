import numpy as np
import numpy.random as rnd
import warnings
import matplotlib.pyplot as plt

#TODO Other functions

class Decay_Functions:

    """A class of methods corresponding to different decay distributions.

    The distributions can be mass and/or z dependent."""

    def __init__(self, m0=1, H=1000, mass_threshold = 1, z_threshold = 1, threshold = 1):

        """"""
        self.m0 = m0
        self.H = H
        self.mt = mass_threshold
        self.zt = z_threshold
        self.t = threshold

    def mass_sqrt(self):
        return lambda m, z : self.mt * np.sqrt(m / self.m0)
    
    def z_linear(self):
        return lambda m, z : self.zt * (z / self.H + 1)
    
    def mass_sqrt_z_linear(self):

        if self.t == 1 and self.mt == 1 and self.zt == 1:
            warnings.warn("If all thresholds are set to 1, the inital particle will always decay.")
        
        return lambda m, z : self.t * (self.z_linear()(z) + self.mass_sqrt()(m)) / (self.mt + self.zt)

    def choose_distribution(self, method_name, plot = False):
        if method_name == None:
            raise NameError("No decay distribution was given.")
        if plot:
            self.plot(method_name)
        if hasattr(Decay_Functions, method_name):
            return getattr(Decay_Functions, method_name)
        else:
            errormsg = " does not exist. Try:\n"
            method_names = [method_name for method_name in dir(Decay_Functions) if callable(getattr(Decay_Functions, method_name))]
            method_names.remove("choose_distribution")
            method_names.remove("plot")
            for name in method_names:
                if name[0] == "_":
                    continue
                errormsg += name + "\n"
            raise NameError("Decay distribution:" + method_name + errormsg)

    def plot(self, method_name):

        """Plot the selected distribution."""

        mass = np.linspace(0, self.m0, 100)
        Z = np.linspace(-self.H, 0, 100)
        if method_name == "mass_sqrt_z_linear":

            for z in np.linspace(-self.H, 0, 10):
                plt.plot(mass, self.mass_sqrt_z_linear()(mass, z), label = "z = {}".format(z))
                plt.xlabel("mass")

            plt.legend()
            plt.show()
            plt.close()
            
            for m in np.linspace(0, self.m0, 10):
                plt.plot(Z, self.mass_sqrt_z_linear()(m, Z), label = "mass = {}".format(m))
                plt.xlabel("z")
            
            plt.legend()
            plt.show()
            plt.close()
        else:
            if "mass" in method_name:
                plt.plot(mass, getattr(Decay_Functions, method_name)(self)(mass, 1))
                plt.xlabel("mass")
            else:
                plt.plot(Z, getattr(Decay_Functions, method_name)(self)(1, Z))
                plt.xlabel("z")
            plt.show()

