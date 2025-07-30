import numpy as np
import os
from scipy.interpolate import UnivariateSpline

class MaterialProperties:
    """Material properties and refractive index calculations"""
    
    def __init__(self):
        self.c = 299792458  # Speed of light in m/s
        self.h = 6.626e-34  # Planck constant
        self.eV = 1.602e-19  # eV to Joules
    
    def load_refractive_index_data(self, filename):
        """Load refractive index data from file"""
        try:
            data = np.loadtxt(filename)
            wavelength = data[:, 0] * 1000  # Convert to nm
            n_real = data[:, 1]
            n_imag = data[:, 2]
            return wavelength, n_real, n_imag
        except Exception as e:
            print(f"Error loading {filename}: {e}")
            return None, None, None
    
    def nComplex_interp(self, filename, wavelengths):
        """Python equivalent of nComplexInterp.m"""
        try:
            data = np.loadtxt(filename)
            wl_data = data[:, 0] * 1000  # Convert to nm
            n_real = data[:, 1]
            n_imag = data[:, 2]
            
            # Use spline interpolation (equivalent to MATLAB spline)
            f_real = UnivariateSpline(wl_data, n_real, s=0, ext=3)
            f_imag = UnivariateSpline(wl_data, n_imag, s=0, ext=3)
            
            n_interp = f_real(wavelengths)
            k_interp = f_imag(wavelengths)
            
            return n_interp + 1j * k_interp
        except Exception as e:
            print(f"Error interpolating {filename}: {e}")
            return np.ones_like(wavelengths, dtype=complex)
    
    def silver_refractive_index(self, wavelengths):
        """Calculate silver refractive index using polynomial fit"""
        wl = np.array(wavelengths)
        eps_real = (-14.463366 + 0.13725*wl - 3.91145e-4*wl**2 + 
                   3.31078e-7*wl**3 - 1.13012e-10*wl**4)
        eps_imag = (2.80407 - 0.01862*wl + 4.73233e-5*wl**2 - 
                   4.93745e-8*wl**3 + 1.83724e-11*wl**4)
        
        return np.sqrt(eps_real + 1j*eps_imag)
    
    def gold_refractive_index(self, wavelengths):
        """Gold refractive index from file or fallback"""
        return self.nComplex_interp('nAu.txt', wavelengths)
    
    def lorentz_ws2(self, wavelengths):
        """WS2 refractive index using exact MATLAB parameters"""
        wl = np.array(wavelengths)
        
        # Energy in J
        E = self.h * self.c / (wl * 1e-9)
        # Energy in eV
        E = E / self.eV
        
        # Parameters from MATLAB code
        fj = np.array([1.59, 0.7, 2.95, 2.8, 12])
        Ej = np.array([2.0195, 2.2379, 2.4087, 2.5996, 2.850])
        yj = np.array([0.028, 0.2, 0.15, 0.3, 0.23])
        
        epsilon = np.ones_like(E, dtype=complex)
        
        for i in range(len(fj)):
            epsilon += fj[i] / (Ej[i]**2 - E**2 - 1j*E*yj[i])
        
        ncomplex = np.sqrt(epsilon)
        return ncomplex
    
    def mos2_refractive_index(self, wavelengths):
        """MoS2 refractive index from file"""
        return self.nComplex_interp('nMoS2.txt', wavelengths)
    
    def hbn_refractive_index(self, wavelengths):
        """h-BN refractive index (constant)"""
        return np.full_like(wavelengths, 2.25, dtype=complex)
    
    def pmma_refractive_index(self, wavelengths):
        """PMMA refractive index (constant)"""
        return np.full_like(wavelengths, 1.5, dtype=complex)
    
    def sio2_refractive_index(self, wavelengths):
        """SiO2 refractive index using Sellmeier equation"""
        wl_um = np.array(wavelengths) * 1e-3  # Convert to micrometers
        n_squared = (1 + 0.6961663*wl_um**2/(wl_um**2 - 0.0684043**2) +
                    0.4079426*wl_um**2/(wl_um**2 - 0.1162414**2) +
                    0.8974794*wl_um**2/(wl_um**2 - 9.896161**2))
        return np.sqrt(n_squared)
    
    def silicon_refractive_index(self, wavelengths):
        """Silicon refractive index from file"""
        return self.nComplex_interp('nSi.txt', wavelengths)
    
    def tdbc_refractive_index(self, wavelengths):
        """TDBC refractive index using multiple Lorentz oscillators"""
        E_eV = 1240.0 / np.array(wavelengths)
        
        eps_base = 2.1
        fj9 = 0.008
        fj1 = 0.113
        
        # Multiple Lorentz oscillators as in the MATLAB code
        oscillators = [
            (fj9, 2.09, 0.02), (fj9, 2.10, 0.02), (fj9, 2.105, 0.016),
            (fj9, 2.11, 0.02), (fj9, 2.36, 0.1), (fj9, 2.34, 0.2),
            (fj9, 2.30, 0.1), (fj9, 2.28, 0.09), (fj1, 2.12, 0.05),
            (fj1, 2.13, 0.09), (fj9, 2.14, 0.09), (fj9, 2.15, 0.09),
            (fj9, 2.16, 0.06), (fj9, 2.17, 0.06), (fj9, 2.18, 0.07)
        ]
        
        eps = np.full_like(E_eV, eps_base, dtype=complex)
        for f, E0, gamma in oscillators:
            eps -= f * E0**2 / (E_eV**2 - E0**2 + 1j*E_eV*gamma)
        
        return np.sqrt(eps)
    
    def get_material_index(self, material, wavelengths):
        """Get refractive index for a specific material"""
        material = material.lower()
        
        if material in ['ag', 'silver']:
            return self.silver_refractive_index(wavelengths)
        elif material in ['au', 'gold']:
            return self.gold_refractive_index(wavelengths)
        elif material in ['hbn', 'h-bn']:
            return self.hbn_refractive_index(wavelengths)
        elif material == 'pmma':
            return self.pmma_refractive_index(wavelengths)
        elif material == 'ws2':
            return self.lorentz_ws2(wavelengths)
        elif material == 'mos2':
            return self.mos2_refractive_index(wavelengths)
        elif material == 'tdbc':
            return self.tdbc_refractive_index(wavelengths)
        elif material == 'sio2':
            return self.sio2_refractive_index(wavelengths)
        elif material in ['silicon', 'si']:
            return self.silicon_refractive_index(wavelengths)
        elif material == 'air':
            return np.ones_like(wavelengths, dtype=complex)
        elif material == 'glass':
            return np.full_like(wavelengths, 1.5, dtype=complex)
        else:
            # Try to load from txt file for new materials
            filename = f"{material}.txt"
            if os.path.exists(filename):
                print(f"Loading refractive index for {material} from {filename}")
                return self.nComplex_interp(filename, wavelengths)
            else:
                print(f"Unknown material: {material}, using n=1")
                return np.ones_like(wavelengths, dtype=complex)