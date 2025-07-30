import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d

class TMMAnalysis:
    """Analysis functions for TMM results"""
    
    def __init__(self):
        pass
    
    def find_rabi_splitting(self, wavelengths, reflection):
        """Find Rabi splitting from reflection spectrum"""
        # Ensure 1D data
        reflection_1d = reflection.flatten() if hasattr(reflection, 'flatten') else reflection
        wavelengths = np.array(wavelengths).flatten()
        
        # Find minima in reflection
        try:
            # Smooth the data slightly to avoid noise
            smoothed = gaussian_filter1d(reflection_1d, sigma=1)
            
            # Find peaks (minima in reflection = peaks in -reflection)
            peaks, properties = find_peaks(-smoothed, height=-0.8, distance=10)
            
            if len(peaks) >= 2:
                # Sort peaks by depth (most negative reflection)
                peak_depths = smoothed[peaks]
                sorted_indices = np.argsort(peak_depths)
                
                # Take the two deepest peaks
                peak1_idx = peaks[sorted_indices[0]]
                peak2_idx = peaks[sorted_indices[1]]
                
                # Sort by wavelength
                if peak1_idx > peak2_idx:
                    peak1_idx, peak2_idx = peak2_idx, peak1_idx
                
                wl1 = wavelengths[peak1_idx]
                wl2 = wavelengths[peak2_idx]
                
                # Convert to energy
                E1 = 1240.0 / wl1  # meV
                E2 = 1240.0 / wl2  # meV
                
                rabi_splitting_wl = abs(wl2 - wl1)
                rabi_splitting_energy = abs(E2 - E1)
                
                return {
                    'splitting_nm': rabi_splitting_wl,
                    'splitting_meV': rabi_splitting_energy,
                    'peak1_nm': wl1,
                    'peak2_nm': wl2,
                    'peak1_meV': E1,
                    'peak2_meV': E2,
                    'found': True
                }
        except Exception as e:
            print(f"Error in Rabi splitting analysis: {e}")
        
        return {'found': False}
    
    def calculate_fwhm(self, wavelengths, reflection):
        """Calculate Full Width at Half Maximum"""
        try:
            reflection_1d = reflection.flatten() if hasattr(reflection, 'flatten') else reflection
            wavelengths = np.array(wavelengths).flatten()
            
            # Find minimum
            min_idx = np.argmin(reflection_1d)
            min_wavelength = wavelengths[min_idx]
            min_reflection = reflection_1d[min_idx]
            max_val = np.max(reflection_1d)
            half_max = (max_val + min_reflection) / 2
            
            # Find points where reflection crosses half maximum
            indices = np.where(reflection_1d <= half_max)[0]
            
            if len(indices) >= 2:
                # Find the range around the minimum
                left_idx = indices[indices < min_idx]
                right_idx = indices[indices > min_idx]
                
                if len(left_idx) > 0 and len(right_idx) > 0:
                    left_wl = wavelengths[left_idx[-1]]
                    right_wl = wavelengths[right_idx[0]]
                    fwhm_wl = right_wl - left_wl
                    fwhm_energy = 1240.0 / left_wl - 1240.0 / right_wl
                    
                    # Calculate Q-factor
                    q_factor = min_wavelength / fwhm_wl if fwhm_wl > 0 else 0
                    
                    return {
                        'found': True,
                        'fwhm_nm': fwhm_wl,
                        'fwhm_meV': fwhm_energy,
                        'center_nm': min_wavelength,
                        'q_factor': q_factor
                    }
        except Exception as e:
            print(f"Error in FWHM calculation: {e}")
        
        return {'found': False}
    
    def find_peaks_and_dips(self, wavelengths, spectrum):
        """Find peaks and dips in spectrum"""
        try:
            spectrum_1d = spectrum.flatten() if hasattr(spectrum, 'flatten') else spectrum
            wavelengths = np.array(wavelengths).flatten()
            
            # Find peaks
            peaks, _ = find_peaks(spectrum_1d, height=0.1, distance=5)
            
            # Find dips (peaks in inverted spectrum)
            dips, _ = find_peaks(-spectrum_1d, height=-0.9, distance=5)
            
            peak_wavelengths = wavelengths[peaks] if len(peaks) > 0 else []
            dip_wavelengths = wavelengths[dips] if len(dips) > 0 else []
            
            return {
                'peaks': peak_wavelengths,
                'dips': dip_wavelengths,
                'peak_values': spectrum_1d[peaks] if len(peaks) > 0 else [],
                'dip_values': spectrum_1d[dips] if len(dips) > 0 else []
            }
        except Exception as e:
            print(f"Error in peak/dip finding: {e}")
            return {'peaks': [], 'dips': [], 'peak_values': [], 'dip_values': []}
    
    def calculate_quality_factor(self, wavelengths, reflection):
        """Calculate quality factor from reflection spectrum"""
        fwhm_result = self.calculate_fwhm(wavelengths, reflection)
        if fwhm_result['found']:
            return fwhm_result['q_factor']
        return None
    
    def analyze_spectrum(self, wavelengths, transmission, reflection, absorption):
        """Comprehensive spectrum analysis"""
        results = {}
        
        try:
            # Ensure all inputs are numpy arrays and 1D
            wavelengths = np.array(wavelengths).flatten()
            transmission = np.array(transmission).flatten()
            reflection = np.array(reflection).flatten()
            absorption = np.array(absorption).flatten()
            
            # Basic statistics
            results['max_transmission'] = np.max(transmission)
            results['min_reflection'] = np.min(reflection)
            results['max_absorption'] = np.max(absorption)
            
            # Find resonances
            results['rabi_splitting'] = self.find_rabi_splitting(wavelengths, reflection)
            results['fwhm'] = self.calculate_fwhm(wavelengths, reflection)
            results['peaks_dips'] = self.find_peaks_and_dips(wavelengths, reflection)
            
            # Quality factor
            results['q_factor'] = self.calculate_quality_factor(wavelengths, reflection)
            
        except Exception as e:
            print(f"Error in spectrum analysis: {e}")
            # Return safe defaults
            results = {
                'max_transmission': 0,
                'min_reflection': 0,
                'max_absorption': 0,
                'rabi_splitting': {'found': False},
                'fwhm': {'found': False},
                'peaks_dips': {'peaks': [], 'dips': [], 'peak_values': [], 'dip_values': []},
                'q_factor': None
            }
        
        return results