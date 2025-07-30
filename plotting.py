import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

class TMMPlotter:
    """Plotting functions for TMM results"""
    
    def __init__(self):
        self.echelle_colormap = self.create_echelle_colormap()
        # Store colorbars to avoid duplication
        self.colorbars = {}
    
    def create_echelle_colormap(self):
        """Create custom colormap exactly matching MATLAB version"""
        colors = [
            [0, 0, 0], [0.06835, 0.008964, 0.1277], [0.1367, 0.01793, 0.2555],
            [0.205, 0.02689, 0.3832], [0.2734, 0.03585, 0.5109],
            [0.3417, 0.04482, 0.6387], [0.4101, 0.05378, 0.7664],
            [0.4784, 0.06275, 0.8941], [0.4485, 0.1213, 0.9007],
            [0.4186, 0.1799, 0.9074], [0.3887, 0.2385, 0.914],
            [0.3588, 0.2971, 0.9206], [0.3289, 0.3556, 0.9272],
            [0.299, 0.4142, 0.9338], [0.2691, 0.4728, 0.9404], [0.4, 0.4, 1],
            [0.35, 0.475, 1], [0.3, 0.55, 1], [0.25, 0.625, 1], [0.2, 0.7, 1],
            [0.15, 0.775, 1], [0.1, 0.85, 1], [0.05, 0.925, 1], [0, 1, 1],
            [0, 0.9373, 0.875], [0, 0.8745, 0.75], [0, 0.8118, 0.625],
            [0, 0.749, 0.5], [0, 0.6863, 0.375], [0, 0.6235, 0.25],
            [0, 0.5608, 0.125], [0, 0.498, 0], [0.125, 0.5608, 0],
            [0.25, 0.6235, 0], [0.375, 0.6863, 0], [0.5, 0.749, 0],
            [0.625, 0.8118, 0], [0.75, 0.8745, 0], [0.875, 0.9373, 0],
            [1, 1, 0], [1, 0.9377, 0], [1, 0.8755, 0], [1, 0.8132, 0],
            [1, 0.751, 0], [1, 0.6887, 0], [1, 0.6265, 0], [1, 0.5642, 0],
            [1, 0.502, 0], [1, 0.4392, 0], [1, 0.3765, 0], [1, 0.3137, 0],
            [1, 0.251, 0], [1, 0.1882, 0], [1, 0.1255, 0], [1, 0.06275, 0],
            [1, 0, 0], [0.9375, 0, 0], [0.875, 0, 0], [0.8125, 0, 0],
            [0.75, 0, 0], [0.6875, 0, 0], [0.625, 0, 0], [0.5625, 0, 0], [0.5, 0, 0]
        ]
        
        from matplotlib.colors import ListedColormap
        return ListedColormap(colors)
    
    def plot_tra_spectrum(self, fig, ax, wavelengths, transmission, reflection, absorption, show_trans=True, show_refl=True, show_abs=True):
        """Plot transmission, reflection, and absorption spectrum with selectable curves"""
        ax.clear()
        labels = []
        if show_trans:
            ax.plot(wavelengths, transmission, 'b-', label='Transmission', linewidth=2)
            labels.append('Transmission')
        if show_refl:
            ax.plot(wavelengths, reflection, 'r-', label='Reflection', linewidth=2)
            labels.append('Reflection')
        if show_abs:
            ax.plot(wavelengths, absorption, 'g-', label='Absorption', linewidth=2)
            labels.append('Absorption')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('T, R, A')
        ax.set_title('Transmission, Reflection, and Absorption')
        if labels:
            ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 1)
        fig.tight_layout()
    
    def plot_field_intensity(self, fig, ax, wavelengths, z, field, thicknesses):
        """Plot electric field intensity with proper colorbar management"""
        # Clear the entire figure to remove old colorbars
        fig.clear()
        ax = fig.add_subplot(111)
        
        # Remove any stored field colorbars
        field_colorbar_keys = [k for k in self.colorbars.keys() if 'field' in k]
        for key in field_colorbar_keys:
            try:
                self.colorbars[key].remove()
            except:
                pass
            del self.colorbars[key]
        
        # Plot field intensity
        intensity = np.abs(field)**2
        
        extent = [wavelengths[0], wavelengths[-1], z[0], z[-1]]
        im = ax.imshow(intensity.T, aspect='auto', extent=extent, 
                      origin='lower', cmap='hot')
        
        # Add layer boundaries
        cumulative_thickness = np.cumsum([0] + list(thicknesses))
        z_start = z[0]
        
        for i, pos in enumerate(cumulative_thickness[:-1]):
            z_pos = z_start + pos
            ax.axhline(z_pos, color='white', linewidth=1, alpha=0.7)
        
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Position (nm)')
        ax.set_title('Electric Field Intensity |E|²')
        
        # Add colorbar only once
        try:
            self.colorbars['field_new'] = fig.colorbar(im, ax=ax, label='|E|²')
        except Exception as e:
            print(f"Error adding colorbar: {e}")
        
        fig.tight_layout()
        return im
    
    def calculate_dispersion_coordinates(self, wavelengths, angles):
        """Calculate k-vectors and energy coordinates exactly like MATLAB"""
        wavelengths = np.array(wavelengths)
        angles = np.array(angles)
        
        # Convert wavelengths to energy (eV) - exact MATLAB formula
        Eev = 1240.0 / wavelengths  # Energy in eV
        
        # Create coordinate matrices matching MATLAB dimensions exactly
        NstepLambda = len(wavelengths)
        NstepTheta = len(angles)
        
        kxx = np.zeros((NstepLambda, NstepTheta))
        EEV = np.zeros((NstepLambda, NstepTheta))
        
        # Replicate MATLAB nested loops exactly
        for i in range(NstepLambda):
            for j in range(NstepTheta):
                # MATLAB: kxx(i,j) = 2*pi/veclambda(i) * sin(vectheta(j)*pi/180) * 1E3
                kxx[i, j] = 2 * np.pi / wavelengths[i] * np.sin(angles[j] * np.pi / 180) * 1e3
                # MATLAB: EEV(i,j) = Eev(i)
                EEV[i, j] = Eev[i]
        
        return kxx, EEV
    
    def plot_multiple_dispersion(self, fig, wavelengths, angles, transmission, reflection, absorption):
        """Plot multiple dispersion plots exactly matching MATLAB implementation"""
        # Clear the entire figure and all stored colorbars
        fig.clear()
        
        # Clear all stored colorbars for this figure
        keys_to_remove = [k for k in self.colorbars.keys() if 'multi' in k or 'disp' in k]
        for key in keys_to_remove:
            try:
                self.colorbars[key].remove()
            except:
                pass
            del self.colorbars[key]
        
        try:
            # Ensure inputs are numpy arrays with correct shapes
            wavelengths = np.array(wavelengths)
            angles = np.array(angles)
            transmission = np.array(transmission)
            reflection = np.array(reflection)
            absorption = np.array(absorption)
            
            # Calculate coordinates exactly like MATLAB
            kxx, EEV = self.calculate_dispersion_coordinates(wavelengths, angles)
            
            # Verify dimensions match exactly
            print(f"MATLAB replication check:")
            print(f"Wavelengths: {len(wavelengths)}, Angles: {len(angles)}")
            print(f"kxx shape: {kxx.shape}, EEV shape: {EEV.shape}")
            print(f"Transmission shape: {transmission.shape}")
            print(f"Reflection shape: {reflection.shape}")
            print(f"Absorption shape: {absorption.shape}")
            
            # Create subplots matching MATLAB layout exactly (2x2)
            ax1 = fig.add_subplot(2, 2, 1)  # Transmission
            ax2 = fig.add_subplot(2, 2, 2)  # Reflection  
            ax3 = fig.add_subplot(2, 2, 3)  # Absorption
            ax4 = fig.add_subplot(2, 2, 4)  # Combined contours
            
            # 1. Transmission Dispersion - exact MATLAB replication
            im1 = ax1.pcolor(kxx, EEV, transmission, shading='auto', cmap=self.echelle_colormap)
            ax1.set_xlabel('In-plane momentum [μm⁻¹]')
            ax1.set_ylabel('Energy [eV]')
            ax1.set_title('Transmission Dispersion')
            ax1.set_ylim([1.9, 2.15])  # Exact MATLAB ylim
            ax1.grid(True, alpha=0.3)
            im1.set_clim([0, 1])  # Exact MATLAB caxis
            
            # Add colorbar for transmission
            try:
                cb1 = fig.colorbar(im1, ax=ax1, shrink=0.8)
                cb1.set_label('Transmission')
                self.colorbars['multi_trans_new'] = cb1
            except Exception as e:
                print(f"Error adding transmission colorbar: {e}")
            
            # 2. Reflection Dispersion - exact MATLAB replication
            im2 = ax2.pcolor(kxx, EEV, reflection, shading='auto', cmap=self.echelle_colormap)
            ax2.set_xlabel('In-plane momentum [μm⁻¹]')
            ax2.set_ylabel('Energy [eV]')
            ax2.set_title('Reflection Dispersion')
            ax2.set_ylim([1.9, 2.15])  # Exact MATLAB ylim
            ax2.grid(True, alpha=0.3)
            im2.set_clim([0, 1])  # Exact MATLAB caxis
            
            # Add colorbar for reflection
            try:
                cb2 = fig.colorbar(im2, ax=ax2, shrink=0.8)
                cb2.set_label('Reflection')
                self.colorbars['multi_refl_new'] = cb2
            except Exception as e:
                print(f"Error adding reflection colorbar: {e}")
            
            # 3. Absorption Dispersion - exact MATLAB replication
            im3 = ax3.pcolor(kxx, EEV, absorption, shading='auto', cmap=self.echelle_colormap)
            ax3.set_xlabel('In-plane momentum [μm⁻¹]')
            ax3.set_ylabel('Energy [eV]')
            ax3.set_title('Absorption Dispersion')
            ax3.set_ylim([1.9, 2.15])  # Exact MATLAB ylim
            ax3.grid(True, alpha=0.3)
            im3.set_clim([0, 1])  # Exact MATLAB caxis
            
            # Add colorbar for absorption
            try:
                cb3 = fig.colorbar(im3, ax=ax3, shrink=0.8)
                cb3.set_label('Absorption')
                self.colorbars['multi_abs_new'] = cb3
            except Exception as e:
                print(f"Error adding absorption colorbar: {e}")
            
            # 4. Combined Contour Plot - exact MATLAB replication
            # MATLAB: contour(kxx, EEV, Reflection, [0.1:0.1:0.9], 'r', 'LineWidth', 1);
            contour_levels = np.arange(0.1, 1.0, 0.1)  # [0.1:0.1:0.9]
            
            # Reflection contours (red)
            cs1 = ax4.contour(kxx, EEV, reflection, levels=contour_levels, 
                             colors='r', linewidths=1)
            
            # Transmission contours (blue) 
            cs2 = ax4.contour(kxx, EEV, transmission, levels=contour_levels,
                             colors='b', linewidths=1)
            
            # Absorption contours (green)
            cs3 = ax4.contour(kxx, EEV, absorption, levels=contour_levels,
                             colors='g', linewidths=1)
            
            ax4.set_xlabel('In-plane momentum [μm⁻¹]')
            ax4.set_ylabel('Energy [eV]')
            ax4.set_title('Combined Contour Plot')
            ax4.set_ylim([1.9, 2.15])  # Exact MATLAB ylim
            ax4.grid(True, alpha=0.3)
            
            # Create legend for contour plot
            import matplotlib.lines as mlines
            red_line = mlines.Line2D([], [], color='r', label='Reflection')
            blue_line = mlines.Line2D([], [], color='b', label='Transmission')  
            green_line = mlines.Line2D([], [], color='g', label='Absorption')
            ax4.legend(handles=[red_line, blue_line, green_line], loc='best')
            
            # Add main title exactly like MATLAB
            # MATLAB: sgtitle('Dispersion Diagrams (Ag/hBN/WS2/hBN/Ag)');
            fig.suptitle('Dispersion Diagrams', fontsize=14, y=0.98)
            
            # Adjust layout to prevent overlap
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            
            return [im1, im2, im3]
            
        except Exception as e:
            # Error handling
            print(f"Error in dispersion plotting: {e}")
            ax_err = fig.add_subplot(1, 1, 1)
            ax_err.text(0.5, 0.5, f'Error plotting multiple dispersions:\n{str(e)}', 
                       transform=ax_err.transAxes, ha='center', va='center')
            fig.tight_layout()
            return []
    
    def clear_all_colorbars(self):
        """Clear all stored colorbars"""
        for key in list(self.colorbars.keys()):
            try:
                self.colorbars[key].remove()
            except:
                pass
            del self.colorbars[key]