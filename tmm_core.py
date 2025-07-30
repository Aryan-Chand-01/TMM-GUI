import numpy as np
from material_properties import MaterialProperties

class TMMCore:
    """Transfer Matrix Method calculation core"""
    
    def __init__(self):
        self.c = 3e8  # Speed of light
        self.materials = MaterialProperties()
    
    def get_refractive_indices(self, wavelengths, materials, substrate):
        """Get refractive index matrix for all materials and wavelengths"""
        wavelengths = np.array(wavelengths)
        n_wl = len(wavelengths)
        n_layers = len(materials) + 2  # Include input and output media
        
        n_matrix = np.ones((n_layers, n_wl), dtype=complex)
        
        # Input medium (air)
        n_matrix[0, :] = 1.0
        
        # Each layer
        for i, material in enumerate(materials):
            n_matrix[i + 1, :] = self.materials.get_material_index(material, wavelengths)
        
        # Output medium (substrate)
        n_matrix[-1, :] = self.materials.get_material_index(substrate, wavelengths)
        
        return n_matrix
    
    def calculate_transfer_matrices(self, indices, theta, omega0, thicknesses, polarization):
        """Calculate transfer matrices for given parameters"""
        N = len(thicknesses)
        
        # Initialize angles and wavevectors
        theta_array = np.zeros(N + 2, dtype=complex)
        kz = np.zeros(N + 2, dtype=complex)
        theta_array[0] = theta
        kz[0] = indices[0] * omega0 * np.cos(theta_array[0]) / self.c
        
        # Transfer matrices
        T_ml_TE = np.eye(2, dtype=complex)
        T_ml_TM = np.eye(2, dtype=complex)
        
        # Store matrices for field calculation
        T_int_TE = []
        T_int_TM = []
        T_prop = []
        
        # Calculate through each layer
        for layer_idx in range(N):
            # Snell's law
            sin_theta_next = (indices[layer_idx] * np.sin(theta_array[layer_idx]) / 
                            indices[layer_idx + 1])
            
            # Handle total internal reflection
            if np.abs(sin_theta_next) > 1:
                theta_array[layer_idx + 1] = (np.pi/2 - 
                                            1j * np.arccosh(np.abs(sin_theta_next)))
            else:
                theta_array[layer_idx + 1] = np.arcsin(sin_theta_next)
            
            kz[layer_idx + 1] = (indices[layer_idx + 1] * omega0 * 
                               np.cos(theta_array[layer_idx + 1]) / self.c)
            
            # Ensure correct sign for imaginary part (absorbing media)
            if np.imag(kz[layer_idx + 1]) < 0:
                kz[layer_idx + 1] = -np.conj(kz[layer_idx + 1])
            
            # Interface matrix coefficients
            kappa = kz[layer_idx + 1] / kz[layer_idx]
            eta = indices[layer_idx + 1] / indices[layer_idx]
            
            # Interface matrices
            T_int_TM_current = 0.5 * np.array([
                [eta + kappa/eta, -eta + kappa/eta],
                [-eta + kappa/eta, eta + kappa/eta]
            ])
            T_int_TE_current = 0.5 * np.array([
                [1 + kappa, 1 - kappa],
                [1 - kappa, 1 + kappa]
            ])
            
            # Propagation matrix
            phase = kz[layer_idx + 1] * thicknesses[layer_idx]
            T_prop_current = np.array([
                [np.exp(-1j * phase), 0],
                [0, np.exp(1j * phase)]
            ])
            
            T_int_TE.append(T_int_TE_current)
            T_int_TM.append(T_int_TM_current)
            T_prop.append(T_prop_current)
            
            # Update total transfer matrix
            T_ml_TM = T_ml_TM @ T_int_TM_current @ T_prop_current
            T_ml_TE = T_ml_TE @ T_int_TE_current @ T_prop_current
        
        # Final interface (to substrate)
        sin_theta_final = (indices[N] * np.sin(theta_array[N]) / indices[N + 1])
        if np.abs(sin_theta_final) > 1:
            theta_array[N + 1] = np.pi/2 - 1j * np.arccosh(np.abs(sin_theta_final))
        else:
            theta_array[N + 1] = np.arcsin(sin_theta_final)
        
        kz[N + 1] = (indices[N + 1] * omega0 * np.cos(theta_array[N + 1]) / self.c)
        
        if np.imag(kz[N + 1]) < 0:
            kz[N + 1] = -np.conj(kz[N + 1])
        
        kappa = kz[N + 1] / kz[N]
        eta = indices[N + 1] / indices[N]
        
        T_int_TM_final = 0.5 * np.array([
            [eta + kappa/eta, -eta + kappa/eta],
            [-eta + kappa/eta, eta + kappa/eta]
        ])
        T_int_TE_final = 0.5 * np.array([
            [1 + kappa, 1 - kappa],
            [1 - kappa, 1 + kappa]
        ])
        
        T_int_TE.append(T_int_TE_final)
        T_int_TM.append(T_int_TM_final)
        
        T_ml_TM = T_ml_TM @ T_int_TM_final
        T_ml_TE = T_ml_TE @ T_int_TE_final
        
        return {
            'T_ml_TM': T_ml_TM,
            'T_ml_TE': T_ml_TE,
            'T_int_TM': T_int_TM,
            'T_int_TE': T_int_TE,
            'T_prop': T_prop,
            'kz': kz,
            'theta': theta_array
        }
    
    def calculate_coefficients(self, transfer_matrices, polarization):
        """Calculate transmission and reflection coefficients"""
        T_ml_TM = transfer_matrices['T_ml_TM']
        T_ml_TE = transfer_matrices['T_ml_TE']
        kz = transfer_matrices['kz']
        
        # Calculate transmission and reflection coefficients
        t_TM = 1 / T_ml_TM[0, 0]
        r_TM = T_ml_TM[1, 0] / T_ml_TM[0, 0]
        t_TE = 1 / T_ml_TE[0, 0]
        r_TE = T_ml_TE[1, 0] / T_ml_TE[0, 0]
        
        # Power transmission and reflection
        cos_factor = np.real(kz[-1]) / np.real(kz[0])
        
        if polarization == 'TM':
            T_power = cos_factor * np.abs(t_TM)**2
            R_power = np.abs(r_TM)**2
            t_coeff = t_TM
            r_coeff = r_TM
        elif polarization == 'TE':
            T_power = cos_factor * np.abs(t_TE)**2
            R_power = np.abs(r_TE)**2
            t_coeff = t_TE
            r_coeff = r_TE
        else:  # Unpolarized
            T_power = 0.5 * cos_factor * (np.abs(t_TM)**2 + np.abs(t_TE)**2)
            R_power = 0.5 * (np.abs(r_TE)**2 + np.abs(r_TM)**2)
            t_coeff = t_TE  # Use TE for field calculation
            r_coeff = r_TE
        
        return {
            'transmission': np.real(T_power),
            'reflection': np.real(R_power),
            't_coeff': t_coeff,
            'r_coeff': r_coeff
        }
    
    def calculate_field(self, transfer_matrices, coeffs, thicknesses, z):
        """Calculate electric field distribution"""
        kz = transfer_matrices['kz']
        T_int_TE = transfer_matrices['T_int_TE']
        T_prop = transfer_matrices['T_prop']
        r = coeffs['r_coeff']
        t = coeffs['t_coeff']
        
        N = len(thicknesses)
        Nz = len(z)
        
        # Field amplitudes: [forward, backward] for each layer
        A = np.zeros((2, N+2), dtype=complex)
        Ef = np.zeros((N+2, Nz), dtype=complex)  # Forward field
        Eb = np.zeros((N+2, Nz), dtype=complex)  # Backward field
        
        # Cumulative layer positions
        cumL = np.cumsum(thicknesses)
        
        # Left side (incident medium)
        A[:, 0] = np.array([1, r])
        mask_left = z < 0
        Ef[0, mask_left] = A[0, 0]
        Eb[0, mask_left] = A[1, 0]
        
        # Right side (substrate)
        A[:, N+1] = np.array([t, 0])
        mask_right = z > np.sum(thicknesses)
        Ef[N+1, mask_right] = A[0, N+1]
        Eb[N+1, mask_right] = A[1, N+1]
        
        # Calculate field amplitudes in each layer (backwards)
        for idx in range(N, 0, -1):  # N down to 1
            A[:, idx] = T_prop[idx-1] @ T_int_TE[idx] @ A[:, idx+1]
            
            if idx > 1:
                z_start = cumL[idx-2]
                z_end = cumL[idx-1]
                mask = (z >= z_start) & (z < z_end)
                phase_ref = z_start
            else:
                z_start = 0
                z_end = thicknesses[0]
                mask = (z >= 0) & (z < thicknesses[0])
                phase_ref = 0
            
            if np.any(mask):
                phase_forward = 1j * kz[idx] * (z[mask] - phase_ref)
                phase_backward = -1j * kz[idx] * (z[mask] - phase_ref)
                
                Ef[idx, mask] = A[0, idx] * np.exp(phase_forward)
                Eb[idx, mask] = A[1, idx] * np.exp(phase_backward)
        
        # Total field
        E_total = np.sum(Ef + Eb, axis=0)
        
        return E_total
    
    def tmm_calculation(self, wavelengths, thicknesses, materials, substrate='air', 
                       polarization='', angles=None, calc_field=False):
        """
        Main TMM calculation with field computation option
        """
        if angles is None:
            angles = np.array([0])
        
        wavelengths = np.array(wavelengths)
        thicknesses = np.array(thicknesses) * 1e-9  # Convert to meters
        N = len(thicknesses)
        
        # Build refractive index matrix
        n_layers = self.get_refractive_indices(wavelengths, materials, substrate)
        
        # Initialize results
        transmission = np.zeros((len(wavelengths), len(angles)))
        reflection = np.zeros((len(wavelengths), len(angles)))
        
        # Field calculation setup
        if calc_field:
            z_min = -50e-9
            z_max = np.sum(thicknesses) + 50e-9
            Nz = 1024
            z = np.linspace(z_min, z_max, Nz)
            E_field = np.zeros((len(wavelengths), Nz), dtype=complex)
        else:
            E_field = None
            z = None
        
        for angle_idx, theta_deg in enumerate(angles):
            theta_rad = np.deg2rad(theta_deg)
            
            for wl_idx, wl in enumerate(wavelengths):
                # Current wavelength parameters
                lambda0 = wl * 1e-9  # Convert to meters
                k0 = 2 * np.pi / lambda0
                omega0 = self.c * k0
                
                # Refractive indices for current wavelength
                indices = n_layers[:, wl_idx]
                
                # Calculate transfer matrices
                transfer_matrices = self.calculate_transfer_matrices(
                    indices, theta_rad, omega0, thicknesses, polarization
                )
                
                # Calculate coefficients
                coeffs = self.calculate_coefficients(transfer_matrices, polarization)
                
                transmission[wl_idx, angle_idx] = coeffs['transmission']
                reflection[wl_idx, angle_idx] = coeffs['reflection']
                
                # Field calculation for normal incidence only
                if calc_field and angle_idx == 0:
                    E_field[wl_idx, :] = self.calculate_field(
                        transfer_matrices, coeffs, thicknesses, z
                    )
        
        absorption = 1 - transmission - reflection
        
        result = {
            'transmission': transmission,
            'reflection': reflection,
            'absorption': absorption
        }
        
        if calc_field:
            result['field'] = E_field
            result['z'] = z
        
        return result