import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import dblquad


class elTCSF:
    # Temporal Contrast Sensitivity Function based on temporal frequency, eccentricity, and luminance (no area)
    def TCSF(self, t_frequency, eccentricity, luminance):
        '''
        :param t_frequency: Temporal Frequency (Hz)
        :param eccentricity: Eccentricity (degree)
        :param luminance: Luminance (cd/m^2)
        :return: Return contrast sensitivity for the given set of parameters.
        !!! Please avoid using elTCSF independently, as it is an integral component of elaTCSF. Isolated use of elTCSF can result in significant coefficient deviations.
        '''
        # Parameters for the TCSF model (based on the structure of TCSF_IDMS)
        TCSF_n1 = 15
        TCSF_n2 = 16
        TCSF_xi = 154.133
        TCSF_tau = 0.00292069
        TCSF_kappa = 2.12547
        TCSF_zeta = 0.721095
        tcsf_lum_k1 = 0.222269
        tcsf_ecc_k1 = 0.0341811
        tcsf_lum_b1 = 0.6678
        ecc_peak_f = -2 # Frequency shift for eccentricity adjustment

        # Adjust temporal frequency based on eccentricity and luminance
        t_frequency = (t_frequency - ecc_peak_f) / ((1 + tcsf_ecc_k1 * eccentricity)) + ecc_peak_f
        t_frequency = t_frequency / (tcsf_lum_b1 + tcsf_lum_k1 * np.log10(luminance))

        # Compute sensitivity using the TCSF_IDMS model formula
        S = np.abs(TCSF_xi * ((1 + 2j * np.pi * t_frequency * TCSF_tau) ** (-TCSF_n1) -
                              TCSF_zeta * (1 + 2j * np.pi * t_frequency * TCSF_kappa * TCSF_tau) ** (-TCSF_n2)))
        return S

    def S_ecc(self, eccentricity):
        # Compute the eccentricity scaling factor
        ecc_k1 = 0.0330933
        S_ecc_factor = 10 ** (-ecc_k1 * eccentricity)
        return S_ecc_factor

    def S_lum(self, luminance):
        # Compute the luminance scaling factor
        lum_k = [1.76801, 1.62402, 0.533781]
        S_luminance_factor = lum_k[0] * ((1 + lum_k[1] / luminance) ** (-lum_k[2]))
        return S_luminance_factor

    def sensitivity(self, eccentricity, luminance, t_frequency):
        # Compute the overall sensitivity combining eccentricity, luminance, and TCSF
        S = self.S_ecc(eccentricity) * self.S_lum(luminance) * self.TCSF(t_frequency, eccentricity, luminance)
        return S


class elaTCSF:
    # Temporal Contrast Sensitivity Function based on temporal frequency, eccentricity, luminance, and area.
    def __init__(self):
        self.elTCSF = elTCSF()  # Instantiate the elTCSF model
        self.beta = 3.80022 # Nonlinearity parameter for spatial summation
        self.E_thr = 6.52801 # Threshold energy parameter

    def S_CSF(self, csf_model, eccentricity, luminance, t_frequency):
        # Wrapper for the sensitivity calculation from the base TCSF model
        S = csf_model.sensitivity(eccentricity, luminance, t_frequency)
        return S

    def spatial_summation_disk(self, eccentricity, luminance, radius, t_frequency):
        '''
        Compute sensitivity for a circular area using spatial summation
        :param eccentricity: Eccentricity (degree)
        :param luminance: Luminance (cd/m^2)
        :param radius: Radius of the Disk (degree)
        :param t_frequency: Temporal Frequency (Hz)
        :return: Return contrast sensitivity for the given set of parameters.
        '''
        # Integration function for polar coordinates
        S_intergration = lambda r, theta: (self.S_CSF(csf_model=self.elTCSF,
                                                      eccentricity=(r ** 2 + eccentricity ** 2 + 2 * eccentricity * r * np.cos(theta)) ** 0.5,
                                                      luminance=luminance, t_frequency=t_frequency) ** self.beta) * r
        # Perform double integration over the disk
        intergration_value, error = dblquad(S_intergration,  0, 2*np.pi, 0, radius)
        # Compute contrast sensitivity from the integrated energy
        contrast = (self.E_thr / intergration_value) ** (1 / self.beta)
        return 1 / contrast

    def spatial_summation_rectangle(self, eccentricity, luminance, width, height, t_frequency):
        '''
        Compute sensitivity for a rectangular area using spatial summation
        :param eccentricity: Eccentricity (degree)
        :param luminance: Luminance (cd/m^2)
        :param width: Width of the Rectangle
        :param height: Height of the Rectangle
        :param t_frequency: Temporal Frequency (Hz)
        :return: Return contrast sensitivity for the given set of parameters.
        '''
        # Integration function for rectangular coordinates
        S_intergration = lambda degree_x, degree_y: self.S_CSF(csf_model=self.elTCSF,
                                                               eccentricity=(degree_x**2 + degree_y**2)**0.5,
                                                               luminance=luminance, t_frequency=t_frequency) ** self.beta
        # Perform double integration over the rectangle
        intergration_value, error = dblquad(S_intergration, -height / 2, height / 2, eccentricity - width / 2, eccentricity + width / 2)
        # Compute contrast sensitivity from the integrated energy
        contrast = (self.E_thr / intergration_value) ** (1 / self.beta)
        return 1 / contrast

    def sensitivity_disk(self, eccentricity, luminance, radius, t_frequency):
        # Compute sensitivity for a circular stimulus
        S = self.spatial_summation_disk(eccentricity, luminance, radius, t_frequency)
        return S

    def sensitivity_rectangle(self, eccentricity, luminance, width, height, t_frequency):
        # Compute sensitivity for a rectangular stimulus
        S = self.spatial_summation_rectangle(eccentricity, luminance, width, height, t_frequency)
        return S


if __name__ == '__main__':
    # Instantiate the extended TCSF model
    CSF_elaTCSF = elaTCSF()
    # Compute sensitivity for a circular area with a given radius
    S_disk_radius_16 = CSF_elaTCSF.sensitivity_disk(eccentricity=10, luminance=3, radius=16, t_frequency=10)
    # Compute sensitivity for a rectangular area with specified width and height
    S_rectangle_width_40_height_32 = CSF_elaTCSF.sensitivity_rectangle(eccentricity=10, luminance=3, width=40, height=32, t_frequency=10)

