import numpy as np
from Z_obs import Pywrap_compute_Z_obs
from E_FF import Pywrap_compute_FFmat


class ObserveField:
    def __init__(self, mesh, coefficients_j, coefficients_m, frequency):
        self.meshData = mesh
        self.coefficients_J = coefficients_j
        self.coefficients_M = coefficients_m
        self.w = 2 * np.pi * frequency

    def computeFields(self, r_obs, eps_r, mu_r, targetSurfaces):
        z_ej, z_em, z_hj, z_hm \
            = Pywrap_compute_Z_obs(r_obs, self.meshData.RWGNumber_trianglesCoord,
                                   self.meshData.triangle_NRWG,
                                   self.meshData.triangle_RWGNumber,
                                   self.meshData.triangle_signInRWG,
                                   self.meshData.triangle_surface, self.w,
                                   eps_r, mu_r, targetSurfaces)
        e_obs = np.dot(z_ej, self.coefficients_J) + np.dot(z_em, self.coefficients_M)
        h_obs = np.dot(z_hj, self.coefficients_J) + np.dot(z_hm, self.coefficients_M)
        return e_obs, h_obs

    def computeFields_sphericalCoord(self, r, theta, phi, eps_r, mu_r, target_surfaces):
        r_cartesian = np.c_[
            r * np.sin(theta) * np.cos(phi), r * np.sin(theta) * np.sin(phi),
            r * np.cos(theta)
        ]
        e_cartesian, h_cartesian = self.computeFields(r_cartesian,
                                                      eps_r, mu_r, target_surfaces)

        e_obs_r = np.zeros(theta.size, dtype=complex)  # not implemented now
        e_obs_th = np.c_[np.cos(theta) * np.cos(phi) * e_cartesian[:, 0]
                         + np.cos(theta) * np.sin(phi) * e_cartesian[:, 1]
                         - np.sin(theta) * e_cartesian[:, 2]]
        e_obs_ph = np.c_[-np.sin(phi) * e_cartesian[:, 0] + np.cos(phi) * e_cartesian[:, 1]]

        h_obs_r = np.zeros(theta.size, dtype=complex)  # not implemented now
        h_obs_th = np.c_[np.cos(theta) * np.cos(phi) * h_cartesian[:, 0]
                         + np.cos(theta) * np.sin(phi) * h_cartesian[:, 1]
                         - np.sin(theta) * h_cartesian[:, 2]]
        h_obs_ph = np.c_[-np.sin(phi) * h_cartesian[:, 0] + np.cos(phi) * h_cartesian[:, 1]]

        return np.c_[e_obs_r, e_obs_th, e_obs_ph], np.c_[h_obs_r, h_obs_th, h_obs_ph]

    # def far_field(self, theta, phi, eps_r, mu_r, target_surface):
    #     ff = Pywrap_compute_FF(theta, phi,
    #                            self.meshData.RWGNumber_trianglesCoord,
    #                            self.coefficients_J, self.coefficients_M,
    #                            self.meshData.triangle_NRWG,
    #                            self.meshData.triangle_RWGNumber,
    #                            self.meshData.triangle_signInRWG,
    #                            self.meshData.triangle_surface, self.w, eps_r, mu_r,
    #                            target_surface)
    #     return ff

    def far_field_matrix(self, theta, phi, eps_r, mu_r):
        matrix = Pywrap_compute_FFmat(theta, phi,
                                      self.meshData.RWGNumber_trianglesCoord,
                                      self.meshData.triangle_NRWG,
                                      self.meshData.triangle_RWGNumber,
                                      self.meshData.triangle_signInRWG,
                                      self.meshData.triangle_surface, self.w, eps_r, mu_r)
        return matrix


def saveObsFile(r, f, filename):
    data_to_save = np.c_[r, f.real[:, 0], f.imag[:, 0], f.real[:, 1], f.imag[:, 1],
                         f.real[:, 2], f.imag[:, 2]]
    header = "r_obs_x[m] r_obs_y[m] r_obs_z[m] " \
             + "re(Ex)[V/m] im(Ex)[V/m] re(Ey)[V/m] im(Ey)[V/m] re(Ez)[V/m] im(Ez)[V/m]"
    np.savetxt(filename, data_to_save, header=header)
