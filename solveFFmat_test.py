from scipy.sparse.linalg import LinearOperator, gmres
from inputParams.simulation_parameters import params_simu
from readRWGfiles import *
import pickle
from MLFMA import computeCurrentsVisualization

mesh_path = "./mesh/"
result_path = "./result/"

N_RWG = readIntFromDisk(mesh_path + "N_RWG.txt")
r_obs = np.load(mesh_path + "r_obs.npy")
Vobs = np.load("mesh/vObs.npy")
Nobs = Vobs.size

S = readIntFromDisk(mesh_path + "S.txt")

FF_J_th, FF_J_ph, FF_M_th, FF_M_ph \
    = np.load("./mat/FF_J_th.npy"), np.load("./mat/FF_J_ph.npy"), \
      np.load("./mat/FF_M_th.npy"), np.load("./mat/FF_M_ph.npy")
alpha_th, alpha_ph = np.load("./mat/alpha_th.npy"), np.load("./mat/alpha_ph.npy")
matvec = lambda vec: np.dot(alpha_th, np.dot(FF_J_th, vec)) \
                     + np.dot(alpha_ph, np.dot(FF_J_ph, vec))

with open("mat/surface_RWG.pickle", "rb") as f:
    surface_RWGnumbers, surface_NRWG = pickle.load(f)

Z_PMCHWT_i = np.load("mat/Z_PMCHWT_i.npy")
Z_PMCHWT_i_r = Z_PMCHWT_i.conj().T

Dsurface = params_simu.dielectricSurface
Ssurface = 1 - Dsurface

if params_simu.CurrentModel == "SEqF_J":

    Z_CFIE = -np.r_[np.load("mat/Z_CFIE_EJ_f.npy")[surface_RWGnumbers[Dsurface][:, np.newaxis],
                                                   surface_RWGnumbers[Ssurface]],
                    np.load("mat/Z_CFIE_HJ_f.npy")[surface_RWGnumbers[Dsurface][:, np.newaxis],
                                                   surface_RWGnumbers[Ssurface]]]
    Z_CFIE_r = Z_CFIE.conj().T

    FF = []
    FF.append(np.r_[FF_J_th[:, surface_RWGnumbers[Ssurface]],
                    FF_J_ph[:, surface_RWGnumbers[Ssurface]]])
    FF.append(np.c_[np.r_[FF_J_th[:, surface_RWGnumbers[Dsurface]],
                          FF_J_ph[:, surface_RWGnumbers[Dsurface]]],
                    np.r_[FF_M_th[:, surface_RWGnumbers[Dsurface]],
                          FF_M_ph[:, surface_RWGnumbers[Dsurface]]]])
    FF_r = []
    FF_r.append(FF[0].conj().T)
    FF_r.append(FF[1].conj().T)

    alpha = np.c_[alpha_th, alpha_ph]
    alpha_r = alpha.conj().T


    def matvec(vec):
        ff1 = np.dot(FF[1], np.dot(Z_PMCHWT_i, np.dot(Z_CFIE, vec)))
        return np.dot(alpha, np.dot(FF[0], vec) + ff1)


    def rmatvec(vec):
        arvec = np.dot(alpha_r, vec)
        ff1 = np.dot(Z_CFIE_r, np.dot(Z_PMCHWT_i_r, np.dot(FF_r[1], arvec)))
        return np.dot(FF_r[0], arvec) + ff1

# x_test = np.zeros(surface_NRWG[Ssurface], dtype=complex)
x_test = np.zeros(N_RWG, dtype=complex)
x_test[0] = 1.
FF_J = np.r_[FF_J_th, FF_J_ph]
y_test = np.dot(alpha, np.dot(FF_J, x_test))

from Z_obs import Pywrap_compute_Z_obs
mesh_data = MeshData(mesh_path)
angular_freq = 2 * np.pi * params_simu.f
Z_obs = Pywrap_compute_Z_obs(
    r_obs, mesh_data.RWGNumber_trianglesCoord, mesh_data.triangle_NRWG,
    mesh_data.triangle_RWGNumber, mesh_data.triangle_signInRWG, mesh_data.triangle_surface,
    angular_freq, params_simu.eps_r_f, params_simu.mu_r_f, np.asarray([0, 1], dtype="i")
)
y_test_2 = Z_obs[0].dot(x_test)[:, 0]
