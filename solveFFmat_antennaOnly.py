import numpy as np
from scipy.sparse.linalg import lsqr, LinearOperator, gmres, bicgstab, aslinearoperator
from MLFMA import computeCurrentsVisualization
import sys, time
from inputParams.simulation_parameters import params_simu
from ReadWriteBlitzArray import readIntFromDisk
from readRWGfiles import *

mesh_path = "./mesh/"
result_path = "./result/"

N_RWG = readIntFromDisk(mesh_path + "N_RWG.txt")
r_obs = np.load(mesh_path + "r_obs.npy")
Vobs = np.load("mesh/E_obs.npy")
Nobs = Vobs.size

S = readIntFromDisk(mesh_path + "S.txt")

FF_J_th, FF_J_ph, FF_M_th, FF_M_ph \
    = np.load("./mat/FF_J_th.npy"), np.load("./mat/FF_J_ph.npy"), \
      np.load("./mat/FF_M_th.npy"), np.load("./mat/FF_M_ph.npy")
alpha_th, alpha_ph = np.load("./mat/alpha_th.npy"), np.load("./mat/alpha_ph.npy")
matvec = lambda vec: np.dot(alpha_th, np.dot(FF_J_th, vec)) \
                     + np.dot(alpha_ph, np.dot(FF_J_ph, vec))

import pickle

with open("mat/Z_obs.pickle", "r") as f:
    Z_obs = pickle.load(f)
with open("mat/surface_RWG.pickle", "r") as f:
    surface_RWGnumbers, surface_NRWG = pickle.load(f)

Z_PMCHWT_i = np.load("mat/Z_PMCHWT_i.npy")
Z_PMCHWT_i_r = Z_PMCHWT_i.conj().T

if params_simu.CurrentModel == "SEqF_J":

    Z_CFIE = -np.r_[np.load("mat/Z_CFIE_EJ_f.npy"), np.load("mat/Z_CFIE_HJ_f.npy")]
    Z_CFIE_r = Z_CFIE.conj().T

    FF = np.r_[FF_J_th, FF_J_ph]
    FF_r = FF.conj().T

    alpha = np.c_[alpha_th, alpha_ph]
    alpha_r = alpha.conj().T

    matvec = lambda vec: np.dot(alpha, np.dot(FF, vec))
    rmatvec = lambda vec: np.dot(FF_r, np.dot(alpha_r, vec))

    Vobs_r = rmatvec(Vobs)

    matvec_all = lambda vec: rmatvec(matvec(vec))
    matvecOp = LinearOperator((N_RWG, N_RWG), matvec=matvec_all)

    I, info = gmres(matvecOp, Vobs_r, tol=1e-3)

    V_test = matvec(I)
    if (info == 0):
        print "Well converged!"
        print "True residual =", np.linalg.norm(V_test - Vobs) / np.linalg.norm(Vobs)
    else:
        print "tolerance is not archived."
        print "Number of interation =", info
        print "Apparent residual =", np.linalg.norm(matvecOp(I) - Vobs_r) / np.linalg.norm(Vobs_r)
        print "True residual =", np.linalg.norm(V_test - Vobs) / np.linalg.norm(Vobs)

    I_J, I_M = I, np.zeros(N_RWG, dtype=complex)

    np.save("./result/I_J.npy", I_J), np.save("./result/I_M.npy", I_M)

elif params_simu.CurrentModel == "DEqF":

    Z_CFIE = np.c_[np.load("mat/Z_CFIE_EJ_f.npy"),
                   np.load("mat/Z_CFIE_EM_f.npy") * params_simu.Z_ext]
    Z_CFIE_r = Z_CFIE.conj().T

    matvec_CFIE = lambda vec: np.dot(Z_CFIE, vec)
    rmatvec_CFIE = lambda vec: np.dot(Z_CFIE_r, vec)
    matvecOp_CFIE = LinearOperator((N_RWG, 2 * N_RWG), matvec=matvec_CFIE, rmatvec=rmatvec_CFIE)

    alpha = np.c_[alpha_th, alpha_ph]
    alpha_r = alpha.conj().T

    FF = np.r_[np.c_[FF_J_th, FF_M_th * params_simu.Z_ext],
               np.c_[FF_J_ph, FF_M_ph * params_simu.Z_ext]]
    FF_r = FF.conj().T

    matvec_obs = lambda vec: np.dot(alpha, np.dot(FF, vec))
    rmatvec_obs = lambda vec: np.dot(FF_r, np.dot(alpha_r, vec))

    matvecOp_obs = LinearOperator((Nobs, 2 * N_RWG),
                                  matvec=matvec_obs, rmatvec=rmatvec_obs)

    from power_iteration import power_iteration

    [obs_norm, converge_info] = power_iteration(matvecOp_obs, 100, 1e-3)
    [CFIE_norm, converge_info] = power_iteration(matvecOp_CFIE, 100, 1e-3)
    factor = obs_norm / CFIE_norm * 10.

    matvec = lambda vec: matvecOp_obs.rmatvec(matvecOp_obs(vec)) \
                         + factor * matvecOp_CFIE.rmatvec(factor * matvecOp_CFIE(vec))
    matvecOp = LinearOperator((2 * N_RWG, 2 * N_RWG), matvec=matvec)

    Vobs_r = matvecOp_obs.rmatvec(Vobs)

    I, info = gmres(matvecOp, Vobs_r, tol=1e-3, maxiter=100)

    V = np.r_[Vobs, np.zeros(N_RWG, dtype=complex)]
    V_test = np.r_[matvecOp_obs(I), factor * matvecOp_CFIE(I)]
    if (info == 0):
        print "Well converged!"
        print "True residual =", np.linalg.norm(V_test - V) / np.linalg.norm(V)
    else:
        print "tolerance is not archived."
        print "Number of interation =", info
        print "Apparent residual =", np.linalg.norm(matvecOp(I) - Vobs_r) / np.linalg.norm(Vobs_r)
        print "True residual =", np.linalg.norm(V_test - V) / np.linalg.norm(V)

    I_J, I_M = I[:N_RWG], I[N_RWG:] * params_simu.Z_ext
    np.save("./result/I_J.npy", I_J), np.save("./result/I_M.npy", I_M)

    # alpha_xyz_th = np.load("./mat/alpha_xyz_th.npy")
    # alpha_xyz_ph = np.load("./mat/alpha_xyz_ph.npy")
    # th_obs, ph_obs = np.arccos(r_obs[:, 2]), np.arctan2(r_obs[:, 1], r_obs[:, 0])
    # th_hat_obs = np.c_[np.cos(th_obs)*np.cos(ph_obs),
    #                    np.cos(th_obs)*np.sin(ph_obs), -np.sin(th_obs)]
    # ph_hat_obs = np.c_[-np.sin(ph_obs), np.cos(ph_obs), np.zeros(ph_obs.size)]
    # alpha_th_th = np.sum(ph_hat_obs[:, :, np.newaxis] * alpha_xyz_th, axis=1)
    # alpha_th_ph = np.sum(ph_hat_obs[:, :, np.newaxis] * alpha_xyz_ph, axis=1)
    # alpha_test = np.c_[alpha_th_th, alpha_th_ph]
    # def matvec_obs_test(vec):
    #     vec1 = np.dot(FF[1], np.dot(Z_PMCHWT_i, np.dot(Z_CFIE[1][0], vec)))
    #     return np.dot(alpha_test, np.dot(FF[0], vec) + vec1)
    # V_test = matvec_obs_test(I)

else:
    I_J, I_M = np.zeros(N_RWG, dtype=complex), np.zeros(N_RWG, dtype=complex)
    print "Appropriate current model is not specified. Terminated"
    sys.exit(1)

from MLFMA import computeCurrentsVisualization

for i in range(S):
    computeCurrentsVisualization(
        "./result/", params_simu, I_J,
        filename_norm="surface" + str(i) + ".norm_J_centroids_triangles.pos",
        filename_vec="surface" + str(i) + ".J_centroids_triangles.pos",
        target_surface=np.asarray([i])
    )
    computeCurrentsVisualization(
        "./result/", params_simu, I_M,
        filename_norm="surface" + str(i) + ".norm_M_centroids_triangles.pos",
        filename_vec="surface" + str(i) + ".M_centroids_triangles.pos",
        target_surface=np.asarray([i])
    )

print "complete!"
