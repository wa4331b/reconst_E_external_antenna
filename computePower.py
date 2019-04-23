import numpy as np
from inputParams.simulation_parameters import params_simu
from readRWGfiles import MeshData
from Z_CFIE import Pywrap_Z_CFIE_computation
from PenaltyMethod import *
from scipy.sparse.linalg import LinearOperator
from scipy import optimize
from ReadWriteBlitzArray import readFloatFromDisk
from MLFMA import computeCurrentsVisualization

meshPath = "./mesh/"
meshData = MeshData(meshPath)
surface_RWGnumbers, surface_NRWG = meshData.surface_rwg()
Ssurface, Dsurface = 0, 1

Z_CFIE = -np.r_[np.load("mat/Z_CFIE_EJ_f.npy")[surface_RWGnumbers[Dsurface][:, np.newaxis],
                                               surface_RWGnumbers[Ssurface]],
                np.load("mat/Z_CFIE_HJ_f.npy")[surface_RWGnumbers[Dsurface][:, np.newaxis],
                                               surface_RWGnumbers[Ssurface]]]
Z_CFIE_r = Z_CFIE.conj().T
Z_PMCHWT_i = np.load("mat/Z_PMCHWT_i.npy")
Z_PMCHWT_i_r = Z_PMCHWT_i.conj().T

FF_J_th, FF_J_ph = np.load("./mat/FF_J_th.npy"), np.load("./mat/FF_J_ph.npy")
FF_M_th, FF_M_ph = np.load("./mat/FF_M_th.npy"), np.load("./mat/FF_M_ph.npy")
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

alpha_th, alpha_ph = np.load("./mat/alpha_th.npy"), np.load("./mat/alpha_ph.npy")
alpha = np.c_[alpha_th, alpha_ph]
alpha_r = alpha.conj().T


def matvec(vec):
    ff1 = np.dot(FF[1], np.dot(Z_PMCHWT_i, np.dot(Z_CFIE, vec)))
    return np.dot(alpha, np.dot(FF[0], vec) + ff1)


def rmatvec(vec):
    arvec = np.dot(alpha_r, vec)
    ff1 = np.dot(Z_CFIE_r, np.dot(Z_PMCHWT_i_r, np.dot(FF_r[1], arvec)))
    return np.dot(FF_r[0], arvec) + ff1


# Vobs = np.load("mesh/vObs.npy")
# Vobs_r = rmatvec(Vobs)
# matvecOp = LinearOperator(
#     (Vobs.size, surface_NRWG[Ssurface]), matvec=matvec, rmatvec=rmatvec
# )
# residual = Residual(matvecOp, Vobs)
#
# print()
# print("Start optimizing...")
# res_opt = optimize.minimize(
#     ScalarComplexFunction(residual.f), np.zeros(matvecOp.shape[1] * 2),
#     method="BFGS", jac=VectorComplexFunction(residual.fp),
#     # options={"gtol": 1e-3 * np.linalg.norm(residual.aT_dot_b)}
#     tol=1e-6 * np.linalg.norm(residual.aT_dot_b)
# )
# print("End optimizing!")
# for key in res_opt.keys():
#     val = res_opt[key]
#     if type(val) != np.ndarray:
#         print(key, ":", val)
#
# print("Relative residual =", np.sqrt(res_opt["fun"] / np.vdot(Vobs, Vobs)))
# I = res_opt["x"][:res_opt["x"].size / 2] + 1j * res_opt["x"][res_opt["x"].size / 2:]
# Idielectric = np.dot(Z_PMCHWT_i, np.dot(Z_CFIE, I))
# I_J, I_M = np.zeros(meshData.N_RWG, dtype=complex), np.zeros(meshData.N_RWG, dtype=complex)
# I_J[surface_RWGnumbers[Ssurface]] = I
# I_J[surface_RWGnumbers[Dsurface]] = Idielectric[:surface_NRWG[Dsurface]]
# I_M[surface_RWGnumbers[Dsurface]] = Idielectric[surface_NRWG[Dsurface]:]
# np.save("./result/I_J.npy", I_J), np.save("./result/I_M.npy", I_M)
#
# for i in range(meshData.numberOfSurfaces):
#     computeCurrentsVisualization(
#         "./result/", params_simu, I_J,
#         filename_norm="surface" + str(i) + ".norm_J_centroids_triangles.pos",
#         filename_vec="surface" + str(i) + ".J_centroids_triangles.pos",
#         target_surface=np.asarray([i])
#     )
#     computeCurrentsVisualization(
#         "./result/", params_simu, I_M,
#         filename_norm="surface" + str(i) + ".norm_M_centroids_triangles.pos",
#         filename_vec="surface" + str(i) + ".M_centroids_triangles.pos",
#         target_surface=np.asarray([i])
#     )

# Making matrices for source surface
w = 2 * np.pi * params_simu.f
coefficient_CFIE = np.asarray([1, 0, 1, 0], dtype=complex)
bool_CFIE = np.asarray([1, 0, 1, 0], dtype="i")
testRWG = surface_RWGnumbers[Ssurface]
sourceRWG = np.arange(meshData.N_RWG, dtype="i")
z_EJ, z_EM, z_HJ, z_HM \
    = Pywrap_Z_CFIE_computation(coefficient_CFIE, bool_CFIE, True, -1., -1., sourceRWG, testRWG,
                                meshData.RWGNumber_surface,
                                meshData.RWGNumber_signedTriangles, meshData.RWGNumber_nodes,
                                meshData.nodesCoord, w, params_simu.eps_r_f, params_simu.mu_r_f,
                                np.asarray(range(meshData.numberOfSurfaces), dtype="i"))
zEJ_eachSurface, zEM_eachSurface = [], []
for i in range(meshData.numberOfSurfaces):
    zEJ_eachSurface.append(z_EJ[:, surface_RWGnumbers[i]])
    zEM_eachSurface.append(z_EM[:, surface_RWGnumbers[i]])
zDielectric = np.c_[zEJ_eachSurface[Dsurface], zEM_eachSurface[Dsurface]]


def calc_e_source(j):
    return np.dot(zEJ_eachSurface[Ssurface], j)


def calc_e_dielectric(j):
    return np.dot(zDielectric, np.dot(Z_PMCHWT_i, np.dot(Z_CFIE, j)))


def calc_e_total(j):
    return calc_e_source(j) + calc_e_dielectric(j)


def calc_e_source_real(j):
    return np.dot(zEJ_eachSurface[Ssurface].real, j)


def calc_e_dielectric_real(j):
    real_part = np.dot(zDielectric, np.dot(Z_PMCHWT_i, np.dot(Z_CFIE, j.real))).real
    imag_part = np.dot(zDielectric, np.dot(Z_PMCHWT_i, np.dot(Z_CFIE, j.imag))).real
    return real_part + 1j * imag_part


def calc_e_total_real(j):
    return calc_e_source_real(j) + calc_e_dielectric_real(j)


matvecOp_calcEtotal = LinearOperator(zEJ_eachSurface[Ssurface].shape, matvec=calc_e_total_real)

# radiationPower = 6109.1303814
# filename = './inputParams/inputFiles/' + params_simu.inputFileDir + "/Power.txt"
# radiationPower = readFloatFromDisk(filename)
quadratic0 = Quadratic(matvecOp_calcEtotal, 0.)
# quadratic = Quadratic(matvecOp_calcEtotal, -radiationPower)
power = np.sqrt(quadratic0.f(np.load("./result/I_J.npy")[surface_RWGnumbers[Ssurface]]))
print("Power =", power)
