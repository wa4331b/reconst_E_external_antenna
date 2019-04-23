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


    Vobs_r = rmatvec(Vobs)

    matvec_all = lambda vec: rmatvec(matvec(vec))
    matvecOp = LinearOperator((surface_NRWG[Ssurface], surface_NRWG[Ssurface]), matvec=matvec_all)

    if params_simu.useInit:
        I0 = np.load(
            "./inputParams/initSolutions/" + params_simu.J_init_FILENAME
        )[surface_RWGnumbers[Ssurface]]
    else:
        I0 = np.zeros(surface_NRWG[Ssurface])
    solutionVector, info \
        = gmres(matvecOp, Vobs_r, tol=params_simu.TOL, maxiter=params_simu.MAXITER, x0=I0)
    dielectricVector = np.dot(Z_PMCHWT_i, np.dot(Z_CFIE, solutionVector))

    V_test = matvec(solutionVector)
    if info == 0:
        print("Well converged!")
        print("True residual =", np.linalg.norm(V_test - Vobs) / np.linalg.norm(Vobs))
    else:
        print("tolerance is not archived.")
        print("Number of iteration =", info)
        print("Apparent residual =",
              np.linalg.norm(matvecOp(solutionVector) - Vobs_r) / np.linalg.norm(Vobs_r))
        print("True residual =", np.linalg.norm(V_test - Vobs) / np.linalg.norm(Vobs))

    I_J, I_M = np.zeros(N_RWG, dtype=complex), np.zeros(N_RWG, dtype=complex)

    I_J[surface_RWGnumbers[Ssurface]] = solutionVector
    I_J[surface_RWGnumbers[Dsurface]] = dielectricVector[:surface_NRWG[Dsurface]]
    I_M[surface_RWGnumbers[Dsurface]] = dielectricVector[surface_NRWG[Dsurface]:]

    np.save("./result/I_J.npy", I_J), np.save("./result/I_M.npy", I_M)

    if params_simu.save_as_Init:
        np.save("./inputParams/initSolutions/" + params_simu.J_init_FILENAME_tosave, I_J)
        np.save("./inputParams/initSolutions/" + params_simu.M_init_FILENAME_tosave, I_M)
else:
    print("Appropriate current model is not specified. Terminated")
    sys.exit(1)

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

print("complete!")
