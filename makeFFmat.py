from E_FF import Pywrap_compute_FFmat
from make_trans import Pywrap_makeTrans, Pywrap_makeFFforTrans
from Z_CFIE import Pywrap_Z_CFIE_computation
import numpy as np
from ReadWriteBlitzArray import *
from inputParams.simulation_parameters import params_simu
from inputParams.EM_constants import *
from readRWGfiles import MeshData
import pickle

def read_ff_from_ffe(filename):
    with open(filename) as file:
        line_indicator_nth = "#No. of Theta Samples"
        line_indicator_nph = "#No. of Phi Samples"
        nth, nph = 0, 0
        for line in file:
            if line[:len(line_indicator_nth)] == line_indicator_nth:
                nth = int(line.split(":")[1])
            if line[:len(line_indicator_nph)] == line_indicator_nph:
                nph = int(line.split(":")[1])
            if (nth != 0) and (nph != 0):
                break

    ffe_data = np.loadtxt(filename, comments=("#", "**"))
    theta = ffe_data[:nth, 0]
    phi = ffe_data[::nth, 1]

    e_theta = ffe_data[:, 2] + 1j * ffe_data[:, 3]
    e_phi = ffe_data[:, 4] + 1j * ffe_data[:, 5]
    e_theta = e_theta.reshape((nph, nth)).T
    e_phi = e_phi.reshape((nph, nth)).T

    e_max = max(np.max(np.abs(e_theta)), np.max(np.abs(e_phi)))
    e_theta /= e_max
    e_phi /= e_max

    return e_theta, e_phi, theta, phi, nth, nph


if __name__ == "__main__":

    mesh_path = "./mesh/"
    result_path = "./result/"
    if not "mat" in os.listdir():
        os.mkdir("mat")
    files = os.listdir("./mat/")
    for matFileName in files:
        os.remove("./mat/" + matFileName)

    meshData = MeshData(mesh_path)
    N_RWG = meshData.RWGNumber_nodes.shape[0]

    surface_RWGNumbers = []
    surface_NRWG = []
    for i in range(meshData.numberOfSurfaces):
        surface_RWGNumbers.append(np.where(meshData.RWGNumber_surface == i)[0].astype("i"))
        surface_NRWG.append(surface_RWGNumbers[i].size)

    with open("mat/surface_RWG.pickle", "wb") as f:
        pickle.dump((surface_RWGNumbers, surface_NRWG), f)

    r_obs = np.load(mesh_path + "r_obs.npy")
    w = 2 * np.pi * params_simu.f

    print("Computing matrix elements of observation...")

    rMax = np.max(np.linalg.norm(meshData.nodesCoord, axis=1))
    vc = 1. / np.sqrt(eps_0 * params_simu.eps_r_f * mu_0 * params_simu.mu_r_f)
    k0 = w / vc
    d0 = 6
    kr = np.abs(k0 * rMax)
    L = int(kr + 1.8 * d0 ** (2. / 3.) * kr ** (1. / 3.)) + 1
    nTheta, nPhi = 2 * L, 4 * L
    nDirection = nTheta * nPhi
    xGauss, wGauss = np.polynomial.legendre.leggauss(nTheta)
    thetas = np.arccos(xGauss)[::-1]
    [phis, dph] = np.linspace(0, 2 * np.pi, nPhi, endpoint=False, retstep=True)

    wThetas = np.tile(wGauss, [nPhi, 1]).T.reshape(nDirection)
    thetas = np.tile(thetas, [nPhi, 1]).T.reshape(nDirection)
    phis = np.tile(phis, [nTheta])

    FF_J_th, FF_J_ph, FF_M_th, FF_M_ph \
        = Pywrap_compute_FFmat(thetas, phis, meshData.RWGNumber_trianglesCoord,
                               meshData.triangle_NRWG, meshData.triangle_RWGNumber,
                               meshData.triangle_signInRWG,
                               meshData.triangle_surface, w,
                               params_simu.eps_r_f, params_simu.mu_r_f)
    np.save("./mat/FF_J_th.npy", FF_J_th), np.save("./mat/FF_J_ph.npy", FF_J_ph)
    np.save("./mat/FF_M_th.npy", FF_M_th), np.save("./mat/FF_M_ph.npy", FF_M_ph)
    del FF_J_th, FF_J_ph, FF_M_th, FF_M_ph

    alpha = Pywrap_makeTrans(r_obs, k0, nTheta, thetas, phis)
    weight = wThetas * dph
    alpha *= weight

    if params_simu.ReadFFfromFile:
        ff0_th, ff0_ph, theta0, phi0, nTheta0, nPhi0 \
            = read_ff_from_ffe("./inputParams/FarField/" + params_simu.FFfilename)
        ff0_th = ff0_th.reshape(nTheta0 * nPhi0)
        ff0_ph = ff0_ph.reshape(nTheta0 * nPhi0)
        theta0 = np.deg2rad(theta0)
        phi0 = np.deg2rad(phi0)
    else:
        nTheta0, nPhi0 = 2, 1
        theta0, dth0 = np.linspace(0, np.pi, nTheta0, endpoint=False, retstep=True)
        theta0 += dth0 * 0.5
        phi0 = np.linspace(0, 2 * np.pi, nPhi0, endpoint=False)
        nDirection0 = nTheta0 * nPhi0
        ff0_th = np.tile(-np.sin(theta0), [nPhi0, 1]).T.reshape(nDirection0).astype(complex)
        ff0_ph = np.zeros(nDirection0, dtype=complex)

    polAngles = np.load("./mesh/polAngles.npy")
    ff_th, ff_ph = Pywrap_makeFFforTrans(polAngles, theta0, phi0, ff0_th, ff0_ph, thetas, phis)

    rObsNumbers = np.load("./mesh/rObsNumbers.npy")
    polNumbers = np.load("./mesh/polNumbers.npy")
    alpha_th = alpha[rObsNumbers, :] * ff_th[polNumbers, :]
    alpha_ph = alpha[rObsNumbers, :] * ff_ph[polNumbers, :]
    np.save("./mat/alpha_th.npy", alpha_th), np.save("./mat/alpha_ph.npy", alpha_ph)
    del alpha_th, alpha_ph

    print("Computing matrix elements of external space...")
    target_surface = np.arange(meshData.numberOfSurfaces, dtype="i")
    CFIE_bool = np.asarray([1, 0, 1, 0], dtype="i")
    RWG_tmp = np.arange(N_RWG, dtype="i")
    Z_CFIE_EJ_f, Z_CFIE_EM_f, Z_CFIE_HJ_f, Z_CFIE_HM_f \
        = Pywrap_Z_CFIE_computation(params_simu.CFIE_ext, CFIE_bool, True, 1., 1.,
                                    RWG_tmp, RWG_tmp,
                                    meshData.RWGNumber_surface,
                                    meshData.RWGNumber_signedTriangles,
                                    meshData.RWGNumber_nodes, meshData.nodesCoord,
                                    w, params_simu.eps_r_f, params_simu.mu_r_f,
                                    target_surface)
    np.save("mat/Z_CFIE_EJ_f.npy", Z_CFIE_EJ_f), np.save("mat/Z_CFIE_EM_f.npy", Z_CFIE_EM_f)
    np.save("mat/Z_CFIE_HJ_f.npy", Z_CFIE_HJ_f), np.save("mat/Z_CFIE_HM_f.npy", Z_CFIE_HM_f)

    if not params_simu.ifSourceOnly:
        print("Computing matrix elements of internal space...")
        Dsurface = params_simu.dielectricSurface
        Z_CFIE_EJ_f = Z_CFIE_EJ_f[surface_RWGNumbers[Dsurface][:, np.newaxis],
                                  surface_RWGNumbers[Dsurface]]
        Z_CFIE_EM_f = Z_CFIE_EM_f[surface_RWGNumbers[Dsurface][:, np.newaxis],
                                  surface_RWGNumbers[Dsurface]]
        Z_CFIE_HJ_f = Z_CFIE_HJ_f[surface_RWGNumbers[Dsurface][:, np.newaxis],
                                  surface_RWGNumbers[Dsurface]]
        Z_CFIE_HM_f = Z_CFIE_HM_f[surface_RWGNumbers[Dsurface][:, np.newaxis],
                                  surface_RWGNumbers[Dsurface]]

        target_surface = np.arange(meshData.numberOfSurfaces, dtype="i")
        CFIE_bool = np.asarray([1, 0, 1, 0], dtype="i")
        Z_CFIE_EJ_d, Z_CFIE_EM_d, Z_CFIE_HJ_d, Z_CFIE_HM_d \
            = Pywrap_Z_CFIE_computation(params_simu.CFIE_int, CFIE_bool, True, -1., -1.,
                                        surface_RWGNumbers[params_simu.dielectricSurface],
                                        surface_RWGNumbers[params_simu.dielectricSurface],
                                        meshData.RWGNumber_surface,
                                        meshData.RWGNumber_signedTriangles,
                                        meshData.RWGNumber_nodes, meshData.nodesCoord,
                                        w, params_simu.eps_r_d, params_simu.mu_r_d, target_surface)

        print("Solving PMCHWT")
        Z_PMCHWT = np.r_[np.c_[Z_CFIE_EJ_f + Z_CFIE_EJ_d, Z_CFIE_EM_f + Z_CFIE_EM_d],
                         np.c_[Z_CFIE_HJ_f + Z_CFIE_HJ_d, Z_CFIE_HM_f + Z_CFIE_HM_d]]
        del Z_CFIE_EJ_f, Z_CFIE_EM_f, Z_CFIE_HJ_f, Z_CFIE_HM_f
        del Z_CFIE_EJ_d, Z_CFIE_EM_d, Z_CFIE_HJ_d, Z_CFIE_HM_d
        Z_PMCHWT_i = np.linalg.inv(Z_PMCHWT)
        np.save("mat/Z_PMCHWT.npy", Z_PMCHWT), np.save("mat/Z_PMCHWT_i", Z_PMCHWT_i)

    print("Complete!")
