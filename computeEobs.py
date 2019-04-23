import numpy as np
from inputParams.simulation_parameters import params_simu
from OvserveField import ObserveField, saveObsFile
from readRWGfiles import MeshData

inputDirName = "./inputParams"
mesh_path = "./mesh/"
result_data_path = "./result/"

I_J, I_M = np.load("./result/I_J.npy"), np.load("./result/I_M.npy")

meshData = MeshData(mesh_path)
observeField = ObserveField(meshData, I_J, I_M, params_simu.f)

print("Computing near-field of free space...")
r_obs = np.loadtxt(inputDirName + "/r_obs_f.txt")
target_surface = np.arange(meshData.numberOfSurfaces, dtype="i")
E_obs, H_obs = observeField.computeFields(r_obs, params_simu.eps_r_f,
                                          params_simu.mu_r_f, target_surface)
saveObsFile(r_obs, E_obs, result_data_path + "E_obs_f.txt")
saveObsFile(r_obs, H_obs, result_data_path + "H_obs_f.txt")

if not params_simu.ifSourceOnly:
    print("Computing near-field of dielectric space...")
    r_obs = np.loadtxt(inputDirName + "/r_obs_d.txt")
    target_surface = np.array([params_simu.dielectricSurface], dtype="i")
    E_obs, H_obs, = observeField.computeFields(r_obs, params_simu.eps_r_d,
                                               params_simu.mu_r_d, target_surface)
    saveObsFile(r_obs, E_obs, result_data_path + "E_obs_d.txt")
    saveObsFile(r_obs, H_obs, result_data_path + "H_obs_d.txt")

with open("PlanarPlot.py") as f:
    code = f.read()
    exec(code)
print("Computation finished!")
