# we import the base parameters
# if you want to modify them, you can do it here
# or in MLFMA_parameters.py
from . import Params_Simu as Params_Simu
from . import EM_constants
import numpy as np

params_simu = Params_Simu.Params_Simu()

# the geometry or target we are going to use in the simulations.
# The list of targets can be viewed by looking up the geo directory
params_simu.pathToTarget = './geo'
# params_simu.targetName = "samPhantom&box_coarse_2"
params_simu.targetName = "cylinder&sphere"
params_simu.meshFileName = params_simu.targetName + '.msh'
params_simu.ifSourceOnly = False
if not params_simu.ifSourceOnly:
    params_simu.dielectricSurface = 1

obsName = "dipole&Dsphere_cartesian"
# obsName = "dipole&Dsphere_HdipoleProbe"
params_simu.inputFileDir = obsName

params_simu.useInit = False
if params_simu.useInit:
    params_simu.J_init_FILENAME = "I_J_cylinderWide&sphere_dipole&DsphereEps42.npy"
    params_simu.M_init_FILENAME = "I_M_cylinderWide&sphere_dipole&DsphereEps42.npy"

params_simu.save_as_Init = False
if params_simu.save_as_Init:
    params_simu.J_init_FILENAME_tosave \
        = "I_J_" + params_simu.targetName + "_" + obsName + ".npy"
    params_simu.M_init_FILENAME_tosave \
        = "I_M_" + params_simu.targetName + "_" + obsName + ".npy"

params_simu.ReadFFfromFile = False
if params_simu.ReadFFfromFile:
    params_simu.FFfilename = "dipole.ffe"

params_simu.TOL = 1e-6
params_simu.MAXITER = 100

# frequency
params_simu.f = 2.5e9

# relative epsilon and mu of the host medium
params_simu.eps_r_f = 1. + 0.j
params_simu.mu_r_f = 1. + 0.j

eps_r_real = 40.
# sigma = 1.
sigma = 0.
w = 2 * np.pi * params_simu.f
eps_r_imag = -sigma / (w * EM_constants.eps_0)
params_simu.eps_r_d = eps_r_real + 1j * eps_r_imag
params_simu.mu_r_d = 1. + 0.j

Z_f = 120 * np.pi
params_simu.Z_ext = Z_f * np.sqrt(params_simu.mu_r_f / params_simu.eps_r_f)
params_simu.Z_int = Z_f * np.sqrt(params_simu.mu_r_d / params_simu.eps_r_d)
params_simu.CFIE_ext = np.asarray([1., 1., params_simu.Z_ext, params_simu.Z_ext], dtype=complex)
params_simu.CFIE_int = np.asarray([1., 1., params_simu.Z_int, params_simu.Z_int], dtype=complex)

CurrentModels = ["DEqF", "SEqF_J"]
CurrentModel_num = 1  # Change here
params_simu.CurrentModel = CurrentModels[CurrentModel_num]

# Bistatic computation settings
# 1. EXCITATION
###############
# type, strenght and phase and direction, origin of the source.
# It will not be used in case of monostatic_RCS computation
# we have 2 types of excitation: 
# - plurality of dipoles (1 dipole is ok)
# - plane wave
# We can have these excitations separately or at the same time, i.e.,
# plane wave and dipoles.
params_simu.BISTATIC_EXCITATION_DIPOLES = 0
params_simu.BISTATIC_EXCITATION_PLANE_WAVE = 1
# now the details of each excitation
# the name (with path) of the user-supplied excitation file. Set to "" if empty
params_simu.BISTATIC_EXCITATION_J_DIPOLES_FILENAME = "J_dipoles.txt"
params_simu.BISTATIC_EXCITATION_M_DIPOLES_FILENAME = ""
# the structure of the excitation file MUST BE AS FOLLOWS:
# 1 line per dipole, as many lines as there are dipoles
# each line has 9 columns that must be arranged as follows:
#
# real(J_x) imag(J_x) real(J_y) imag(J_y) real(J_z) imag(J_z) r_x r_y r_z
#
# where J = [J_x J_y J_z] is the dipole and r = [r_x r_y r_z] its origin.

if params_simu.BISTATIC_EXCITATION_PLANE_WAVE == 1:
    # origin, strength, phase and polarization of the plane wave
    params_simu.theta_inc = np.pi / 2.0
    params_simu.phi_inc = 0.
    params_simu.E_inc_theta = 1.0 + 0.j  # the theta component
    params_simu.E_inc_phi = 0.0 + 0.j  # the phi component

# 2. OBSERVATION
################
# sampling points: sampling of the resulting field at user-specified points in space.
# It will be used only for BISTATIC
params_simu.BISTATIC_R_OBS = 1
# the name (with path) of the user-supplied r_obs file. Set to "" if empty
params_simu.BISTATIC_R_OBS_FILENAME = "r_obs.txt"
params_simu.DIPOLE_POLARIZATION_FILENAME = "DIPOLE_POLARIZATION_ANGLES.txt"
# the structure of the r_obs file MUST BE AS FOLLOWS:
# 1 line per observation point, as many lines as there are points
# each line has 3 columns that must be arranged as follows:
#
# r_obs_x r_obs_y r_obs_z

# we can also have sampling angles from a user-input file
params_simu.BISTATIC_ANGLES_OBS = 0
# the name (with path) of the user-supplied r_obs file. Set to "" if empty
params_simu.BISTATIC_ANGLES_OBS_FILENAME = "bistatic_angles_obs.txt"
# the structure of the bistatic _angles_obs file MUST BE AS FOLLOWS:
# 1 line per observation angle, as many lines as there are angles
# each line has 2 columns which are the angles in degrees (easier for human reading).
# We must have: 0 <= theta <= 180 degrees (0 is the z axis, 180 is -z).
# likewise, 0 <= phi <= 360 degrees. For theta = 90 degrees, phi = 0 is the x axis, 
# 180 is -x, 360 is x.
#
# theta_obs phi_obs

# the angles for the monostatic RCS or the bistatic far-field data.
# Normally the code provides "best angles of observation", best
# from a "field spatial information for minimal sampling size" point of view.
# if you want the program to choose the best points, set AUTOMATIC to True.
# if you want to provide your own sampling points (because you need less/more angles,
# for example), set AUTOMATIC to False. If you chose 1 point only, and AUTOMATIC is False,
# then the point will correspond to START_THETA for thetas and START_PHI for phis.
# thetas
params_simu.START_THETA = np.pi / 2.0
params_simu.STOP_THETA = np.pi
params_simu.AUTOMATIC_THETAS = False
params_simu.USER_DEFINED_NB_THETA = 1
# phis
params_simu.START_PHI = 0.
params_simu.STOP_PHI = 2 * np.pi
params_simu.AUTOMATIC_PHIS = False
params_simu.USER_DEFINED_NB_PHI = 46
# monostatic angles can be defined from a file. 
# They will take precedence over the previous definition.
params_simu.ANGLES_FROM_FILE = 0
if params_simu.ANGLES_FROM_FILE == 1:
    # the name (with path) of the user-supplied r_obs file. Set to "" if empty
    params_simu.ANGLES_FILENAME = "monostatic_RCS_angles.txt"
    # the structure of the angles file MUST BE AS FOLLOWS:
    # 1 line per angle, as many lines as there are angles
    # each line has 2 columns  which are the angles in degrees (easier for human reading):
    #
    # theta phi

# figure that shows the far field or the RCS of the target
params_simu.SHOW_FIGURE = 1

# CURRENTS_VISUALIZATION tells if we want to create a view of the currents in
# GMSH or not. if 1, you can use GMSH to view the currents that have been created.
# open with GMSH the *.pos files that are in the result directory
# only works if BISTATIC == 1
params_simu.CURRENTS_VISUALIZATION = 1
# we can also save the currents at the centroids of the triangles.
# It will be stored in the result directory
params_simu.SAVE_CURRENTS_CENTROIDS = 1

# MLFMA parameters
params_simu.a_factor = 0.25
params_simu.NB_DIGITS = 6
params_simu.int_method_theta = "GAUSSL"
params_simu.int_method_phi = "PONCELET"
params_simu.INCLUDE_BOUNDARIES = 0
params_simu.alphaTranslation_smoothing_factor = 1.4
params_simu.alphaTranslation_thresholdRelValueMax = 1.0e-3
params_simu.alphaTranslation_RelativeCountAboveThreshold = 0.6
