import sys, os, argparse
import numpy as np

try:
    import cPickle
except ImportError:
    import pickle as cPickle
from scipy import array, sqrt
from ReadWriteBlitzArray import readIntFromDisk, writeScalarToDisk, readFloatFromDisk
from inputParams.simulation_parameters import params_simu
from inputParams.EM_constants import eps_0, mu_0, pi


def setup_mesh(meshPath):
    # size of cube at finest level

    average_RWG_length = readFloatFromDisk(os.path.join(meshPath, 'average_RWG_length.txt'))
    a = average_RWG_length * params_simu.a_factor

    w = 2. * pi * params_simu.f
    k_f = (w * sqrt(eps_0 * params_simu.eps_r_f * mu_0 * params_simu.mu_r_f)).astype(complex)

    writeScalarToDisk(a, meshPath + 'octtree_data/leaf_side_length.txt')
    writeScalarToDisk(2.0 * pi * params_simu.f, os.path.join(meshPath, 'octtree_data/w.txt'))
    writeScalarToDisk(
        params_simu.eps_r_f, os.path.join(meshPath, 'octtree_data/free_space/eps_r.txt')
    )
    writeScalarToDisk(
        params_simu.mu_r_f, os.path.join(meshPath, 'octtree_data/free_space/mu_r.txt')
    )
    writeScalarToDisk(k_f, os.path.join(meshPath, 'octtree_data/free_space/k.txt'))

    np.save(
        os.path.join(meshPath, 'octtree_data/free_space/CFIEcoeffs'),
        np.asarray(params_simu.CFIE_ext).astype(complex)
    )


if __name__ == '__main__':
    setup_mesh("./mesh/")
