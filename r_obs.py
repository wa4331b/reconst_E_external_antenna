import numpy as np
from xyz import xyz
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def threed_plot(r_obs):
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.axis("off")
    ax.scatter3D(r_obs[:, 0], r_obs[:, 1], r_obs[:, 2], marker=".")
    maxval = np.max(abs(r_obs))
    plt.ylim(-maxval, maxval)
    plt.xlim(-maxval, maxval)
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.savefig("SphScan.eps", bbox_inches="tight")
    plt.show()


input_dir_name = "./inputParams"

obs_f = xyz()
obs_d = xyz()

# for dielectric sphere measurement
offset_f = np.asarray([0., 0., 0.05])
offset_d = offset_f
obs_d.plane(yzmax=0.2, yzstep=100, x=0)
obs_f.plane(yzmax=0.2, yzstep=100, x=0)
obs_d.xyz = np.roll(obs_d.xyz, 2, axis=1)
obs_f.xyz = np.roll(obs_f.xyz, 2, axis=1)
obs_f.xyz += offset_f
obs_d.xyz += offset_d

# obs_d.xyz[:, 0] += 0.02 + 0.05

# for dielectric plate measurement
# obs_f.plane(yzmax=0.1, yzstep=40, x=0)
# obs_d.plane(yzmax=0.05, yzstep=50, x=(0.06 + 0.0675) * 0.5)

np.savetxt(input_dir_name + '/r_obs_f.txt', obs_f.xyz)
np.savetxt(input_dir_name + '/r_obs_d.txt', obs_d.xyz)
# np.savetxt('r_obs.asc', obs.xyz*1000.)

print('Complete!')
