import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys, os


def ReadData(filename):
    data = np.loadtxt(filename, skiprows=1)
    xyz = data[:, :3]
    fx = data[:, 3] + 1j * data[:, 4]
    fy = data[:, 5] + 1j * data[:, 6]
    fz = data[:, 7] + 1j * data[:, 8]

    nyz = np.sqrt(fz.size).astype("i")
    fx = np.reshape(fx, [nyz, nyz])
    fy = np.reshape(fy, [nyz, nyz])
    fz = np.reshape(fz, [nyz, nyz])

    return xyz, [fx, fy, fz]


def PlanarPlot(r_obs, f, clim=None, title=None, patch=None, if_clip=False):
    fig, axs = plt.subplots(1, 1, figsize=(5, 4))
    plt.axis("off")
    l = np.max(r_obs)
    im = axs.imshow(f.T[::-1, :], clim=clim, interpolation="nearest")

    if patch != None:
        axs.add_patch(patch)
        if if_clip: im.set_clip_path(patch)
    axs.set_title(title)
    # axs.set_xlim(-l, l), axs.set_ylim(l, -l)
    axs.set_xlabel("y [m]"), axs.set_ylabel("x [m]")
    fig.colorbar(im, ax=axs)
    plt.tight_layout()


if __name__ == "__main__":
    filepath = "./result"

    r_obs, E = ReadData(filepath + "/E_obs_f.txt")
#    PlanarPlot(r_obs, np.linalg.norm(np.asarray(E), axis=0), clim=[0, 1e+4], title="E [V/m]")
    PlanarPlot(r_obs, np.linalg.norm(np.asarray(E), axis=0), clim=[0, 10], title="E [V/m]")
    plt.savefig(filepath + "/E_obs_f.png")

    r_obs, E = ReadData(filepath + "/E_obs_d.txt")
    sphere = patches.Circle((E[0].shape[0] * 0.5 - 0.5, E[0].shape[0] * 0.5 - 0.5),
                            E[0].shape[0] * 0.5, linewidth=1, edgecolor='w',
                            facecolor='none')
    factor = (E[0].shape[0] - 1) / 0.1
    height = 0.068 * factor
    width = 0.088 * factor
    center = ((E[0].shape[0] - 1) * 0.5, (E[0].shape[1] - 1) * 0.5)
    lower_corner = (center[0] - width * 0.5, center[1] - height * 0.5)
    rect = patches.Rectangle(lower_corner, width, height,
                             linewidth=1, edgecolor="w", facecolor="none")

    Enorm = np.linalg.norm(np.asarray(E), axis=0)
    rnorm = np.linalg.norm(r_obs[:, 1:], axis=1).reshape(Enorm.shape)
    y_obs = r_obs[:, 1].reshape(Enorm.shape)
    z_obs = r_obs[:, 2].reshape(Enorm.shape)
    ind = np.where(
        np.logical_and(np.logical_and(-0.03 < z_obs, z_obs < 0.03),
                       np.logical_and(-0.04 < y_obs, y_obs < 0.04))
    )
    # Enorm = Enorm/np.max(Enorm[np.where(rnorm<0.03)])
#    PlanarPlot(r_obs, Enorm, clim=[0, 1e+4], title="E (normalized)")
    PlanarPlot(r_obs, Enorm, clim=[0, 0.1], title="E (normalized)")
    # PlanarPlot(r_obs, np.linalg.norm(np.asarray(E), axis=0), clim=[0, 0.08],
    #            title="E [V/m]", if_clip=True)
    plt.savefig(filepath + "/E_obs_d.png")

    plt.show()
