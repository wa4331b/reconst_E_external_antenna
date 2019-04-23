import numpy as np
import scipy.special as sp


class xyz:

    def line(self, ymax=0, ystep=0, x=0, z=0):

        dy = ymax / ystep
        ystep = 2 * ystep + 1

        y = np.arange(ystep) * dy - ymax

        self.xyz = np.zeros([y.size, 3])
        self.xyz[:, 0], self.xyz[:, 1], self.xyz[:, 2] = x, y, z

    def plane(self, yzmax=0, yzstep=0, x=0):

        dyz = yzmax / yzstep
        yzstep2 = 2 * yzstep + 1
        samstep = yzstep2 ** 2

        yz = np.arange(yzstep2) * dyz - yzmax
        self.xyz = np.zeros([samstep, 3])

        m = 0
        for n in range(yzstep2):
            self.xyz[m:m + yzstep2, 1] = yz[n]
            self.xyz[m:m + yzstep2, 2] = yz
            m = m + yzstep2

        self.xyz[:, 0] = x

    def sphere_horz(self, r, step, start=0, stop=360,
                    endpoint=False, random=False):

        phi = np.linspace(start, stop, step, endpoint=endpoint)
        if (random):
            phi = np.random.rand(step) * 360.
            r = np.random.randn(step) * 0.05 * r + r
        phi = np.radians(phi)

        self.xyz = np.zeros([step, 3])
        self.xyz[:, 0] = r * np.cos(phi)
        self.xyz[:, 1] = r * np.sin(phi)
        self.xyz[:, 2] = 0

    def sphere(self, r=1., thstep=1, phstep=1):

        dth, dph = np.pi / thstep, 2 * np.pi / phstep
        [th, ph] = np.mgrid[0:thstep, 0:phstep] + 0.5

        th, ph = th * dth, ph * dph

        x = np.reshape(r * np.sin(th) * np.cos(ph), th.size)
        y = np.reshape(r * np.sin(th) * np.sin(ph), th.size)
        z = np.reshape(r * np.cos(th), th.size)

        self.xyz = np.c_[x, y, z]

    def sphere_random(self, r=0., thstep=1, phstep=1, disp=0.05):

        step = thstep * phstep

        r = np.random.randn(step) * disp * r + r
        th = np.random.rand(step) * np.pi
        ph = np.random.rand(step) * 2 * np.pi

        x = np.reshape(r * np.sin(th) * np.cos(ph), th.size)
        y = np.reshape(r * np.sin(th) * np.sin(ph), th.size)
        z = np.reshape(r * np.cos(th), th.size)

        self.xyz = np.c_[x, y, z]

    # def ellipse(self, r, num, m):
    #
    #     pi = np.pi
    #     xi = np.linspace(0, 2 * pi, num, endpoint=False)
    #
    #     # Oblate case
    #     t = 2 / pi * xi * sp.ellipe(m)
    #     if (m != 1):
    #         phi = elp.Einv(t, m)
    #         print('Ellipse.')
    #     else:
    #         phi = elp.Einv_plane(xi)
    #         print('Plane.')
    #     # Prolate case
    #     # t=(2/pi*xi-1)*sp.ellipe(m)
    #     # phi=elp.Einv(t,m)+pi/2
    #
    #     x = r * np.cos(phi)
    #     y = r * np.sin(phi)
    #     z = np.zeros(np.shape(phi))
    #
    #     self.xyz = np.c_[x, y, z]

    def cube(self, SideLength=0, Step=0):

        d = SideLength / (Step - 1)

        [x, y, z] = np.mgrid[0:Step, 0:Step, 0:Step]
        x, y, z = np.reshape(x, x.size), np.reshape(y, y.size), \
                  np.reshape(z, z.size)
        self.xyz = np.c_[x, y, z]
        self.xyz = self.xyz * d - SideLength / 2.
        self.xyz = self.xyz[np.where(np.any(np.abs(self.xyz)
                                            > (SideLength - d) / 2., axis=1))]
