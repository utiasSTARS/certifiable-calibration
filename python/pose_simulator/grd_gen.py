import numpy as np

class gnd_gen():
    def __init__(self, constants):
        if constants is None:
            self.constants = [0.08, 0.1, 3 * np.random.randn(), 5 * np.random.randn()]
        else:
            assert len(constants) == 4, 'Incorrect Number of Ground Constants'
            self.constants = constants
        return

    def calc_z(self, x, y):
        z = self.constants[2] * np.cos(self.constants[0] * (x + y)) + self.constants[3] * np.sin(self.constants[1] * (x - y)) + 3
        return z

    def calc_dz_dx(self, x, y):
        dz_dx = -self.constants[0] * self.constants[2] * np.sin(self.constants[0] * (x + y)) + self.constants[1] * self.constants[3] * np.cos(self.constants[1] * (x - y))
        return dz_dx

    def calc_dz_dy(self, x, y):
        dz_dy = -self.constants[0] * self.constants[2] * np.sin(self.constants[0] * (x + y)) - self.constants[1] * self.constants[3] * np.cos(self.constants[1] * (x - y))
        return dz_dy

    def calc_ddz_dxdx(self, x, y):
        ddz_dxdx = -self.constants[0] ** 2 * self.constants[2] * np.cos(self.constants[0] * (x + y)) - self.constants[1] ** 2 * self.constants[3] * np.sin(self.constants[1] * (x - y))
        return ddz_dxdx

    def calc_ddz_dydy(self, x, y):
        ddz_dydy = self.calc_ddz_dxdx(x, y)
        return ddz_dydy

    def calc_ddz_dxdy(self, x, y):
        ddz_dxdy = -self.constants[0] ** 2 * self.constants[2] * np.cos(self.constants[0] * (x + y)) + self.constants[1] ** 2 * self.constants[3] * np.sin(self.constants[1] * (x - y))
        return ddz_dxdy