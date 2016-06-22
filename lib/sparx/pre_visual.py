import matplotlib.pyplot as plt

def plot(mesh,phys):
        r = mesh.R_c
        plt.plot(r,phys.n_H2)
        plt.ylabel('H2 density (m^-3)')
        plt.show()