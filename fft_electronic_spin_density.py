# include cubetools.py in folder utils/

from utils.cubetools import read_cube, write_cube
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy.fft as fft
# import bohr radius
from scipy.constants import physical_constants

fin_cube = './cube_files/Mn2GeO4_rho_sz.cube'


def main(verbose=True):
    # cube[0] is the scalar field numpy array
    cube = read_cube(fin_cube)

    nx, ny, nz = cube[0].shape

    # cube[1] contains dictionary with metadata - keys are 'org', 'xvec', 'yvec', 'zvec', 'atoms'
    # get unit cell size: 'xvec' gives 
    #    the >>spacing<< between grid points 
    #     along the first index of the array (--> a lattice vector of the unit celll) 
    #      in units of Bohr radii
    

    # ================== UNITS ==================

    # unit conversion: see http://publish.illinois.edu/yubo-paul-yang/tutorials/quantum-espresso/understand-fast-fourier-transform/
    #    and   https://en.wikipedia.org/wiki/Reciprocal_lattice (formulas for reciprocal lattice vectors in 3D)

    #--- REAL SPACE ---
    # real-space grid spacing in Angstrom
    da = np.array(list(cube[1]['xvec'])) * physical_constants['Bohr radius'][0] * 1e10  # Angstrom
    db = np.array(list(cube[1]['yvec'])) * physical_constants['Bohr radius'][0] * 1e10  # Angstrom
    dc = np.array(list(cube[1]['zvec'])) * physical_constants['Bohr radius'][0] * 1e10  # Angstrom

    # real-space lattice vectors (nx, ny, nz are the number of grid points in each direction - the dimensions of the cube numpy array)
    a =  da * nx # Angstrom
    b =  db * ny # Angstrom
    c =  dc * nz # Angstrom

    # volume of the unit cell
    V = np.dot(a, np.cross(b, c))  # Angstrom^3

    # lattice parameter matrix (Angstrom) - should match the lattice parameters in scf.in file - first row is a, second row is b, third row is c
    A = np.vstack((a, b, c))

    #--- RECIPROCAL SPACE ---

    # get reciprocal lattice vectors
    ka = 2 * np.pi * np.cross(b,c) / V  # 1/Angstrom
    kb = 2 * np.pi * np.cross(c,a) / V  # 1/Angstrom
    kc = 2 * np.pi * np.cross(a,b) / V  # 1/Angstrom

    # reciprocal spacings
    dka = ka / nx  # 1/Angstrom
    dkb = kb / ny  # 1/Angstrom
    dkc = kc / nz  # 1/Angstrom

    # reciprocal lattice parameter matrix (1/Angstrom) - first row is ka, second row is kb, third row is kc
    B = np.vstack((ka, kb, kc))

    if verbose:
        print('A (Angstrom)\n', A)
        print('\nB (1/Angstrom)\n', B)

    # plot_cube_file(cube)

    # get FFT
    # no prefactor applied
    cube_fft = fft.fftn(cube[0], norm='backward')
    cube_fft_abs = np.abs(cube_fft)

    print('cube_fft numpy array shape', cube_fft.shape)
    print('type(cube_fft[0,0,])', type(cube_fft_abs[0,0,0]))

    # sum over c axis
    cube_fft_xy = np.sum(cube_fft_abs, axis=2)
    plot_2D_fft(cube_fft_xy)


def plot_cube_file(cube, fout_name='rho_sz.png'):
    # plot numpy array as a scalar field
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Create a 3D grid
    x = np.arange(cube[0].shape[0])
    y = np.arange(cube[0].shape[1])
    z = np.arange(cube[0].shape[2])
    X, Y, Z = np.meshgrid(x, y, z)

    # Plot the scalar field
    plot = ax.scatter(X, Y, Z, c=cube[0].flatten(), cmap='viridis')
    ax.set_aspect('equal', adjustable='box')
    # plot colorbar the colorbar
    plt.colorbar(plot)
    plt.tight_layout()
    plt.savefig(fout_name)
    plt.close()


def plot_2D_fft(cube_fft_xy, fout_name='rho_sz_fft.png'):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Create a 2D grid
    x = np.arange(cube_fft_xy.shape[0])
    y = np.arange(cube_fft_xy.shape[1])
    X, Y = np.meshgrid(x, y)

    # Plot the scalar field
    plot = ax.scatter(X, Y, c=cube_fft_xy.flatten(), cmap='Greys')
    # plot colorbar the colorbar
    plt.colorbar(plot)
    plt.tight_layout()
    plt.savefig(fout_name)
    plt.close()


if __name__ == '__main__':
    main()