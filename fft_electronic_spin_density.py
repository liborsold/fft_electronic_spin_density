# include cubetools.py in folder utils/

from utils.cubetools import read_cube, write_cube
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy.fft as fft
# import bohr radius
from scipy.constants import physical_constants

fin_cube = './cube_files/Mn2GeO4_rho_sz.cube'


def main():
    # cube[0] is the scalar field numpy array
    cube = read_cube(fin_cube)

    nx, ny, nz = cube[0].shape

    # cube[1] contains dictionary with metadata - keys are 'org', 'xvec', 'yvec', 'zvec', 'atoms'
    # get unit cell size: 'xvec' gives 
    #    the >>spacing<< between grid points 
    #     along the first index of the array (--> a lattice vector of the unit celll) 
    #      in units of Bohr radii
    

    # ================== UNITS ==================

    #--- REAL SPACE ---
    # get real-space grid spacing in Angstrom
    da = np.array(list(cube[1]['xvec'])) * physical_constants['Bohr radius'][0] * 1e10  # Angstrom
    db = np.array(list(cube[1]['yvec'])) * physical_constants['Bohr radius'][0] * 1e10  # Angstrom
    dc = np.array(list(cube[1]['zvec'])) * physical_constants['Bohr radius'][0] * 1e10  # Angstrom

    # get lattice vectors
    a =  da * nx # Angstrom
    b =  db * ny # Angstrom
    c =  dc * nz # Angstrom

    # stack into single matrix - should match the lattice parameters in scf.in file
    A_latt = np.vstack((a, b, c))
    print('A_latt', A_latt)

    #--- RECIPROCAL SPACE ---
    # get reciprocal spacings
    dka = 2 * np.pi / a # 1/Angstrom
    dkb = 2 * np.pi / b # 1/Angstrom
    dkc = 2 * np.pi / c # 1/Angstrom

    # get reciprocal lattice vectors
    ka = dka * nx # 1/Angstrom
    kb = dkb * ny # 1/Angstrom
    kc = dkc * nz # 1/Angstrom

    # stack into single matrix
    B_latt = np.vstack((ka, kb, kc))
    print('B_latt', B_latt)


    # unit conversion: see http://publish.illinois.edu/yubo-paul-yang/tutorials/quantum-espresso/understand-fast-fourier-transform/:
    #    In condensed matter, consider a 1D box of length L. For an FFT grid of size N, the real space spacing is dx = L/N. The reciprocal space box length is 2πN/L, so the spacing is dk = 2π/L.
    # -> dka = 2*pi / a_latt = 2*pi / (da * nx)

    exit()

    print('cube[0] numpy array shape', cube[0].shape)

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