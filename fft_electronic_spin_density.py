# include cubetools.py in folder utils/

from utils.cubetools import read_cube, write_cube
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

fin_cube = './cube_files/Mn2GeO4_rho_sz.cube'


def main():
    # cube[0] is the scalar field numpy array, cube[1] contains metadata - info about atoms etc.
    cube = read_cube(fin_cube)

    print('cube[0] numpy array shape', cube[0].shape)

    # plot_cube_file(cube)

    # get FFT
    cube_fft = np.fft.fftn(cube[0])
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