# include cubetools.py in folder utils/

from utils.cubetools import read_cube, write_cube
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy.fft as fft
# import bohr radius
from scipy.constants import physical_constants

fin_cube = './cube_files/Mn2GeO4_rho_1024.cube'


class Density:

    def __init__(self, kz_target=5, verbose=True, plot_real_space_spin_density=False):
        """_summary_

        Args:
            kz (int, optional): The selected cut at k_z in 1/Angstrom. Defaults to 5.
            verbose (bool, optional): _description_. Defaults to True.
            plot_real_space_spin_density (bool, optional): _description_. Defaults to False.
        """

        # ================== READ CUBE FILE ==================

        # cube[0] is the scalar field numpy array
        cube_data = read_cube(fin_cube)

        # cube[1] contains dictionary with metadata - keys are 'org', 'xvec', 'yvec', 'zvec', 'atoms'
        # get unit cell size: 'xvec' gives 
        #    the >>spacing<< between grid points 
        #     along the first index of the array (--> a lattice vector of the unit celll) 
        #      in units of Bohr radii
        # ... being used below
        
        # numpy array dimensions
        self.array = cube_data[0]
        self.nx, self.ny, self.nz = self.array.shape


        # ================== UNITS ==================

        # unit conversion: see http://publish.illinois.edu/yubo-paul-yang/tutorials/quantum-espresso/understand-fast-fourier-transform/
        #    and   https://en.wikipedia.org/wiki/Reciprocal_lattice (formulas for reciprocal lattice vectors in 3D)

        #--- REAL SPACE ---
        # real-space grid spacing in Angstrom
        self.da = np.array(list(cube_data[1]['xvec'])) * physical_constants['Bohr radius'][0] * 1e10  # Angstrom
        self.db = np.array(list(cube_data[1]['yvec'])) * physical_constants['Bohr radius'][0] * 1e10  # Angstrom
        self.dc = np.array(list(cube_data[1]['zvec'])) * physical_constants['Bohr radius'][0] * 1e10  # Angstrom

        # real-space lattice vectors (nx, ny, nz are the number of grid points in each direction - the dimensions of the cube numpy array)
        self.a =  self.da * self.nx # Angstrom
        self.b =  self.db * self.ny # Angstrom
        self.c =  self.dc * self.nz # Angstrom

        # volume of the unit cell
        self.V = np.dot(self.a, np.cross(self.b, self.c))  # Angstrom^3

        # lattice parameter matrix (Angstrom) - should match the lattice parameters in scf.in file - first row is a, second row is b, third row is c
        self.A = np.vstack((self.a, self.b, self.c))

        #--- RECIPROCAL SPACE ---

        # get reciprocal lattice spacings
        self.dka = 2 * np.pi * np.cross(self.b,self.c) / self.V  # 1/Angstrom
        self.dkb = 2 * np.pi * np.cross(self.c,self.a) / self.V  # 1/Angstrom
        self.dkc = 2 * np.pi * np.cross(self.a,self.b) / self.V  # 1/Angstrom

        # reciprocal vectors
        self.ka = self.dka * self.nx  # 1/Angstrom
        self.kb = self.dkb * self.ny  # 1/Angstrom
        self.kc = self.dkc * self.nz  # 1/Angstrom

        # reciprocal lattice parameter matrix (1/Angstrom) - first row is ka, second row is kb, third row is kc
        self.B = np.vstack((self.ka, self.kb, self.kc))

        if verbose:
            print('A (Angstrom)\n', self.A)
            print('\nB (1/Angstrom)\n', self.B)


        # need to convert also units of the scalar field >>contained<< in the numpy array
        # spin density in units of Bohr magnetons per Angstrom^3 ??


        # ================== FFT ==================

        # norm='backward' means no prefactor applied
        self.F = fft.fftshift(fft.fftn(self.array, norm='backward'))
        self.F_abs_sq = np.square(np.abs(self.F))


        # ================== PLOTTING ==================

        # ----------------- REAL SPACE -----------------
        if plot_real_space_spin_density:
            self.plot_cube_file(cube_data)


        # ----------------- RECIPROCAL SPACE -----------------
        # sum all projections into plane (defined by a vector normal to the plane)
        # n_vec_plane = np.array([0, 0, 1])

        #    SIMPLE FIRST: just assume c is along z and sum along c axis
        # sum along c
        # !!!!! stupid coordinate system of MnGeO4 -- need to sum along x axis

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # !!!!!!!!!!!!!!!!!!!!!!!1 sum along z for the next
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # find a momentum along z !!! 
        # along a for now change <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        kc_array_max = abs(self.ka[2])
        print(f'kz_target must be between 0.00 and {kc_array_max:.6f} 1/Angstrom')
        if kz_target > kc_array_max:
            raise ValueError(f'kz_target must be between 0.00 and {kc_array_max:.6f} 1/Angstrom')
        kc_array = np.linspace(0, kc_array_max, self.nx)
        i_kz = np.argmin(np.abs(kc_array - kz_target))

        # F_abs_sq_sum_a = np.sum(F_abs_sq, axis=0)
        
        for i_kz in range(0, self.nx):
            F_abs_sq_sum_a = self.F_abs_sq[i_kz, :, :]
            self.plot_2D_fft(F_abs_sq_sum_a, k1=self.kb, k2=self.kc, fout_name=f'./Mn2GeO4_kz_tomography/log_scale/F_abs_squared_log-scale_kz_at_index_{i_kz}.png')


    def plot_cube_file(self, fout_name='rho_sz.png'):
        # plot numpy array as a scalar field
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Create a 3D grid
        x = np.arange(self.nx)
        y = np.arange(self.ny)
        z = np.arange(self.nz)
        X, Y, Z = np.meshgrid(x, y, z)

        # Plot the scalar field
        plot = ax.scatter(X, Y, Z, c=self.array.flatten(), cmap='viridis')
        ax.set_aspect('equal', adjustable='box')
        # plot colorbar the colorbar
        plt.colorbar(plot)
        plt.tight_layout()
        plt.savefig(fout_name)
        plt.close()


    def plot_2D_fft(self, scalar_2D_array, k1=[1,0,0], k2=[-0.5,1,0], fout_name='colormap_2D_out.png', verbose=True):

        i_relevant = 0
        j_relevant = 1


        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)

        k1 = np.array(k1)
        k2 = np.array(k2)

        na = scalar_2D_array.shape[0]
        nb = scalar_2D_array.shape[1]

        if verbose:
            print(k1, k2)
            print(scalar_2D_array.shape)

        # 2D grid with correct units but no dimensionality
        i_vals = (np.arange(na)-na//2) / na
        j_vals = (np.arange(nb)-nb//2) / nb
        I, J = np.meshgrid(i_vals, j_vals, indexing='ij')

        # Compute the actual coordinates in 2D space
        X = I * k1[i_relevant] + J * k2[i_relevant]
        Y = I * k1[j_relevant] + J * k2[j_relevant]    

        plt.pcolormesh(X, Y, np.log(np.abs(scalar_2D_array)), shading='auto', cmap='viridis', )

        # colorbar
        plt.colorbar(label="log(|F|^2)_xy")

        # Overlay grid points
        # plt.scatter(X, Y, color='black', s=1)

        # plot lattice vectors
        head_width = 0.04*np.linalg.norm(k1)
        head_length = 0.07*np.linalg.norm(k1)
        arrow_line_color = 'k'
        linestyle = (5, (5, 5))
        ax.arrow(0, 0, k1[0], k1[1], head_width=head_width, head_length=head_length, fc=arrow_line_color, ec=arrow_line_color)
        ax.arrow(0, 0, k2[0], k2[1], head_width=head_width, head_length=head_length, fc=arrow_line_color, ec=arrow_line_color)
        ax.arrow(-k1[0]/2-k2[0]/2, -k1[1]/2-k2[1]/2, k1[0], k1[1], head_width=0, head_length=0, fc=arrow_line_color, ec=arrow_line_color, linestyle=linestyle)
        ax.arrow(-k1[0]/2-k2[0]/2, -k1[1]/2-k2[1]/2, k2[0], k2[1], head_width=0, head_length=0, fc=arrow_line_color, ec=arrow_line_color, linestyle=linestyle)
        ax.arrow(k1[0]/2-k2[0]/2, k1[1]/2-k2[1]/2, k2[0], k2[1], head_width=0, head_length=0, fc=arrow_line_color, ec=arrow_line_color, linestyle=linestyle)
        ax.arrow(-k1[0]/2+k2[0]/2, -k1[1]/2+k2[1]/2, k1[0], k1[1], head_width=0, head_length=0, fc=arrow_line_color, ec=arrow_line_color, linestyle=linestyle)


        # Formatting
        plt.xlabel(r"$k_x$ ($\mathrm{\AA}^{-1}$)", fontsize=12)
        plt.ylabel(r"$k_y$ ($\mathrm{\AA}^{-1}$)", fontsize=12)
        plt.title("F")

        plt.axis('equal')  # Keep aspect ratio
        plt.tight_layout()
        plt.savefig(fout_name, dpi=400)
        plt.close()


if __name__ == '__main__':

    kz = 3

    density = Density(kz_target=kz, verbose=True, plot_real_space_spin_density=False)

    # twoD_data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    # plot_2D_fft(twoD_data)