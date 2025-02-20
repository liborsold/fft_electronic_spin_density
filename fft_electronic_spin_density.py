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

    def __init__(self, verbose=True):
        """_summary_

        Args:
            kz (int, optional): The selected cut at k_z in 1/Angstrom. Defaults to 5.
            verbose (bool, optional): _description_. Defaults to True.
            plot_real_space_spin_density (bool, optional): _description_. Defaults to False.
        """

        # ================== READ CUBE FILE ==================

        # cube[0] is the scalar field numpy array
        cube_data = read_cube(fin_cube)

        cube_data = self.coordinate_permutation(cube_data, permutation=[2,1,0])

        # cube[1] contains dictionary with metadata - keys are 'org', 'xvec', 'yvec', 'zvec', 'atoms'
        # get unit cell size: 'xvec' gives 
        #    the >>spacing<< between grid points 
        #     along the first index of the array (--> a lattice vector of the unit celll) 
        #      in units of Bohr radii
        # ... being used below
        
        # numpy array dimensions
        self.array = cube_data[0]
        self.na, self.nb, self.nc = self.array.shape


        # ================== UNITS ==================

        # unit conversion: see http://publish.illinois.edu/yubo-paul-yang/tutorials/quantum-espresso/understand-fast-fourier-transform/
        #    and   https://en.wikipedia.org/wiki/Reciprocal_lattice (formulas for reciprocal lattice vectors in 3D)

        #--- REAL SPACE ---
        # real-space grid spacing in Angstrom
        self.da = np.array(list(cube_data[1]['xvec'])) * physical_constants['Bohr radius'][0] * 1e10  # Angstrom
        self.db = np.array(list(cube_data[1]['yvec'])) * physical_constants['Bohr radius'][0] * 1e10  # Angstrom
        self.dc = np.array(list(cube_data[1]['zvec'])) * physical_constants['Bohr radius'][0] * 1e10  # Angstrom

        # real-space lattice vectors (nx, ny, nz are the number of grid points in each direction - the dimensions of the cube numpy array)
        self.a =  self.da * self.na # Angstrom
        self.b =  self.db * self.nb # Angstrom
        self.c =  self.dc * self.nc # Angstrom

        # volume of the unit cell
        self.V = np.dot(self.a, np.cross(self.b, self.c))  # Angstrom^3

        # lattice parameter matrix (Angstrom) - should match the lattice parameters in scf.in file - first row is a, second row is b, third row is c
        self.A = np.vstack((self.a, self.b, self.c))

        # make a grid in x,y,z cartesian coordinates of the real-space grid points
        a_idx = np.arange(self.na) / self.na
        b_idx = np.arange(self.nb) / self.nb
        c_idx = np.arange(self.nc) / self.nc

        # mesh arrays are 3D arrays
        a_idx_mesh, b_idx_mesh, c_idx_mesh = np.meshgrid(a_idx, b_idx, c_idx, indexing='ij')
        
        # flatten them into (nx*ny*nz, 3) array of reciprocal coordinates
        a_idx_mesh_flat = a_idx_mesh.flatten()
        b_idx_mesh_flat = b_idx_mesh.flatten()
        c_idx_mesh_flat = c_idx_mesh.flatten()

        # check that we can reshape back to the original shape
        # a_idx_mesh_reshaped = a_idx_mesh_flat.reshape((self.na, self.nb, self.nc))
        # assert np.allclose(a_idx_mesh, a_idx_mesh_reshaped), 'Meshed arrays are not the same after flattening and reshaping back again!'

        # concatenate the flattened arrays
        r_rec_mesh_flat = np.vstack((a_idx_mesh_flat, b_idx_mesh_flat, c_idx_mesh_flat)).T

        # convert to cartesian coordinates
        r_cart_mesh_flat = r_rec_mesh_flat @ self.A

        # cartesian coordinates 3D mesh arrays - ready for plotting (in Angstrom)
        self.x_cart_mesh = r_cart_mesh_flat[:, 0].reshape((self.na, self.nb, self.nc))
        self.y_cart_mesh = r_cart_mesh_flat[:, 1].reshape((self.na, self.nb, self.nc))
        self.z_cart_mesh = r_cart_mesh_flat[:, 2].reshape((self.na, self.nb, self.nc))

        #--- RECIPROCAL SPACE ---

        # get reciprocal lattice spacings
        self.dka = 2 * np.pi * np.cross(self.b,self.c) / self.V  # 1/Angstrom
        self.dkb = 2 * np.pi * np.cross(self.c,self.a) / self.V  # 1/Angstrom
        self.dkc = 2 * np.pi * np.cross(self.a,self.b) / self.V  # 1/Angstrom

        # reciprocal vectors
        self.ka = self.dka * self.na  # 1/Angstrom
        self.kb = self.dkb * self.nb  # 1/Angstrom
        self.kc = self.dkc * self.nc  # 1/Angstrom

        # reciprocal lattice parameter matrix (1/Angstrom) - first row is ka, second row is kb, third row is kc
        self.B = np.vstack((self.ka, self.kb, self.kc))

        if verbose:
            print('A (Angstrom)\n', self.A)
            print('\nB (1/Angstrom)\n', self.B)


        # need to convert also units of the scalar field >>contained<< in the numpy array
        # spin density in units of Bohr magnetons per Angstrom^3 ??


    def coordinate_permutation(self, cube_data, permutation=[2,1,0]):
        """Swap in a cyclic way (x,y,z) -> (y,z,x) -> (z,x,y) depending on number of steps (1 or 2).

        Args:
            cube_data (_type_): _description_
            steps (int, optional): _description_. Defaults to 1.
        """
        array = cube_data[0]
        xvec = cube_data[1]['xvec']
        yvec = cube_data[1]['yvec']
        zvec = cube_data[1]['zvec']

        array = np.moveaxis(array, [0, 1, 2], permutation)
        
        xyz_vec = list(np.array([xvec, yvec, zvec], dtype=object)[permutation])
    
        cube_data[1]['xvec'] = xyz_vec[0]
        cube_data[1]['yvec'] = xyz_vec[1]
        cube_data[1]['zvec'] = xyz_vec[2]

        return (array, cube_data[1])
        

    def plot_cube_file_outer_surface(self, fout_name='rho_sz.png'):
        # plot numpy array as a scalar field
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Create a 3D grid
        x = np.arange(self.na)
        y = np.arange(self.nb)
        z = np.arange(self.nc)
        X, Y, Z = np.meshgrid(x, y, z)

        # Plot the scalar field
        plot = ax.scatter(X, Y, Z, c=self.array.flatten(), cmap='viridis')
        ax.set_aspect('equal', adjustable='box')
        # plot colorbar the colorbar
        plt.colorbar(plot)
        plt.tight_layout()
        plt.savefig(fout_name)
        plt.close()


    def plot_cube_file(self, c_idx_arr=[0,1,-1], fout_name='rho_sz.png'):
        """For an array of inices, plot a 2D map as contourf at that z index of the 3D scalar field into a 3D plot at the height given by the z value.

        Args:
            c_idx_arr (list, optional): _description_. Defaults to [0,1,-1].
            fout_name (str, optional): _description_. Defaults to 'rho_sz.png'.
        """
        scale_down_data = 0.0001

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Create a 3D grid
        X = self.x_cart_mesh[:,:,0]
        Y = self.y_cart_mesh[:,:,0]
        z_cart = np.arange(self.nc)/self.nc * self.c[2]

        z_max_abs = np.max(np.abs(self.array))*scale_down_data

        # Plot the scalar field
        for c_idx in c_idx_arr:
            z_curr = z_cart[c_idx]
            Z_arr = self.array[:, :, c_idx]
            levels = np.linspace(np.min(Z_arr), np.max(Z_arr), 100)*scale_down_data
            plot = ax.contourf(X, Y, z_curr+Z_arr*scale_down_data, cmap='coolwarm', zdir='z', levels=z_curr+levels, vmin=z_curr-z_max_abs, vmax=z_curr+z_max_abs)

        # get the colorbar and adjust its limits
        cbar = plt.colorbar(plot)
        # cbar.set_clim(-z_max_abs, z_max_abs)
        # plt.colorbar(plot)

        # color
        
        # margin = 0.05
        # plt.xlim(0, np.max(self.x_cart_mesh*(1+margin)))
        # plt.ylim(0, np.max(self.y_cart_mesh*(1+margin)))
        # ax.set_zlim(0, np.max(self.z_cart_mesh))

        # ax.set_aspect('equal', adjustable='box')
        # plot colorbar the colorbar

        # plt.tight_layout()
        plt.savefig(fout_name)
        plt.close()

    def FFT(self, verbose=True):
        # norm='backward' means no prefactor applied
        self.F = fft.fftshift(fft.fftn(self.array, norm='backward'))
        self.F_abs_sq = np.square(np.abs(self.F))

    def get_i_kz(self, kz_target):
                #    SIMPLE FIRST: just assume c is along z and sum along c axis
        # sum along c
        # !!!!! stupid coordinate system of MnGeO4 -- need to sum along x axis

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # !!!!!!!!!!!!!!!!!!!!!!!1 sum along z for the next
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # find a momentum along z !!! 
        # along a for now change <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        kc_array_max = abs(self.kc[2])
        print(f'- by the way, kz_target must be between 0.00 and {kc_array_max:.6f} 1/Angstrom')
        if kz_target > kc_array_max:
            raise ValueError(f'kz_target must be between 0.00 and {kc_array_max:.6f} 1/Angstrom')
        kc_array = np.linspace(0, kc_array_max, self.nc)
        i_kz = np.argmin(np.abs(kc_array - kz_target))
        return i_kz


    def plot_2D_fft(self, i_kz, k1_idx=0, k2_idx=1, fout_name='colormap_2D_out.png', verbose=True):

        # ----------------- RECIPROCAL SPACE PLOTTING -----------------
        # sum all projections into plane (defined by a vector normal to the plane)
        # n_vec_plane = np.array([0, 0, 1])

        n1 = np.array((self.na, self.nb, self.nc))[k1_idx]
        n2 = np.array((self.na, self.nb, self.nc))[k2_idx]

        take_idx = [0, 1, 2]
        take_idx.remove(k1_idx)
        take_idx.remove(k2_idx)
        take_idx = take_idx[0]
        print('take_idx', take_idx)

        # PREPARE 2D array 
        #    - sum?
        # F_abs_sq_sum_a = np.sum(F_abs_sq, axis=take_idx)

        #    - cut
        F_abs_sq_cut = self.F_abs_sq.take(i_kz, axis=take_idx)

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)

        k1 = [self.ka, self.kb, self.kc][k1_idx]
        k2 = [self.ka, self.kb, self.kc][k2_idx]

        if verbose:
            print(k1, k2)
            print(self.array.shape)

        # 2D grid with correct units but no dimensionality
        i_vals = (np.arange(n1)-n1//2) / n1
        j_vals = (np.arange(n2)-n2//2) / n2
        I, J = np.meshgrid(i_vals, j_vals, indexing='ij')

        # Compute the actual coordinates in 2D space
        X = I * k1[0] + J * k2[0]
        Y = I * k1[1] + J * k2[1]  

        plt.pcolormesh(X, Y, np.log(np.abs(F_abs_sq_cut)), shading='auto', cmap='viridis', )

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

    density = Density(verbose=True)

    for i in np.arange(0, density.nc, 4):
        c_idx_array = np.array([i, 124]) #,11,12])+25 #np.arange(0, density.nc, density.nc//4)
        density.plot_cube_file(c_idx_arr=c_idx_array, fout_name=f'./rho_sz_exploded_{i}.png')
    exit()
    density.FFT(verbose=True)

    # kz = 30
    # i_kz = density.get_i_kz(kz_target=kz)
    # density.plot_2D_fft(i_kz=i_kz, k1_idx=k1_idx, k2_idx=k2_idx, fout_name=f'./test_fft.png')

    for i_kz in range(0, density.nc):
        density.plot_2D_fft(i_kz=i_kz, fout_name=f'./Mn2GeO4_kz_tomography/log_scale/F_abs_squared_log-scale_kz_at_index_{i_kz}.png')

    # twoD_data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    # plot_2D_fft(twoD_data)