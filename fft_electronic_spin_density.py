# include cubetools.py in folder utils/

from utils.cubetools import read_cube, write_cube
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy.fft as fft
# import bohr radius
from scipy.constants import physical_constants
import matplotlib as mpl
from matplotlib import pylab
from matplotlib.colors import ListedColormap
import os
from matplotlib.ticker import FuncFormatter
from copy import deepcopy

class Density:
    """Read, visualize and fourier transform (spin) density from gaussian .cube files.
        Replace by a model function if required.
    """

    def __init__(self, fname_cube_file='./seedname.cube', permutation=None, verbose=True):
        """_summary_

        Args:
            kz (int, optional): The selected cut at k_z in 1/Angstrom. Defaults to 5.
            verbose (bool, optional): _description_. Defaults to True.
            plot_real_space_spin_density (bool, optional): _description_. Defaults to False.
        """

        # ================== READ CUBE FILE ==================

        # cube[0] is the scalar field numpy array
        cube_data = read_cube(fname_cube_file)
        if permutation:
            cube_data = self.coordinate_permutation(cube_data, permutation=permutation)

        # cube[1] contains dictionary with metadata - keys are 'org', 'xvec', 'yvec', 'zvec', 'atoms'
        # get unit cell size: 'xvec' gives 
        #    the >>spacing<< between grid points 
        #     along the first index of the array (--> a lattice vector of the unit celll) 
        #      in units of Bohr radii
        # ... being used below

        # metadata
        self.metadata = cube_data[1]
        
        # numpy array dimensions
        self.array = cube_data[0]
        self.na, self.nb, self.nc = self.array.shape

        self.metadata_orig = deepcopy(self.metadata)
        self.array_orig = deepcopy(self.array)


        # ================== UNITS ==================

        # unit conversion: see http://publish.illinois.edu/yubo-paul-yang/tutorials/quantum-espresso/understand-fast-fourier-transform/
        #    and   https://en.wikipedia.org/wiki/Reciprocal_lattice (formulas for reciprocal lattice vectors in 3D)

        #--- REAL SPACE ---
        # real-space grid spacing in Angstrom
        self.da = np.array(list(cube_data[1]['xvec'])) * physical_constants['Bohr radius'][0] * 1e10  # Angstrom
        self.db = np.array(list(cube_data[1]['yvec'])) * physical_constants['Bohr radius'][0] * 1e10  # Angstrom
        self.dc = np.array(list(cube_data[1]['zvec'])) * physical_constants['Bohr radius'][0] * 1e10  # Angstrom

        if verbose: print('da, db, dc', self.da, self.db, self.dc, 'Angstrom')

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

        # cartesian coordinates 3D mesh arrays - ready for use in plotting (in Angstrom)
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

        if verbose: print('dka, dkb, dkc', self.dka, self.dkb, self.dkc, '1/Angstrom')

        # reciprocal lattice parameter matrix (1/Angstrom) - first row is ka, second row is kb, third row is kc
        self.B = np.vstack((self.ka, self.kb, self.kc))

        if verbose:
            print('A (Angstrom)\n', self.A)
            print('\nB (1/Angstrom)\n', self.B)


        # need to convert also units of the scalar field >>contained<< in the numpy array
        # spin density in units of Bohr magnetons per Angstrom^3 ??

    def get_sites_of_atoms(self, site_idx):
        """Return the site centers of the atoms at the given indices.

        Args:
            site_idx (list of integers): integers of the required sites

        Returns:
            list of tuples: tuples are cartesian coordinates (xyz) in Angstrom of the site centers
        """
        site_centers = []
        for idx in site_idx:
            # self.metadata['atoms'][idx] is a tuple with the atomic mass and the coordinates of the atom in the unit cell
            #   coordinates is a map object -> convert to list: first element is again the atomic mass -> take the rest: units are in Bohr radii -> convert to Angstrom -> convert to tuple -> add to list
            site_centers.append(tuple(np.array(list(self.metadata['atoms'][idx][1])[1:])*physical_constants['Bohr radius'][0]*1e10))
        return site_centers

    def mask_except_sites(self, leave_sites):
        # initialize true array
        mask = np.zeros_like(self.x_cart_mesh, dtype=bool)
        for center, radius in zip(leave_sites['site_centers'], leave_sites['site_radii']):
            mask_i = np.sqrt((self.x_cart_mesh - center[0])**2 + (self.y_cart_mesh - center[1])**2 + (self.z_cart_mesh - center[2])**2) < radius
            mask = np.logical_or(mask, mask_i)
        self.array[~mask] = 0

    def get_kz_at_index(self, kz_index=30):
        """Return the k_z value (in Angstrom^-1) at a given index: first or last indices are -+ k_max/2, the middle index is k_z=0 (data is zero-centered; fftshifted after fft was performed).
        Check that index is in range.
        """
        assert kz_index < self.nc, f'kz_index must be between 0 and {self.nc-1} (inclusive)'
        kz = (kz_index - self.nc//2) * self.dkc[2]
        print(f'kz at index {kz_index} is {kz:.6f} 1/Angstrom')
        return kz
    

    def get_index_at_kz(self, kz_target=15):
        """Return the index at a given k_z value (in Angstrom^-1). Check that the k_z value is in range.

        Args:
            kz_target (int, optional): _description_. Defaults to 15.
        """
        if kz_target > np.abs(self.kc[2]):
            raise ValueError(f'kz_target must be between 0.00 and {np.abs(self.kc[2]):.6f} 1/Angstrom')
        i_kz = np.argmin(np.abs((np.arange(self.nc) - self.nc//2) * self.dkc[2] - kz_target))
        print(f'index for kz_target {(i_kz - self.nc//2) * self.dkc[2]:.6f} 1/Angstrom is {i_kz}')
        return i_kz


    def replace_by_model(self, type='gaussian', parameters={'sigmas':[0.5], 'centers':[(0.5, 0.5, 0.5)], 'signs':[1]}):
        """Replace the scalar field in the numpy array by a model function.

        Args:
            type (str, optional): Type of the model function. Defaults to 'gaussian'.
            parameters (dict, optional): Parameters of the model function. Defaults to {'sigmas':[0.5], 'centers':[(0.5, 0.5, 0.5)], 'signs':[1]}.
        """

        # create models
        models = {}
        def gaussian(x, y, z, sigma=0.5, center=(3,3,3), sign=1):
            """Gaussian distribution in 3D space - https://math.stackexchange.com/questions/434629/3-d-generalization-of-the-gaussian-point-spread-function

            Args:
                x (_type_): Cartesian x coordinate in Angstrom.
                y (_type_): Cartesian y coordinate in Angstrom.
                z (_type_): Cartesian z coordinate in Angstrom.
                sigma (float, optional): _description_. Defaults to 0.5.
                center (tuple, optional): _description_. Defaults to (3,3,3).
                sign (int, optional): _description_. Defaults to 1.

            Returns:
                _type_: _description_
            """
            # normalized gaussian function - see 
            return sign * np.exp(-((x-center[0])**2 + (y-center[1])**2 + (z-center[2])**2)/(2*sigma**2)) # * 1/(sigma**3 * (2*np.pi)**(3/2)) 
        models['gaussian'] = gaussian

        def dz2(x, y, z, sigma=0.5, center=(3,3,3), sign=1):
            """dz2 orbital distribution in 3D space - https://math.stackexchange.com/questions/434629/3-d-generalization-of-the-gaussian-point-spread-function

            Args:
                x (_type_): Cartesian x coordinate in Angstrom.
                y (_type_): Cartesian y coordinate in Angstrom.
                z (_type_): Cartesian z coordinate in Angstrom.
                sigma (float, optional): _description_. Defaults to 0.5.
                center (tuple, optional): _description_. Defaults to (3,3,3).
                sign (int, optional): _description_. Defaults to 1.

            Returns:
                _type_: _description_
            """
            # normalized gaussian function - see
            r2 = (x-center[0])**2 + (y-center[1])**2 + (z-center[2])**2
            return sign * np.abs(3*(z-center[2])**2 - r2)/r2 * np.exp(-(r2/(2*sigma**2)))
        models['dz2'] = dz2
        
        # choose the 3D scalar field function from models
        f = models[type]

        # construct an equivalent of the scalar array loaded from the .cube file but
        #   feeding in the model
        model_density = np.zeros_like(self.array)

        # place all the e.g. gaussians in the space
        for i in range(len(parameters['sigmas'])):
            # create a function in 3D space that gives a Gaussian density distribution around point centers[i] with standard deviation sigmas[i]
            # the sign of the density is given by signs[i]
            sigma = parameters['sigmas'][i]
            center = parameters['centers'][i]
            sign = parameters['signs'][i]

            model_density += f(self.x_cart_mesh, self.y_cart_mesh, self.z_cart_mesh, sigma=sigma, center=center, sign=sign)
        
        # replace the scalar field in the numpy array by the model
        self.array = model_density


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


    # def plot_cube_file_original(self, c_idx_arr=[0,1,-1], fout_name='rho_sz.png', alpha=0.2, figsize=(8.0, 6), dpi=300, zeros_transparent=True,
    #                    xlims=None, ylims=None, zlims=None, show_plot=False):
    #     """THE ORIGINAL
        
    #     For an array of indices, plot a 2D map as contourf at that z index of the 3D scalar field into a 3D plot at the height given by the z value.

    #     Args:
    #         c_idx_arr (list, optional): show cuts at these indeces. Defaults to [0,1,-1].
    #         fout_name (str, optional): _description_. Defaults to 'rho_sz.png'.
    #     """
    #     scale_down_data = 0.02

    #     if zeros_transparent:
    #         transparent_sigma = 0.15
    #         alpha_baseline = 0.50
    #         print(f"Plotting with transparency near the middle of the colormap: alpha = {alpha_baseline:.3f} (1 - exp(-(x/sigma)^2) with sigma={transparent_sigma:.3f})")
    #         x = np.linspace(-0.5, 0.5, pylab.cm.coolwarm.N)
    #         alpha = alpha_baseline*(1 - np.exp(-(x/transparent_sigma)**2))
    #         plt.plot(x,alpha)
    #         plt.xlabel('Transparency')
    #         plt.ylabel('rho_sz')
    #         plt.savefig('/'.join(fout_name.split('/')[:-1]) + '/transparency_profile_rho_sz_all-in-one.png', dpi=400)
    #         plt.close()
    #     else:
    #         # uniform 
    #         alpha = np.ones(pylab.cm.coolwarm.N)*alpha

    #     my_cmap = pylab.cm.coolwarm(np.arange(pylab.cm.coolwarm.N))
    #     my_cmap[:,-1] = alpha
    #     cmap = ListedColormap(my_cmap)

    #     fig = plt.figure(figsize=figsize)
    #     ax = fig.add_subplot(111, projection='3d')

    #     # Create a 3D grid
    #     X = self.x_cart_mesh[:,:,0] # !!!(1)
    #     Y = self.y_cart_mesh[:,:,0]# !!!(2)
    #     z_cart = np.arange(self.nc)/self.nc * self.c[2] # !!!(3)

    #     z_max_abs_unscaled = np.max(np.abs(self.array))
    #     z_max_abs = z_max_abs_unscaled*scale_down_data

    #     # Plot the scalar field
    #     for c_idx in c_idx_arr:
    #         z_curr = z_cart[c_idx]
    #         Z_arr = self.array[:, :, c_idx]
            
    #         # get levels
    #         min_Z_arr = np.min(Z_arr)
    #         max_Z_arr = np.max(Z_arr)
    #         if abs(min_Z_arr - max_Z_arr) < 1e-10:
    #             levels = np.linspace(min_Z_arr-1e-10, max_Z_arr+1e-10, 100)*scale_down_data
    #         else:
    #             levels = np.linspace(min_Z_arr, max_Z_arr, 100)*scale_down_data
    #         ax.contourf(X, Y, z_curr+Z_arr*scale_down_data, cmap=cmap, zdir='z', levels=z_curr+levels, vmin=z_curr-z_max_abs, vmax=z_curr+z_max_abs)

    #     fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-z_max_abs_unscaled, z_max_abs_unscaled), cmap=cmap),
    #          ax=ax, orientation='vertical', label='spin density')
    #     # manually make a colorbar with limits -z_max_abs, z_max_abs and cmap coolwarm
        
    #     # cbar.set_clim(-z_max_abs, z_max_abs)
    #     # plt.colorbar(plot)

    #     # color
        
    #     # margin = 0.05
    #     # plt.xlim(0, np.max(self.x_cart_mesh*(1+margin)))
    #     # plt.ylim(0, np.max(self.y_cart_mesh*(1+margin)))

    #     plt.title(f'Spin density ranging from {np.min(self.array):.3f} to {np.max(self.array):.3f}')

    #     ax.set_zlim(min(0, self.c[2]), max(0, self.c[2]))

    #     ax.set_xlabel(r'$x$ ($\mathrm{\AA}$)', fontsize=11)
    #     ax.set_ylabel(r'$y$ ($\mathrm{\AA}$)', fontsize=11)
    #     ax.set_zlabel(r'$z$ ($\mathrm{\AA}$)', fontsize=11)

    #     if xlims:
    #         ax.set_xlim(xlims)
    #     if ylims:
    #         ax.set_ylim(ylims)
    #     if zlims:
    #         ax.set_zlim(zlims)

    #     ax.set_aspect('equal', adjustable='box')
    #     # plot colorbar the colorbar
    #     plt.tight_layout()

    #     if show_plot:
    #         plt.show()
    #     else:
    #         plt.savefig(fout_name, dpi=dpi)
    #         plt.close()


    def plot_cube_file_general(self, X, Y, z_levels_cart, scalar3D_data, c_idx_arr=[0,1,-1], fout_name='rho_sz.png', alpha=0.2, figsize=(8.0, 6), dpi=300, zeros_transparent=True,
                       xlims=None, ylims=None, zlims=None, show_plot=False, xlabel=None, ylabel=None, zlabel=None, colors_centered=True, cmap='coolwarm', alpha_baseline = 0.50, transparent_sigma=0.15, 
                       colorbar_label='spin density'):
        """For an array of indices, plot a 2D map as contourf at that z index of the 3D scalar field into a 3D plot at the height given by the z value.

        Args:
            c_idx_arr (list, optional): show cuts at these indeces. Defaults to [0,1,-1].
            fout_name (str, optional): _description_. Defaults to 'rho_sz.png'.
        """
        scale_down_data = 0.02 * 1/np.max(np.abs(scalar3D_data))

        cmap = plt.get_cmap(cmap)
        if zeros_transparent:        
            if colors_centered:
                t = np.linspace(-0.5, 0.5, cmap.N)
                alpha = alpha_baseline*(1 - np.exp(-(t/transparent_sigma)**2))
                print(f"Plotting with transparency near the middle of the colormap: alpha = {alpha_baseline:.3f} (1 - exp(-(x/sigma)^2) with sigma={transparent_sigma:.3f})")
            else:
                t = np.linspace(0, 1, cmap.N)
                alpha = alpha_baseline * np.exp(-((1-t)/transparent_sigma)**3)
                print(f"Plotting with transparency everywhere except around largest values: alpha = {alpha_baseline:.3f} exp(-((1-x)/sigma)^3) with sigma={transparent_sigma:.3f})")
            
            plt.plot(t,alpha)
            plt.xlabel('Transparency')
            plt.ylabel('rho_sz')
            fsplit = fout_name.split('/')
            transparency_fout = '/'.join(fsplit[:-1]) + '/transparency_profile_' + fsplit[-1]
            plt.savefig(transparency_fout, dpi=400)
            plt.close()
        else:
            # uniform 
            alpha = np.ones(cmap.N)*alpha

        my_cmap = cmap(np.arange(cmap.N))
        my_cmap[:,-1] = alpha
        cmap = ListedColormap(my_cmap)

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')

        # Create a 3D grid

        z_max_abs_unscaled = np.max(np.abs(scalar3D_data))
        z_max_abs = z_max_abs_unscaled*scale_down_data

        # Plot the scalar field
        for c_idx in c_idx_arr:
            z_curr = z_levels_cart[c_idx]
            Z_arr = scalar3D_data[:, :, c_idx]
            
            # get levels
            min_Z_arr = np.min(Z_arr)
            max_Z_arr = np.max(Z_arr)
            if abs(min_Z_arr - max_Z_arr) < 1e-10:
                # if all levels in the array are close to a single value (zero typically)
                levels = np.linspace(min_Z_arr-(1e-10/scale_down_data), max_Z_arr+(1e-10/scale_down_data), 100)*scale_down_data
            else:
                levels = np.linspace(min_Z_arr, max_Z_arr, 100)*scale_down_data
            
            if np.max(levels) < 1e-10:
                levels *= 1e-10 / np.max(np.abs(levels))

            ax.contourf(X, Y, z_curr+Z_arr*scale_down_data, cmap=cmap, zdir='z', levels=z_curr+levels, vmin=z_curr-z_max_abs, vmax=z_curr+z_max_abs)

        if colors_centered:
            colorbar_min = -z_max_abs_unscaled
            colorbar_max = z_max_abs_unscaled
        else:
            colorbar_min = 0
            colorbar_max = z_max_abs_unscaled
        fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(colorbar_min, colorbar_max), cmap=cmap),
             ax=ax, orientation='vertical', label=colorbar_label)
        # manually make a colorbar with limits -z_max_abs, z_max_abs and cmap coolwarm
        
        # cbar.set_clim(-z_max_abs, z_max_abs)
        # plt.colorbar(plot)

        # color
        
        # margin = 0.05
        # plt.xlim(0, np.max(self.x_cart_mesh*(1+margin)))
        # plt.ylim(0, np.max(self.y_cart_mesh*(1+margin)))

        plt.title(f'{colorbar_label} ranging from {np.min(scalar3D_data):.3f} to {np.max(scalar3D_data):.3f}')

        ax.set_zlim(min(0, self.c[2]), max(0, self.c[2]))  #!!!

        # axes labels
        if not xlabel:
            xlabel = r'$x$ ($\mathrm{\AA}$)'
        if not ylabel:
            ylabel = r'$y$ ($\mathrm{\AA}$)'
        if not zlabel:
            zlabel = r'$z$ ($\mathrm{\AA}$)'

        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_zlabel(zlabel, fontsize=11)


        if xlims:
            ax.set_xlim(min(xlims), max(xlims))
        if ylims:
            ax.set_ylim(min(ylims), max(ylims))
        if zlims:
            ax.set_zlim(min(zlims), max(zlims))

        ax.set_aspect('equal', adjustable='box')
        # plot colorbar the colorbar
        plt.tight_layout()

        if show_plot:
            plt.show()
        else:
            plt.savefig(fout_name, dpi=dpi)
            plt.close()


    def plot_cube_rho_sz(self, c_idx_arr=[0,1,-1], fout_name='rho_sz.png', alpha=0.2, figsize=(8.0, 6), dpi=300, zeros_transparent=True,
                       xlims=None, ylims=None, zlims=None, show_plot=False):
        
        """Concrete use of plot_cube_file_general for spin density files.
        """
        
        X = self.x_cart_mesh[:,:,0]
        Y = self.y_cart_mesh[:,:,0]
        z_levels_cart = np.arange(self.nc)/self.nc * self.c[2]
        scalar3D_data = self.array

        xlabel = r'$x$ ($\mathrm{\AA}$)'
        ylabel = r'$y$ ($\mathrm{\AA}$)'
        zlabel = r'$z$ ($\mathrm{\AA}$)'

        if not zlims:
            zlims = [min(0, self.c[2]), max(0, self.c[2])]

        self.plot_cube_file_general(X, Y, z_levels_cart, scalar3D_data, c_idx_arr=c_idx_arr, fout_name=fout_name, alpha=alpha, figsize=figsize, dpi=dpi, 
                                    zeros_transparent=zeros_transparent, xlims=xlims, ylims=ylims, zlims=zlims, show_plot=show_plot, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel)
        

    def plot_cube_fft(self, c_idx_arr=[0,1,-1], fout_name='rho_sz.png', alpha=0.2, figsize=(8.0, 6), dpi=300, zeros_transparent=True,
                       xlims=None, ylims=None, zlims=None, show_plot=False):
        
        """Concrete use of plot_cube_file_general for spin density files.
        """

        # 2D grid with correct units but no dimensionality
        i_vals = (np.arange(self.na)-self.na//2) / self.na
        j_vals = (np.arange(self.nb)-self.nb//2) / self.nb
        I, J = np.meshgrid(i_vals, j_vals, indexing='ij')

        X = I * self.ka[0] + J * self.kb[0]
        Y = I * self.ka[1] + J * self.kb[1]

        scalar3D_data = self.F_abs_sq

        z_levels_cart = (np.arange(self.nc)/self.nc - 0.5) * self.kc[2]
        xlabel = r'$k_x$ ($\mathrm{\AA}^{-1}$)'
        ylabel = r'$k_y$ ($\mathrm{\AA}^{-1}$)'
        zlabel = r'$k_z$ ($\mathrm{\AA}^{-1}$)'

        if not zlims:
            zlims = (min(z_levels_cart), max(z_levels_cart))

        self.plot_cube_file_general(X, Y, z_levels_cart, scalar3D_data, c_idx_arr=c_idx_arr, fout_name=fout_name, alpha=alpha, figsize=figsize, dpi=dpi, 
                                    zeros_transparent=zeros_transparent, xlims=xlims, ylims=ylims, zlims=zlims, show_plot=show_plot, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel, 
                                    colors_centered=False, cmap='viridis', 
                                    alpha_baseline=0.5,
                                    transparent_sigma=0.25, 
                                    colorbar_label=r'$|F|^2$')


    def FFT(self, verbose=True):
        # norm='backward' means no prefactor applied
        self.F = fft.fftshift(fft.fftn(self.array, norm='backward'), )
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


    def plot_fft_2D(self, i_kz, fft_as_log=False, k1_idx=0, k2_idx=1, fout_name='colormap_2D_out.png', verbose=True, figsize=(8.0, 6.0), 
                    dpi=500,
                    fixed_z_scale=True, 
                    xlims=None, ylims=None,zlims=None):

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

        fig = plt.figure(figsize=figsize)
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

        plot_array = np.abs(F_abs_sq_cut)
        if fft_as_log:
            plot_array = np.log(plot_array)
        plt.pcolormesh(X, Y, plot_array, shading='auto', cmap='viridis', )

        # colorbar
        label = r'$\mathrm{log}\(|F|^2\)_{xy}$' if fft_as_log else '$|F|^2$'
        def fmt(x, pos): 
            base, exponent = f"{x:2.1e}".split('e')
            exponent = f"{int(exponent):+01d}"  # Format exponent with a sign and 3 digits
            return f"{base}e{exponent}"
        # format string which keeps fixed length of the number, scientific format 
        plt.colorbar(label=label, format=FuncFormatter(fmt))

        if fixed_z_scale:
            plt.clim(0, np.max(self.F_abs_sq))

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
        plt.title(r"$|F|^2$($k_x$, $k_y$; $k_z$ = " + f'{self.get_kz_at_index(i_kz):.4f} '+r'$\mathrm{\AA}^{-1})$', fontsize=12)

        ax.set_aspect('equal', adjustable='box')  # Keep aspect ratio

        if xlims:
            plt.xlim(xlims)
        if ylims:
            plt.ylim(ylims)
        if zlims:
            plt.clim(zlims)

        plt.tight_layout()
        if fixed_z_scale:
            # add appendix to name (while keeping original file format)
            fsplit = fout_name.split('.')
            fout_name = '.'.join(fsplit[:-1]) + '_fix-scale.' + fsplit[-1]
        plt.savefig(fout_name, dpi=dpi)
        plt.close()

    def write_cube_file_rho_sz(self, fout='rho_sz_modified.cube'):
        """Write out the modified rho_sz to a cube file.

        Args:
            fout (str, optional): _description_. Defaults to 'rho_sz_modified.cube'.
        """
        write_cube(self.array, self.metadata_orig, fout)

    def write_cube_file_fft(self, fout='rho_sz_fft.cube'):
        """Write out the modified rho_sz to a cube file.

        Args:
            fout (str, optional): _description_. Defaults to 'rho_sz_modified.cube'.
        """
        meta_fft = deepcopy(self.metadata_orig)
        meta_fft['xvec'] = self.ka
        meta_fft['yvec'] = self.kb
        meta_fft['zvec'] = self.kc
        write_cube(self.F_abs_sq, meta_fft, fout)


def test_shift():
    array_3D = np.zeros((9,9,9), dtype=np.float_)
    idx_3D = np.zeros((9,9,9), dtype=np.bool_)

    for i in range(9):
        for j in range(9):
            for k in range(9):
                idx_3D[i,j,k] = i > 4 and j > 4 and k > 4
    # print(array3D)

    array_3D[idx_3D] = 1
    array_3D[~idx_3D] = -1

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # Create a 3D grid
    x = np.arange(9)
    y = np.arange(9)
    z = np.arange(9)
    X, Y, Z = np.meshgrid(x, y, z)
    # Plot the scalar field
    plot = ax.scatter(X, Y, Z, c=idx_3D.flatten(), cmap='coolwarm')
    ax.set_aspect('equal', adjustable='box')
    # plot colorbar the colorbar
    plt.colorbar(plot)
    plt.tight_layout()
    plt.savefig('fake_3D.png')
    plt.close()

    # shift
    array3D_shifted = np.fft.fftshift(idx_3D)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # Create a 3D grid
    x = np.arange(9)
    y = np.arange(9)
    z = np.arange(9)
    X, Y, Z = np.meshgrid(x, y, z)
    # Plot the scalar field
    plot = ax.scatter(X, Y, Z, c=array3D_shifted.flatten(), cmap='coolwarm')
    ax.set_aspect('equal', adjustable='box')
    # plot colorbar the colorbar
    plt.colorbar(plot)
    plt.tight_layout()
    plt.savefig('fake_3D_shifted.png')
    plt.close()


if __name__ == '__main__':

    # --- INPUT ----

    fname_cube_file = './cube_files/Cu2AC4_rho_sz_512.cube' #'./cube_files/Mn2GeO4_rho_sz.cube'
    
    permutation = None #!! for Mn2GeO4 need to use [2,1,0] to swap x,y,z -> z,y,x

    output_folder = './outputs/Cu2AC4/512/masked_Cu0' #_and_oxygens' # 'Mn2GeO4_kz_tomography_64' #'./gaussian/sigma_0.3_distance_1.0' # Mn2GeO4_kz_tomography_64

    # ---- CALCULATION CONTROL ----

    replace_DFT_by_model = False

    density_3D = False
    density_slices = False
    
    fft_3D = False
    full_range_fft_spectrum_cuts = False
    zoom_in_fft_spectrum_cuts = True

    write_cube_files = False

    # ---- PARAMETERS -----
    r_mt_Cu = 1.1 #Angstrom
    r_mt_O = 0.9 #Angstrom

    site_idx = [0] #, 1] # 16, 25, 9, 40, 1, 41, 8, 24, 17] #  16, 25, 9, 40] #[0]# None #[0,  16, 25, 9, 40, 1, 41, 8, 24, 17] #[1, 41, 8, 24, 17, 0,  16, 25, 9, 40] #[0] #, 16, 25, 9, 40]  # [0] #,  16, 25, 9, 40] #
    site_radii = [r_mt_Cu] # + 4*[r_mt_O] + [r_mt_Cu]+ 4*[r_mt_O]
    # !!!! [0, 3.5e6] for a single site and [0, 7e6] for double !
    fft_zlims = [0, 3.5e6] # arb. units

    density_figsize = (6.0, 4.5)
    dpi_rho = 500
    density_slice_each_n_images = 4

    all_in_one_xlims = (1.5, 6.5) #None
    all_in_one_ylims = (3.5, 8.5) #None
    all_in_one_zlims = (1.0, 7.5) #None

    # all_in_one_xlims = None
    # all_in_one_ylims = None
    # all_in_one_zlims = None

    fft_figsize = (4.5, 4.5)
    fft_dpi = 400
    fft_xlims = [-19, 19] # 1/Angstrom
    fft_ylims = [-19, 19] # 1/Angstrom
    
    fft_as_log = False
    fft_slice_each_n_images = 4
    all_in_one_density_total_slices = 300

    # -- create folder if does not exist --
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # ---- READ CUBE FILE -----
    density = Density(permutation=permutation, verbose=True, fname_cube_file=fname_cube_file)


    # ---- MASKING -----

    # leave_sites = {'site_centers':[(0.5, 0.5, 0.5)], 'site_radii':[0.5]}
    if site_idx and site_radii:
        site_centers = density.get_sites_of_atoms(site_idx)
        print('site_centers', site_centers)
        leave_sites = {'site_centers':site_centers, 'site_radii':site_radii}
        density.mask_except_sites(leave_sites)

    # get kz at index
    # density.get_kz_at_index(80)
    # density.get_index_at_kz(-14.67785)
    # density.get_index_at_kz(14.67785)

    # ---- INSERT MODEL -----
    # model_type = 'gaussian'
    model_type = 'dz2'
    
    # parameters = {'sigmas':[0.2, 0.2], 'centers':[(-3.0, -3, -5), (-2.0, -3, -5)], 'signs':[1, -1]}
    # parameters = {'sigmas':[0.2, 0.2], 'centers':[(-3.25, -3, -5), (-1.75, -3, -5)], 'signs':[1, -1]}
    # parameters = {'sigmas':[0.5, 0.5], 'centers':[(-3.00, -3, -5), (-2.00, -3, -5)], 'signs':[1, -1]}
    parameters = {'sigmas':[0.3, 0.3], 'centers':site_centers, 'signs':[1,-1]}

    if replace_DFT_by_model:
        density.replace_by_model(type=model_type, parameters=parameters)


    # ---- VISUALIZE DENSITY -----
    if density_slices:
        for i in np.arange(0, density.nc, density_slice_each_n_images):
            c_idx_array = np.array([i, 0]) #np.array([i, -1]
            density.plot_cube_rho_sz(c_idx_arr=c_idx_array, fout_name=f'{output_folder}/rho_sz_exploded_masked_{i}.jpg', alpha=0.8, figsize=density_figsize, dpi=dpi_rho, zeros_transparent=False)  # rho_sz_gauss_exploded

    if density_3D:
        c_idx_array = np.arange(0, density.nc, max(1, density.nc//all_in_one_density_total_slices))
        density.plot_cube_rho_sz(c_idx_arr=c_idx_array, fout_name=f'{output_folder}/rho_sz_exploded_masked_all.jpg', alpha=0.05, figsize=(5.5,5.5), dpi=dpi_rho, zeros_transparent=True,
                            xlims=all_in_one_xlims, ylims=all_in_one_ylims, zlims=all_in_one_zlims, show_plot=False)  # rho_sz_gauss_exploded_all

    # single cut
        # kz = 30
    # i_kz = density.get_i_kz(kz_target=kz)
    # density.plot_2D_fft(i_kz=i_kz, k1_idx=k1_idx, k2_idx=k2_idx, fout_name=f'./test_fft.png')

    # ---- WRITE MODIFIED DENSITY TO CUBE FILE -----
    if write_cube_files:
        density.write_cube_file_rho_sz(fout=f'{output_folder}/rho_sz_modified.cube')

    # ---- FFT -----
    density.FFT(verbose=True)

    # ---- WRITE MODIFIED FFT TO CUBE FILE -----
    if write_cube_files:
        density.write_cube_file_fft(fout=f'{output_folder}/fft.cube')

    # ---- VISUALIZE FFT -----
    
    # 3D
    if fft_3D:
        c_idx_array = np.arange(0, density.nc, max(1, density.nc//all_in_one_density_total_slices))
        xlims = [-9, 9] #None
        ylims = xlims
        zlims = xlims
        density.plot_cube_fft(c_idx_arr=c_idx_array, fout_name=f'{output_folder}/F_abs_sq_all.png', figsize=(5.5,5.5), dpi=dpi_rho, zeros_transparent=True,
                                xlims=xlims, ylims=ylims, zlims=zlims, show_plot=False)

    # (1) variable scale, full reciprocal space
    if full_range_fft_spectrum_cuts:
        for i_kz in range(0, density.nc, fft_slice_each_n_images):
            appendix = '_log' if fft_as_log else ''
            density.plot_fft_2D(i_kz=i_kz, fft_as_log=fft_as_log, 
                                fout_name=f'{output_folder}/F_abs_sq{appendix}-scale_kz_at_idx_{i_kz}.png', 
                                figsize=fft_figsize,
                                dpi=fft_dpi, 
                                fixed_z_scale=False,
                                xlims=None,
                                ylims=None)
        
    # (2) fixed scale, zoom-in
    if zoom_in_fft_spectrum_cuts:
        for i_kz in range(80, 121, 1):
            appendix = '_zoom_log' if fft_as_log else '_zoom'
            density.plot_fft_2D(i_kz=i_kz, fft_as_log=fft_as_log, 
                                fout_name=f'{output_folder}/F_abs_sq{appendix}-scale_kz_at_idx_{i_kz}.png', 
                                figsize=(5.5, 4.5),
                                dpi=fft_dpi,
                                fixed_z_scale=True,
                                xlims=fft_xlims,
                                ylims=fft_ylims, 
                                zlims=fft_zlims)

    # test_shift()
    # exit()

    # test plotting
    # twoD_data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    # plot_2D_fft(twoD_data)