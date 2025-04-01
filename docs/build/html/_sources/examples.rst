==========================
Examples
==========================

Example of use is given in ``./examples/example_of_use.ipynb`` in the `examples folder <https://github.com/liborsold/fft_electronic_spin_density/tree/master/examples>`_.

.. It uses the provided input files from the `cube_files folder <https://github.com/liborsold/fft_electronic_spin_density/tree/master/cube_files>`_.

.. fft_electronic_spin_density example image
.. figure::
   ./_static/images/example_of_use.png
   :width: 800px
   :align: center

   Figure: Example of use of the ``fft_electronic_spin_density`` package. Left: The input is a Gaussian .cube file of charge density. Middle: The output is a .cube file of its Fourier transform (---> visualize in VESTA). Right: The Fourier transform is visualized in 3D, 2D, and 1D cuts. The original DFT density can be manipulated by filtering it out only around selected sites, replacing by a model and even fitting this model to the original density.


Import the ``Density`` class
-------------------------------------------------------------------

.. code-block:: python

    from fft_electronic_spin_density.classes import Density


Load the (spin) density .cube file
-------------------------------------------------------------------

.. code-block:: python

    density = Density(fname_cube_file='../cube_files/Cu2AC4_rho_sz_256.cube')


Visualize the density in 3D
-------------------------------------------------------------------

.. code-block:: python

    density.plot_cube_rho_sz(
                        c_idx_arr=np.arange(0, density.nc),
                        fout_name='rho_sz_exploded.jpg', 
                        alpha=0.05, 
                        figsize=(5.5,5.5), 
                        dpi=400, 
                        zeros_transparent=True, 
                        show_plot=True
                        )

.. 3D density
.. image::
   ./_static/images/rho_sz_exploded.jpg
   :width: 500px
   :align: center


Filter out :math:`\rho_\mathrm{s} (\mathbf{r})` around selected sites
---------------------------------------------------------------------------------

.. code-block:: python

    # selected site indices
    site_idx = [0, 1] # atom 0 - Cu0, atom 1 - Cu1

    # muffin-tin radii around the selected sites where density will be kept
    site_radii = [1.1]*2 # Angstrom

    density.mask_except_sites(leave_sites={
                                'site_centers':density.get_sites_of_atoms(site_idx), 
                                'site_radii':site_radii
                                })
   
and visualize...

.. code-block:: python
   
   density.plot_cube_rho_sz(
                    c_idx_arr=np.arange(0, density.nc, 1), 
                    fout_name='rho_sz_exploded_filtered.jpg', 
                    alpha=0.05, 
                    figsize=(5.5,5.5), 
                    dpi=400, 
                    zeros_transparent=True,
                    show_plot=True,
                    xlims=[0, 6], 
                    ylims=[4,10],
                    zlims=[2,5]
                    )

.. filtered density
.. image::
   ./_static/images/rho_sz_exploded_filtered.jpg
   :width: 400px
   :align: center


Perform FFT, visualize in 2D and 1D
-------------------------------------------------------------------

.. code-block:: python

    density.FFT()


.. code-block:: python

    fft_along_line_data = density.plot_fft_along_line(
                                    i_kz=density.nkc//2, 
                                    cut_along='both', 
                                    kx_ky_fun=None, 
                                    k_dist_lim=12, 
                                    N_points=3001, 
                                    fout_name='cut_1D_both.png', 
                                    cax_saturation=0.5,
                                    )

    kx_arr_along, ky_arr_along, F_abs_sq_interp_along, \
    kx_arr_perp, ky_arr_perp, F_abs_sq_interp_perp = fft_along_line_data

    density.plot_fft_2D(
                i_kz=density.nkc//2, 
                fft_as_log=False, 
                fout_name=f'F_abs_sq-scale_kz_at_idx_{density.nkc//2}_cut_both.png', 
                figsize=(5.5, 4.5),
                dpi=400,
                fixed_z_scale=True,
                cax_saturation=0.5,
                xlims=[-19, 19],
                ylims=[-19, 19],
                zlims=[0, 1.6e6],
                plot_line_cut=True, kx_arr_along=kx_arr_along, ky_arr_along=ky_arr_along,
                kx_arr_perp=kx_arr_perp, ky_arr_perp=ky_arr_perp,
                cut_along='both'
                )

.. FFT 2D plot
.. image::
   ./_static/images/F_abs_sq-scale_kz_at_idx_72_cut_both_fix-scale.png
   :width: 500px
   :align: center

.. FFT 1D cuts
.. image::
   ./_static/images/cut_1D_both.png
   :width: 450px
   :align: center

Write out :math:`\mathcal{F}\{ \rho_\mathrm{s} \}` as a .cube file
--------------------------------------------------------------------------------------------

.. code-block:: python

    density.write_cube_file_fft(fout='fft_rho_sz.cube')

→ visualize the .cube file in VESTA

.. FFT 3D VESTA
.. image::
   ./_static/images/FFT_from_VESTA.png
   :width: 250px
   :align: center


Replace :math:`\rho_\mathrm{s} (\mathbf{r})` by a model
-------------------------------------------------------------------

The model is defined as a d\ :sub:`x2y2`\  orbital centered on Cu sites. 
Possibly, sp orbitals centered on oxygen sites can be added. 
Any parameterized function (e.g., a Gaussian) can be defined as a model.

.. code-block:: python

    site_idx = [0, 1]

    parameters_model = {'type':['dx2y2_neat']*2, 
                        'sigmas':[None]*2, 
                        'centers':density.get_sites_of_atoms(site_idx),
                        'spin_down_orbital_all':[False, True],
                        'fit_params_init_all':{
                            'amplitude':[0.3604531, 0.3604531], 
                            'theta0':   [-1.011437, -1.011437], 
                            'phi0':     [-0.598554, -0.598554], 
                            'Z_eff':    [12.848173, 12.848173],
                            'C':        [0.0000000, 0.0000000]
                            }
                        }

    density.replace_by_model(
                        fit=False, 
                        parameters=parameters_model
                        )


.. code-block:: python

    density.plot_cube_rho_sz(
                        c_idx_arr=np.arange(0, density.nc, 1), 
                        fout_name='rho_sz_exploded_model.jpg', 
                        alpha=0.05, 
                        figsize=(5.5,5.5), 
                        dpi=400, 
                        zeros_transparent=True,
                        show_plot=True,
                        xlims=[0, 6], 
                        ylims=[4,10],
                        zlims=[2,5]
                        )

.. filtered density
.. image::
   ./_static/images/rho_sz_exploded_model.jpg
   :width: 400px
   :align: center


*Fit* the model
-------------------------------------------------------------------

.. code-block:: python

    site_idx = [0, 1]

    parameters_model = {'type':['dx2y2_neat']*2, 
                        'sigmas':[None]*2, 
                        'centers':density.get_sites_of_atoms(site_idx),
                        'spin_down_orbital_all':[False, True],
                        'fit_params_init_all':{
                            'amplitude':[0.3604531, 0.3604531], 
                            'theta0':   [-1.011437, -1.011437], 
                            'phi0':     [-0.598554, -0.598554], 
                            'Z_eff':    [12.848173, 12.848173],
                            'C':        [0.0000000, 0.0000000]
                            }
                        }

    density.replace_by_model(
                        fit=True, 
                        parameters=parameters_model
                        )

| **call 1:**   params [ 0.379 0.361 -1.011 -1.011 -0.599 -0.599 12.848 12.848 0. 0.] **R^2 0.800**
| **call 2:**   params [ 0.361 0.361 -1.011 -1.011 -0.599 -0.599 12.842 12.848 0. 0.] **R^2 0.805**
| **call 3:**   ...


Write out the *modified* :math:`\rho_\mathrm{s} (\mathbf{r})`
-------------------------------------------------------------------

.. code-block:: python

    density.write_cube_file_rho_sz(fout='rho_sz_modified.cube')

... to be visualized in VESTA


Integrate :math:`\rho_\mathrm{s} (\mathbf{r})` over the unit cell
-------------------------------------------------------------------

.. code-block:: python

   rho_tot_unitcell, rho_abs_tot_unitcell = density.integrate_cube_file(verbose=False)

   print(f"""Total charge in the unit cell {rho_tot_unitcell:.4f} e.
   Total absolute charge in the unit cell {rho_abs_tot_unitcell:.4f} e.""")

| Total charge in the unit cell 0.0000 e.
| Total absolute charge in the unit cell 8.1414 e.


| (clearly a magnetically compensated antiferromagnetic spin density in this example unit cell)


Visualize the density as 2D slices
-------------------------------------------------------------------

at selected heights (z-coordinates) along the $c$ lattice vector:

.. code-block:: python

    # z position of atom 0
    site_coordinates = density.get_sites_of_atoms(site_idx=[0])
    atom_0_z_coordinate = site_coordinates[0][2]

    # indices along the c lattice vector where density cuts should be plotted
    c_idx = density.get_c_idx_at_z_coordinates(z_coordinates=[0.0, atom_0_z_coordinate])

    density.plot_cube_rho_sz(
                        c_idx_arr=c_idx, 
                        fout_name='rho_sz_exploded_masked.jpg', 
                        alpha=0.8, 
                        figsize=(6.0, 4.5), 
                        dpi=400, 
                        zeros_transparent=False, 
                        show_plot=True
                        )

.. 2D slices
.. image::
   ./_static/images/plot_2D_example_figure.png
   :width: 500px
   :align: center


    



