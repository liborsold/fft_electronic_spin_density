Examples
==========================

An example of use is given in ``./examples/fft_example.ipynb``. It uses input files from the `cube_files folder <https://github.com/liborsold/fft_electronic_spin_density/tree/master/cube_files>`_.

.. fft_electronic_spin_density example image
.. image::
   ./_static/images/example_of_use.png
   :width: 800px
   :align: center



Load the .cube file charge or spin density
-------------------------------------------------------------------

.. code-block:: python

    from fft_electronic_spin_density.utils import Density
    density = Density(fname_cube_file='../cube_files/Cu2AC4_rho_sz_256.cube')


Visualize the density
-------------------------------------------------------------------

.. code-block:: python

    # plot 2D

    # ---> show plot

    # plot 3D

    # ---> show plot


Replace by a model
-------------------------------------------------------------------

.. code-block:: python

    # replace by model


*or even* Fit model to the original density 
-------------------------------------------------------------------

.. code-block:: python

    # fit model


Integrate density to get the total charge (spin)
-------------------------------------------------------------------

.. code-block:: python

    rho_sz_tot, rho_sz_abs_tot = density.integrate_cube_file()

    # show the output
    

Perform FFT, plot and write out as a .cube file
-------------------------------------------------------------------

.. code-block:: python

    density.FFT()

    density.plot_fft_2D(i_kz=0)

    # ---> show plot

    # plot along cuts in the 2D map
    
    # ---> show plot

    # PLOT 3D
    # ---> show plot

    # WRITE OUT

    density.write_cube_file_fft(fout='fft_rho_sz.cube')

    # ----> show how it's visualized in VESTA



    



