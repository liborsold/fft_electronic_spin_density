Examples
==========================

An example of use is given in ``./examples/fft_example.ipynb``. It uses input files from ...

.. fft_electronic_spin_density example image
.. image::
   ./_static/images/example_of_use.png
   :width: 800px
   :align: center



Load the .cube file charge or spin density
-------------------------------------------------------------------

.. code-block:: python

    from fft_electronic_spin_density.utils import Density
    density = Density(fname_cube_file='./cube_files/rho_sz.cube')


Integrate density to get the total charge (spin)
-------------------------------------------------------------------

.. code-block:: python

    rho_sz_tot, rho_sz_abs_tot = density.integrate_cube_file()
    

Perform FFT, plot and write out as a .cube file
-------------------------------------------------------------------

.. code-block:: python

    density.FFT()
    density.plot_fft_2D(i_kz=0)
    density.write_cube_file_fft(fout='fft_rho_sz.cube')



    



