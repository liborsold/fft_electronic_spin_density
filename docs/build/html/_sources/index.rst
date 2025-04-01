.. fft_electronic_spin_density documentation master file, created by
   sphinx-quickstart on Thu Mar 27 15:43:55 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ``fft_electronic_spin_density``!
=================================================================

The ``fft_electronic_spin_density`` Python package performs fast Fourier transform (FFT) on Gaussian .cube files of charge or spin density, primarily to obtain the (magnetic) form factor for neutron scattering.


.. fft_electronic_spin_density example image
.. image::
   ./_static/images/example_of_use.png
   :width: 500px
   :align: center

\ 

See the `project on GitHub <https://github.com/liborsold/fft_electronic_spin_density/tree/master>`_.

\ 

Quick start
================

.. code-block:: python
   
   pip install fft_electronic_spin_density

then 

.. code-block:: python
   
   git clone https://github.com/liborsold/fft_electronic_spin_density.git
   
and execute the ``./examples/example_of_use.ipynb`` Jupyter notebook in the `examples folder <https://github.com/liborsold/fft_electronic_spin_density/tree/master/examples>`_.

See the `Examples page <https://liborsold.github.io/fft_electronic_spin_density/build/html/examples.html>`_.

Citation
===============
If you find this package useful, please cite *L. Spitz, L. Vojáček, et al., under preparation.*


Navigation
===============

.. toctree::
   :maxdepth: 1

   examples
   density_class
   api
   theory

Indices
==================

* :ref:`genindex`
* :ref:`modindex`
