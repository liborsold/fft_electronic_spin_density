Theory
===============================

High-temperature superconductivity in cuprates can be studied by inelastic neutron scattering (INS) experiments.


Neutron scattering
----------------------------

The INS spectra of a copper dimer with non-overlapping magnetic orbitals

.. math::
    \begin{equation}
        S_{zz}=\frac{1}{2} |F(\vec{k})|^2 \left(1-\cos \left(\vec{k} \cdot \vec{r}_{a b}\right)\right) \, \delta\left(\hbar \omega+E_S-E_{T_0}\right).
    \end{equation}

includes the *magnetic form factor*

.. math::
    \begin{equation}
        F(\mathbf{\kappa}) \equiv e^{i \mathbf{\kappa} \cdot \mathbf{r}} \rho_\mathrm{s} (\mathbf{r}) d\mathbf{r} = \int e^{i \mathbf{\kappa} \cdot \mathbf{r}} \rho_\mathrm{s} (\mathbf{r})
    \end{equation}

which is the **Fourier transform of the spin density** :math:`\rho_\mathrm{s} (\mathbf{r})`.

The spin density is obtained from a density functional theory (DFT) calculation and usually output as a Gaussian .cube file. 

**The (spin) density from a .cube file can then be loaded, filtered out, Fourier transformed, and visualized by the present** ``fft_electronic_spin_density`` **package.**

| 

Keeping the spin density only around selected sites is important to get rid of spurious spectra and evaluate the effect of the oxygen ligands.

To assess the influence of a possible overlap (obtaining the spectral function :math:`E_\perp` under the Heitler-London approximation), the ``fft_electronic_spin_density`` package also allows to replace the DFT-calculated density by a model atomic orbitals, which can be fitted to the original density.
Such model is very useful to, e.g., study the :math:`E_\perp` dependence on the Cu-Cu separation :math:`|r_{ab}|`.


Fourier transform
==================================================

Note that the dominant feature in the INS spectra, which is the stratification 
due to the :math:`\left(1-\cos \left(\vec{k} \cdot \vec{r}_{a b}\right)\right)` term, 
arises in general for any Fourier transform of a repeated displaced object. 


.. FFT general
.. image::
   ./_static/images/FFT_general.png
   :width: 500px
   :align: center


- resolution and system size are reciprocal quantities in real and reciprocal space
- show how the zero-padding is used to increase the resolution of the FFT