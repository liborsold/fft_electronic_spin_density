# fft_electronic_spin_density
**<i>Perform FFT on Gaussian .cube files of charge or spin density, primarily to obtain the (magnetic) form factor for neutron scattering.</i>**

```
pip install fft_electronic_spin_density
```

See the [documentation website](https://liborsold.github.io/fft_electronic_spin_density/build/html/index.html).

## Usage

```python
from fft_electronic_spin_density.utils import Density

density = Density(fname_cube_file='./cube_files/rho_sz.cube')

rho_sz_tot, rho_sz_abs_tot = density.integrate_cube_file()

density.FFT()
density.plot_fft_2D(i_kz=0)
density.write_cube_file_fft(fout='fft_rho_sz.cube')
```

<center><img src="https://liborsold.github.io/fft_electronic_spin_density/build/html/_images/example_of_use.png" alt="fft_spin_density_example" width="700" /></center>




