# libmadam
Library version of the Madam mapmaking code

The Madam mapmaking algorithm is discussed in the following research papers:

E. Keihänen, H. Kurki-Suonio and T. Poutanen:  
*Madam - A Map-making method for CMB experiments*,  
MNRAS **360** (2005) 390, [astro-ph/0412517](https://arxiv.org/abs/astro-ph/0412517)

E. Keihänen, R. Keskitalo, H. Kurki-Suonio, T. Poutanen, A.-S. Sirviö:  
*Making CMB temperature and polarization maps with Madam*,  
A&A **510** (2010) A57, [arXiv:0907.0367](https://arxiv.org/abs/0907.0367)

## Installation

To compile, test and install the Fortran library with C-bindings:
```bash
./autogen.sh
./configure
make && make check && make install
```

To install and test the Python wrapper
```bash
cd python
python setup.py install
python setup.py test
```
