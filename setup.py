#!/usr/bin/env python

from distutils.core import setup

setup (name = "DensityFit",
       version = "0.1.4",
       description = "Flexible fit of protein structures into low-resolution electronic density maps",
       long_description =
"""This program takes as input a low-resolution electronic density map
(e.g. from electron microscopy) for a protein and a structure of the
same protein in some other conformation. It then deforms the conformation
to fit into the electronic density while preserving the protein structure
through the use of low-energy normal modes.
""",
       author = "K. Hinsen",
       author_email = "hinsencnrs-orleans.fr",
       url = "http://dirac.cnrs-orleans.fr/",
       license = "CeCILL",
       packages = ['DensityFit'],
       scripts = ['density_fit'],
       )

