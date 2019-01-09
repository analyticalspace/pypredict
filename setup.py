#!/usr/bin/env python

from distutils.core import setup, Extension
setup(
    name='pypredict',
    version='2.1',
    author="Jesse Trutna",
    author_email="jesse@spire.com",
    maintainer="Ben Gaudiosi",
    maintainer_email="ben.gaudiosi@analyticalspace.com",
    url="https://github.com/analyticalspace/pypredict/",
    py_modules=['predict'],
    ext_modules=[Extension('cpredict', ['predict.c', 'pypredict.c'])]
    )
