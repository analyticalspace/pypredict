#!/usr/bin/env python

from distutils.core import setup, Extension
setup(
    name='pypredict',
    version='2.0',
    author="Jesse Trutna",
    author_email="jesse@spire.com",
    url="https://github.com/analyticalspace/pypredict",
    py_modules=['predict'],
    ext_modules=[Extension('cpredict', ['predict.c', 'pypredict.c'])]
    )
