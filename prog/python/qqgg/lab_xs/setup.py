from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize(
        "lab_xs.pyx", annotate=True, compiler_directives=dict(cdivision=True)
    )
)
