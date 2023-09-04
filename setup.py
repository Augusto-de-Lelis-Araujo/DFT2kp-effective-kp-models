from setuptools import setup, find_packages

python_version = '>=3.7'

install_requires = [
    'numpy',
    'scipy',
    'sympy==1.10',
    'tabulate',
    'irrep==1.7.1',
    'qsymm==1.3.0',
    'matplotlib',
    'irreptables'
]

# uses README.md as the package long description
with open("README.md") as f:
    long_description = f.read()

from pydft2kp.__version import __version__ as version

setup(
    name="dft2kp",
    description="Calculates kp model and coefficients from DFT ab initio data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.com/dft2kp/dft2kp",
    author="Gerson J. Ferreira",
    author_email="gersonjferreira@ufu.br",
    version=version,
    packages=find_packages('.'),
    install_requires=install_requires,
    python_requires=python_version)
