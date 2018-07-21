from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize

from codecs import open
from os import path
import re


here = path.abspath(path.dirname(__file__))
package_name = 'Schottky'
version_file = path.join(here, package_name, '_version.py')
with open(version_file, 'rt') as f:
    version_file_line = f.read()
version_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(version_re, version_file_line, re.M)
if mo:
    version_string = mo.group(1)
else:
    raise RuntimeError('Unable to find version string in %s.' % (version_file,))

readme_file = path.join(here, 'README.md')
with open(readme_file, encoding='utf-8') as f:
    long_description = f.read()

extensions = [
    Extension(
        'Schottky.Constants',
        ['Schottky/Constants.pyx'],
        depends=['Schottky/Constants.pxd'],
    ),
    Extension(
        'Schottky.Metal',
        ['Schottky/Metal/__init__.pyx'],
        depends=['Schottky/Metal/__init__.pxd'],
    ),
]

setup(
    name=package_name,
    version=version_string,

    description='Schottky diode simulator',
    long_description=long_description,
    long_description_content_type='text/markdown',

    url='https://github.com/bond-anton/Schottky',

    author='Anton Bondarenko',
    author_email='bond.anton@gmail.com',

    license='Apache Software License',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],

    keywords='Semiconductors, Schottky diode, electronics',

    packages=find_packages(exclude=['demo', 'tests', 'docs', 'contrib', 'venv']),
    ext_modules=cythonize(extensions),
    package_data={'Schottky.Constants': ['Constants.pxd'],
                  'Schottky.Metal': ['Schottky/Metal/__init__.pxd']},
    install_requires=['numpy', 'Cython', 'scipy', 'matplotlib', 'BDMesh', 'BDPoisson1D'],
    test_suite='nose.collector',
    tests_require=['nose']
)
