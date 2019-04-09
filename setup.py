from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.extension import Extension
from Cython.Build import cythonize

from codecs import open
import sys
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
    Extension(
        'Schottky.Trap',
        ['Schottky/Trap/__init__.pyx'],
        depends=['Schottky/Trap/__init__.pxd'],
    ),
    Extension(
        'Schottky.Dopant',
        ['Schottky/Dopant/__init__.pyx'],
        depends=['Schottky/Dopant/__init__.pxd'],
    ),
    Extension(
        'Schottky.Semiconductor',
        ['Schottky/Semiconductor/__init__.pyx'],
        depends=['Schottky/Semiconductor/__init__.pxd'],
    ),
    Extension(
        'Schottky.Potential.ExternalField',
        ['Schottky/Potential/ExternalField.pyx'],
        depends=['Schottky/Potential/ExternalField.pxd'],
    ),
    Extension(
        'Schottky.Potential.PointLike',
        ['Schottky/Potential/PointLike.pyx'],
        depends=['Schottky/Potential/PointLike.pxd'],
    ),
    Extension(
        'Schottky.Potential.TrapPotential',
        ['Schottky/Potential/TrapPotential.pyx'],
        depends=['Schottky/Potential/TrapPotential.pxd'],
    ),
    Extension(
        'Schottky.Helpers.array',
        ['Schottky/Helpers/array.pyx'],
        depends=['Schottky/Helpers/array.pxd'],
    ),
    Extension(
        'Schottky.Helpers.Cache',
        ['Schottky/Helpers/Cache.pyx'],
        depends=['Schottky/Helpers/Cache.pxd'],
    ),
]

copt = {'msvc': ['/openmp', '/Ox', '/fp:fast', '/favor:INTEL64', '/Og'],
        'mingw32': ['-fopenmp', '-O3', '-ffast-math', '-march=native'],
        'unix': ['-fopenmp', '-O3', '-ffast-math', '-march=native']}
lopt = {'mingw32': ['-fopenmp'],
        'unix': ['-fopenmp']}


# check whether compiler supports a flag
def has_flag(compiler, flagname):
    import tempfile
    from distutils.errors import CompileError
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except CompileError:
            return False
    return True


# filter flags, returns list of accepted flags
def flag_filter(compiler, flags):
    result = []
    for flag in flags:
        if has_flag(compiler, flag):
            result.append(flag)
    return result


class CustomBuildExt(build_ext):
    def build_extensions(self):
        c = self.compiler.compiler_type
        print('Compiler:', c)
        opts = flag_filter(self.compiler, copt.get(c, []))
        lopts = flag_filter(self.compiler, lopt.get(c, []))
        for e in self.extensions:
            e.extra_compile_args = opts
            e.extra_link_args = lopts
        build_ext.build_extensions(self)

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
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],

    keywords='Semiconductors, Schottky diode, electronics',

    packages=find_packages(exclude=['demo', 'tests', 'docs', 'contrib', 'venv']),
    ext_modules=cythonize(extensions, compiler_directives={'language_level': sys.version_info[0]}),
    package_data={'Schottky': ['*.pxd'],
                  'Schottky/Helpers': ['*.pxd'],
                  'Schottky/Metal': ['*.pxd'],
                  'Schottky/Trap': ['*.pxd'],
                  'Schottky/Dopant': ['*.pxd'],
                  'Schottky/Semiconductor': ['*.pxd'],
                  'Schottky/Potential': ['*.pxd'],
                  },
    install_requires=['numpy', 'Cython', 'scipy', 'matplotlib', 'BDSpace', 'BDMesh', 'BDPoisson1D'],
    test_suite='nose.collector',
    tests_require=['nose'],
    cmdclass={'build_ext': CustomBuildExt},
    zip_safe=False
)
