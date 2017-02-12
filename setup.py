#   md - Molecular Dynamics Applied to ionic solids.
#   Copyright (C) 2017 Nils Harmening, Marco Manni,
#   Darian Steven Viezzer, Steffi, Hendrik
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.


from setuptools import setup
from setuptools import Extension
import versioneer, sys, os

class lazy_cythonize(list):
    """evaluates extension list lazyly.
    pattern taken from http://tinyurl.com/qb8478q"""
    def __init__(self, callback):
        self._list, self.callback = None, callback
    def c_list(self):
        if self._list is None: self._list = self.callback()
        return self._list
    def __iter__(self):
        for e in self.c_list(): yield e
    def __getitem__(self, ii): return self.c_list()[ii]
    def __len__(self): return len(self.c_list())
            
def extensions():
    from numpy import get_include
    from Cython.Build import cythonize
    ext_fast_neighbourlist = Extension(
        "md._ext.fast_neighbourlist",
        sources=["md/_ext/fast_neighbourlist.pyx", "md/_ext/src_fast_neighbourlist.c"],
        include_dirs=[get_include()],
        extra_compile_args=["-O3", "-std=c99"])
    exts = [ext_fast_neighbourlist]
    return cythonize(exts)


def get_cmdclass():
    versioneer_cmds = versioneer.get_cmdclass()
    class sdist(versioneer_cmds['sdist']):
        """ensure cython files are compiled to c, when distributing"""
        def run(self):
            # only run if .git is present
            if not os.path.exists('.git'):
                print("Not on git, can not create source distribution")
                return
            try:
                from Cython.Build import cythonize
                print("cythonizing sources")
                cythonize(extensions())
            except ImportError:
                warnings.warn('sdist cythonize failed')
                return versioneer_cmds['sdist'].run(self)
    versioneer_cmds['sdist'] = sdist
    from setuptools.command.test import test as TestCommand
    class PyTest(TestCommand):
        user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]
        def initialize_options(self):
            TestCommand.initialize_options(self)
            self.pytest_args = ['md']
        def run_tests(self):
            # import here, cause outside the eggs aren't loaded
            import pytest
            errno = pytest.main(self.pytest_args)
            sys.exit(errno)
    versioneer_cmds['test'] = PyTest
    return versioneer_cmds


setup(
    cmdclass=get_cmdclass(),
    ext_modules=lazy_cythonize(extensions),
    name='md',
    version=versioneer.get_version(),
    description="Simulation of molecular dynamics applied to ionic solids",
    classifiers=[
        'Development Status :: a - Planning',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        #'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics'],
    keywords=[],
    url='https://github.com/mannimaster/compsci-proj-md',
    author='Nils Harmening, Marco Manni, Darian Steven Viezzer, Steffi, Hendrik',
    author_email='nils.harmening@fu-berlin.de',
    license='GPLv3+',
    packages=['md', 'md._ext'],
    setup_requires=['numpy>=1.7.0', 'cython>=0.22', 'setuptools>=0.6'],
    install_requires=['numpy>=1.7.0'],
    tests_require=['pytest', 'cython>=0.22'])
