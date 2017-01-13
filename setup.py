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
import versioneer

setup(
    cmdclass=versioneer.get_cmdclass(),
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
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics'],
    keywords=[],
    url='https://github.com/cwehmeyer/pysor',
    author='Nils Harmening, Marco Manni, Darian Steven Viezzer, Steffi, Hendrik',
    author_email='nils.harmening@fu-berlin.de',
    license='GPLv3+',
    packages=['md'],
    #install_requires=['numpy>=1.7.0', 'cython>=0.22'],
    )
