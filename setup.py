# coding: utf-8
"""
Parts of this are borrowed from Kenneth Reitz's
setup.py project: https://github.com/kennethreitz/setup.py
"""
from __future__ import with_statement
import os, os.path as osp
import sys
from shutil import rmtree

from setuptools import find_packages, setup, Command

# Package meta-data.
NAME = 'clark-biom'
DESCRIPTION = "Create BIOM-format tables from CLARK abundance tables."
URL = 'https://github.com/smdabdoub/clark-biom'
EMAIL = 'dabdoub.2@osu.edu'
AUTHOR = 'Shareef M. Dabdoub'
REQUIRES_PYTHON = '>=2.7.14'

# What packages are required for this module to be executed?
REQUIRED = ["biom-format >= 2.1.5"]

here = osp.abspath(osp.dirname(__file__))

def get_version():
    with open(osp.join(here,'clark_biom.py')) as f:
        for line in f:
            if line.startswith('__version__'):
                return eval(line.split('=')[-1])
 
 
def get_long_description():
    descr = []
    for fname in (osp.join(here,'README.rst'), 
                  osp.join(here,'CHANGELOG.rst')):
        with open(fname) as f:
            descr.append(f.read())
    return '\n\n'.join(descr)
 

class UploadCommand(Command):
    """Support setup.py upload."""

    description = 'Build and publish the package.'
    user_options = []
    
    @staticmethod
    def status(s):
        """Prints things in bold."""
        print('\033[1m{0}\033[0m'.format(s))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status('Removing previous builds…')
            rmtree(os.path.join(here, 'dist'))
        except OSError:
            pass

        self.status('Building Source and Wheel (universal) distribution…')
        os.system('{0} setup.py sdist bdist_wheel --universal'.format(sys.executable))

        self.status('Uploading the package to PyPi via Twine…')
        os.system('twine upload dist/*')

        sys.exit()
 
setup(
    name=NAME,
    version=get_version(),
    description=DESCRIPTION,
    long_description=get_long_description(),
    keywords='CLARK, BIOM, metagenomics, bioinformatics, taxonomy, taxonomic-classification',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    license='MIT',
    py_modules=['clark_biom'],
    namespace_packages=[],
    include_package_data=True,
    zip_safe=False,
    install_requires=REQUIRED,
    entry_points={
        'console_scripts': [
            'clark-biom = clark_biom:main',
        ],
    },
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    # $ setup.py publish support. 
    cmdclass={
        'upload': UploadCommand,
    },
)