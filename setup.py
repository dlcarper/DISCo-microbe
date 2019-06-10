import re
import os
from setuptools import setup, Extension
from codecs import open
from os import path

version_file = open("disco_microbe/_version.py", "r").read()
version_match = re.match(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file)
if (version_match):
    version = version_match.group(1)
else:
    raise RuntimeError("Unable to find version string in _version.py")

here = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()

setup(name = "disco-microbe",
      setup_requires=['pytest-runner'],
      tests_require=['pytest'],
      install_requires=['biopython'],
      packages = ["disco_microbe"],
      python_requires='~=3.5',
      entry_points = {
          "console_scripts": ['disco = disco_microbe.disco_microbe:main']},
      version = version,
      author="Dana L. Carper, Travis J. Lawrence, Alyssa A. Carrell, Dale A. Pelletier, and David J. Weston",
      author_email="carperdl@ornl.gov",
      description = "Design of identifiable synthetic communities",
      long_description=long_description,
      license='LGPLv3',
      url = "https://github.com/dlcarper/DISCo-Design-of-an-identifiable-synthetic-community",
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
                   'Natural Language :: English',
                   'Operating System :: MacOS :: MacOS X',
                   'Operating System :: POSIX :: Linux',
                   'Programming Language :: Python :: 3.5',
                   'Programming Language :: Python :: 3.6',
                   'Programming Language :: Python :: 3.7',
                   'Programming Language :: Python :: Implementation :: CPython',
                   'Topic :: Scientific/Engineering :: Bio-Informatics'],)
