from setuptools import setup
import os

about = {}
# Read version number from __version__.py (see PEP 396)
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'abnumber', '__version__.py'), encoding='utf-8') as f:
    exec(f.read(), about)

# Read contents of readme file into string
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='abnumber',
    version=about['__version__'],
    description='AbNumber - Antibody numbering using ANARCI',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='David Prihoda',
    packages=['abnumber'],
    author_email='david.prihoda@gmail.com',
    license='MIT',
    python_requires=">=3.6",
    keywords='antibody numbering,immunoglobulin,ab,ANARCI,chain,imgt,chothia',
    include_package_data=True,
    url='https://github.com/prihoda/abnumber',
    install_requires=[
        "biopython",
        "pandas",
        "anarcii",
        # ANARCI dependency optional - can be installed separately from bioconda or github
    ],
)
