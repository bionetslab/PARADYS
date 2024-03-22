
from setuptools import setup, find_packages

setup(
    name='paradys',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'graph-tools',
        'scipy'
    ],
    author='Sandra Goizueta',
    description='Patient-specific ranking of genes driving dysregulations',
)
