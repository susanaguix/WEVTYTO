from setuptools import setup, find_packages

setup(
    name='WEVTYTO',
    version='0.1',
    packages=find_packages(),
    scripts=['WEVTYTO/ev_typing_nix.py', 'WEVTYTO/ev_typing_shaw.py'],
    package_data={'': ['WEVTYTO/ev_reference_sequences.fasta']},
)