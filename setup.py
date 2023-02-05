from setuptools import setup, find_packages

setup(
    name="SeqPanther",
    version="0.0.1",
    install_requires=[
        'biopython>=1.80', 'click>=7.1.2', 'numpy>=1.22.1', 'pandas>=1.5.2',
        'pyfaidx>=0.6.3.1', 'pysam>=0.18.0', 'matplotlib>=3.6.2'
    ],
    packages=find_packages(include=["seqPanther", "seqPanther.*"]),
    entry_points={
        'console_scripts': ['seqpanther=seqPanther.seqPanther:run'],
    },
    url="https://github.com/krisp-kwazulu-natal/seqPatcher",
    license="GPLv3",
    author="Anmol Kiran; San James Emmanuel",
    author_email="anmol.kiran@gmail.com;sanemmanueljames@gmail.com",
    description="A set of sequence manipulation tools",
)
