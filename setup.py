#!/usr/bin/env python3
"""
ImmunOmics — Multi-omics integration of immune cell states.
"""
from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="immunomics",
    version="0.1.0",
    description="Multi-omics integration of immune cell states from matched scRNA-seq + scATAC-seq",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Patrick Grady",
    url="https://github.com/pgrady3/immunomics",
    packages=find_packages(),
    include_package_data=True,
    python_requires=">=3.10",
    install_requires=[
        "click>=8.0",
        "scanpy>=1.9",
        "anndata>=0.10",
        "muon>=0.1.5",
        "numpy>=1.24",
        "pandas>=2.0",
        "scipy>=1.10",
        "scikit-learn>=1.3",
        "scvi-tools>=1.0",
        "matplotlib>=3.7",
        "seaborn>=0.12",
        "pyyaml",
        "tqdm",
        "rich",
    ],
    extras_require={
        "r": ["rpy2>=3.5"],
        "dev": ["pytest>=7.0", "pytest-cov", "black", "isort"],
    },
    entry_points={
        "console_scripts": [
            "immunomics=immunomics.cli:main",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
)
