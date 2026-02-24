#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ImmunOmics v0.1.0

Package installation and dependency configuration.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
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
        "numpy>=1.24",
        "pandas>=2.0",
        "scipy>=1.10",
        "scikit-learn>=1.3",
        "matplotlib>=3.7",
        "seaborn>=0.12",
        "pyyaml",
        "tqdm",
        "rich",
        "statsmodels>=0.14",
        "scikit-misc>=0.3",
    ],
    extras_require={
        "ml": ["scvi-tools>=1.0"],
        "r": ["rpy2>=3.5"],
        "multimodal": ["muon>=0.1.5", "pybiomart>=0.2"],
        "dev": ["pytest>=7.0", "pytest-cov", "ruff>=0.4"],
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

# ImmunOmics v0.1.0
# Any usage is subject to this software's license.
