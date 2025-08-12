from setuptools import setup, find_packages
import os

setup(
    name="tbnexplorer2",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
    ],
    entry_points={
        "console_scripts": [
            "tbnexplorer2=tbnexplorer2.cli:main",
            "tbnexplorer2-filter=tbnexplorer2.filter_cli:main",
        ],
    },
    python_requires=">=3.8",
    author="TBN Explorer Team",
    description="A Python library and CLI tool for thermodynamics of binding networks (TBN) analysis",
    long_description=open("README.md").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
)