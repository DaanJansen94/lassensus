from setuptools import setup, find_packages

setup(
    name="lassensus",
    version="0.0.5",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.79",
        "pandas>=1.3.0",
        "requests>=2.26.0",
        "numpy<2.0",
    ],
    entry_points={
        "console_scripts": [
            "lassensus=lassensus.lassensus:main",
        ],
    },
    python_requires=">=3.6",
    author="Daan Jansen",
    author_email="your.email@example.com",
    description="A tool for Lassa virus consensus sequence generation",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/DaanJansen94/lassensus",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL License",
        "Operating System :: OS Independent",
    ],
) 
