from setuptools import setup, find_packages

setup(
    name="genal-python",
    version="0.0",
    description="Cleaning and processing of GWAS data, polygenic risk scoring, mendelian randomization including a parallel implementation of MR-PRESSO",
    author="Cyprien Rivier",
    author_email="riviercyprien@gmail.com",
    url="github...",
    license="GPL",
    packages=find_packages(),
    install_requires=[
        'numpy==1.21.0',
        'pandas>=1.3',
        'requests<3.0',
        'scipy>=1.6,<=1.7',
    ],
    # ... other metadata ...
)