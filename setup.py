#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

requirements = [
    'h5py >= 2.7.0', 'numpy >=1.8.2', 'SciPy >= 0.13.3',
    'scikit-learn >=0.19', 'pandas >= 0.18', 'POT>=0.4.0'
]

extras_require = {
    'GRN': ["'gslrandom>=0.1.4'"],
}

setup_requirements = [
    # put setup requirements (distutils extensions, etc.) here
]

test_requirements = [
    'unittest'
]

setuptools.setup(
    name='wot',
    version='0.1.0',
    description="Optimal transport for time-course single cell data",
    author="Geoffrey Schiebinger, Jian Shu, Marcin Tabaka, Brian Cleary",
    author_email='wot@broadinstitute.org',
    url='https://github.com/broadinstitute/wot',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(include=['wot']),
    include_package_data=True,
    install_requires=requirements,
    license="BSD license",
    zip_safe=False,
    keywords='wot',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
    python_requires='~=3.0',
    entry_points={
        'console_scripts': [
            'wot=wot.__main__:main',
        ]
    }
)
