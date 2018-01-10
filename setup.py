#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

requirements = [
    'dask >= 0.16.0', 'h5py >= 2.7.0', 'numpy >=1.8.2', 'SciPy >= 0.13.3',
    'scikit-learn >=0.19',
    'pandas >= 0.18'
]

setup_requirements = [
    # put setup requirements (distutils extensions, etc.) here
]

test_requirements = [
    'unittest'
]

setup(
    name='wot',
    version='0.1.0',
    description="Optimal transport for time-course single cell data",
    long_description="Uses time-course data to infer how the probability "
                     "distribution of cells "
                     "evolves over time, by using the mathematical approach "
                     "of Optimal Transport (OT)",
    author="Geoffrey Schiebinger, Jian Shu, Marcin Tabaka, Brian Cleary",
    author_email='wot@broadinstitute.org',
    url='https://github.com/broadinstitute/wot',
    packages=find_packages(include=['wot']),
    entry_points={
    },
    include_package_data=True,
    install_requires=requirements,
    license="BSD license",
    zip_safe=False,
    keywords='wot',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
)
