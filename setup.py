#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

# with open('README.rst') as readme_file:
#     readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['numpy', 'scipy >= 0.13', 'pandas >= 0.24', 'statsmodels',
                'KDEpy', 'pybedtools', 'pyBigWig', 'h5py', 'matplotlib']

setup_requirements = ['pytest-runner']

test_requirements = ['pytest>=3']

setup(
    author="Nanxiang Zhao (Samuel)",
    author_email='samzhao@umich.edu',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="A feature density estimator for high-throughput sequence tags.",
    scripts=['bin/fseq2'],
    install_requires=requirements,
    license="GNU General Public License v3",
    include_package_data=True,
    keywords='fseq2',
    name='fseq2',
    packages=find_packages(include=['fseq2', 'fseq2.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/Boyle-Lab/F-Seq2',
    version='2.0',
    zip_safe=False,
)
