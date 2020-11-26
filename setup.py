#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages
import sys
from pathlib import Path

# with open('README.rst') as readme_file:
#     readme = readme_file.read()

def readme():
    if sys.version_info > (3, ):
        with open(Path(__file__).parent.resolve() / 'README.md', encoding='utf-8') as md:
            return md.read()
    else:
        with open('README.md') as md:
            return md.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['numpy >= 1.15.4', 'scipy >= 1.1.0', 'pandas >= 0.24', 'statsmodels',
                'KDEpy', 'pybedtools', 'pyBigWig', 'h5py', 'matplotlib']

setup_requirements = ['pytest-runner']

test_requirements = ['pytest>=3']

setup(
    author="Nanxiang Zhao (Samuel)",
    author_email='samzhao@umich.edu',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Improving the feature density based peak caller with dynamic statistics.",
    long_description=readme(),
    long_description_content_type="text/markdown",
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
    version='2.0.1',
    zip_safe=False,
)
