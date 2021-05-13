#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'rich',
    'typer',
    'jinja2',
    'pandas',
    'biopython',
    'requests'
]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    author="Peter Kruczkiewicz",
    author_email='peter.kruczkiewicz@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description="Generate a standalone HTML file with an interactive phylogenetic tree using PhyloCanvas",
    entry_points={
        'console_scripts': [
            'shiptv=shiptv.cli:app',
        ],
    },
    install_requires=requirements,
    license="Apache Software License 2.0",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='shiptv',
    name='shiptv',
    packages=find_packages(include=['shiptv']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/peterk87/shiptv',
    version='0.4.1',
    zip_safe=False,
)
