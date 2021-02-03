# Copyright (C) 2020 Samuel Baker

DESCRIPTION = "High level API for genetic pipelines"
LONG_DESCRIPTION = """
# pyGenicPipeline

High level API for genetic pipelines

## About The Project

pyGenicPipeline is simply a wrapper around pre-existing projects that have been altered to make them (hopefully) easier
to use and generalised. 

All the source code can be found at the [pyGenicPipeline git repository](https://github.com/sbaker-dev/pyGenicPipeline)
with additional information found on the [Docs Page](https://sbaker-dev.github.io/pyGenicPipeline/.)

"""
LONG_DESCRIPTION_CONTENT_TYPE = "text/markdown"

DISTNAME = 'pyGenicPipeline'
MAINTAINER = 'Samuel Baker'
MAINTAINER_EMAIL = 'samuelbaker.researcher@gmail.com'
LICENSE = 'MIT'
DOWNLOAD_URL = "https://github.com/sbaker-dev/pyGenicPipeline"
VERSION = "0.09.6"
PYTHON_REQUIRES = ">=3.7"

INSTALL_REQUIRES = [

    'pysnptools',
    'csvObject',
    'numpy==1.19.3',
    'colorama',
    'zstd',
    'scipy',
    'PyYAML',
    'miscSupports',
    'pyGenicParser']

CLASSIFIERS = [
    'Programming Language :: Python :: 3.7',
    'License :: OSI Approved :: MIT License',
]

if __name__ == "__main__":

    from setuptools import setup, find_packages

    import sys

    if sys.version_info[:2] < (3, 7):
        raise RuntimeError("pyGenicPipeline requires python >= 3.7.")

    setup(
        name=DISTNAME,
        author=MAINTAINER,
        author_email=MAINTAINER_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        long_description_content_type=LONG_DESCRIPTION_CONTENT_TYPE,
        license=LICENSE,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        python_requires=PYTHON_REQUIRES,
        install_requires=INSTALL_REQUIRES,
        include_package_data=True,
        packages=find_packages(),
        classifiers=CLASSIFIERS
    )
