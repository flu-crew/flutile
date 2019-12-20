from setuptools import setup

from flutile.version import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="flutile",
    version=__version__,
    description="sequence analysis tools for flu research",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/arendsee/flutile",
    author="Zebulun Arendsee",
    author_email="zbwrnz@gmail.com",
    packages=["flutile"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={"console_scripts": ["flutile=flutile.main:main"]},
    py_modules=["flutile"],
    zip_safe=False,
)
