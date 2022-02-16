from setuptools import setup

# This value will be reassigned when version.py is parsed
__version__ = "x.x.x"
exec(open("flutile/version.py", "r").read())

# Read the requirements from the requirements.txt file
with open("requirements.txt", "r") as fh:
    requirements = [r.strip() for r in fh.readlines()]


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="flutile",
    version=__version__,
    description="sequence analysis tools for flu research",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/flu-crew/flutile",
    author="Zebulun Arendsee",
    author_email="zebulun.arendsee@usda.gov",
    packages=["flutile"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points={"console_scripts": ["flutile=flutile.ui:main"]},
    install_requires=requirements,
    py_modules=["flutile"],
    zip_safe=False,
    package_data={"flutile": ["py.typed"]},
    include_package_data=True,
)
