import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="afino",
    version="0.5",
    author="Andrew Inglis",
    author_email="andrew.inglis@nasa.gov",
    description="A package for finding oscillations in timeseries",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/afino_release_version/",
    packages=setuptools.find_packages(),
)
