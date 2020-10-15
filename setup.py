from setuptools import setup

version = open("earthtools/__version__.py").read().strip('"\n')

setup_args = {
    "name": "earthtools",
    "version": version,
    "url": "https://earth-chris.github.io",
    "license": "MIT",
    "author": "Christopher Anderson",
    "author_email": "cbanders@stanford.edu",
    "description": "Personal tools for image processing and ecological analysis",
    "keywords": ["biogeography", "SDM", "species distribution modeling", "ecologyy", "conservation"],
    "packages": ["earthtools"],
    "include_package_data": True,
    "platforms": "any",
    "scripts": [
        "bin/et-align.py",
        "bin/et-las.py",
        "bin/et-mask.py",
        "bin/et-moment.py",
        "bin/et-palsar.py",
        "bin/et-planet.py",
        "bin/et-sample-raster.py",
        "bin/et-unmix.py",
    ],
}

setup(**setup_args)
