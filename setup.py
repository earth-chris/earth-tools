try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'personal image processing and ecological analysis tools',
    'author': 'Christopher B. Anderson',
    'url': 'https://earth-chris.github.io',
    'download_url': 'https://github.com/christobal54/earth-tools.git',
    'author_email': 'cbanders@stanford.edu',
    'version': '0.1',
    'install_requires': ['gdal', 'ogr', 'osr', 'numpy', 'spectral', 'psutil'],
    'packages': ['earthtools','bin'],
    'name': 'earthtools'
}

setup(**config)
