try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'AEI image processing and ecological analysis tools',
    'author': 'Christopher B. Anderson',
    'url': 'https://github.com/christobal54/aei-py',
    'download_url': 'https://github.com/christobal54/aei-py.git',
    'author_email': 'cbanders@stanford.edu',
    'version': '0.1',
    'install_requires': ['gdal', 'ogr', 'osr', 'numpy', 'spectral'],
    'packages': ['aei','bin', 'test'],
    'scripts': [],
    'name': 'aei-py'
}

setup(**config)
