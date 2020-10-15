# earthtools

Contains functions and classes for image processing and ecological analysis - it's a personal platform, and fairly hard to parse—I effectively learned how to write in python as I developed this package.

## Install

The easiest way is to install via [conda][home-conda]. Make sure you have `conda` installed (either via minconda (recommended) or anaconda).

### via conda

```bash
git clone https://github.com/earth-chris/earth-tools.git
cd earth-tools/
conda env update
```

Once you've created the environment, activate it and install `earthtools`.

```bash
conda activate earthtools
pip install .
```

Then you should have a conda environment you can actiave with `conda activate earthtools`. You can then e.g. run the executable `et-align.py -h`, or `import earthtools as et` in python from this environment.

If you're interested in using the ccb default `ipython` profile, you can set an environment variable to do this for you. From the base `earthtools` directory, run the following:

```bash
conda activate earthtools
conda env config vars set IPYTHONDIR=$PWD/ipython
```

You'll have to run `conda deactivate` then `conda activate earthtools` for the changes to take effect. After that you'll be able to run `ipython` with our default settings.

### Required binary packages
**For all users**
- GDAL/OGR - from [source](http://download.osgeo.org/gdal/) or [QGIS](http://qgis.org/en/site/forusers/download.html)
- LAStools - from [Rapidlasso](https://rapidlasso.com/lastools/) - Don't forget the [license!](http://www.cs.unc.edu/~isenburg/lastools/LICENSE.txt)

**For linux/unix users**
- Wine - from [WineHQ](https://www.winehq.org/download)

### Environment variables
Check `earthtools/params.py` for the environment variables to set. If all of the executables for the above binary packages are accessible from the command line, no additional changes are necessary.
- `$LTBASE` - path to the LAStools executables (e.g. ~/src/lastools/bin)
- `$GDALBASE` - path to GDAL executables
- `$OUTPUT_DIR` - the default output directory for scripts
- `$SCRATCH_DIR` - the output directory for temp files

## Contact info:

Christopher Anderson
- [E-mail](mailto:cbanders@stanford.edu)
- [Twitter](http://twitter.com/earth_chris)
- [Google Scholar](https://scholar.google.com/citations?user=LoGxS40AAAAJ&hl=encba@anderson-ubuntu:)
- [Library Thing](http://www.librarything.com/catalog/anderzen)
