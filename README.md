# aei
Contains functions and classes for image processing and ecological analysis - a personal platform. The acronym stands for anderson ecological imaging.

## Installation
No pip/easy_install/etc. method is currently available. Update your `PYTHONPATH`.

### Required binary packages
**For all users**
- GDAL/OGR - from [source](http://download.osgeo.org/gdal/) or [QGIS](http://qgis.org/en/site/forusers/download.html)
- LAStools - from [Rapidlasso](https://rapidlasso.com/lastools/) - Don't forget the [license!](http://www.cs.unc.edu/~isenburg/lastools/LICENSE.txt)

**For linux/unix users**
- Wine - from [WineHQ](https://www.winehq.org/download)

### Environment variables
Check `aei/params.py` for the environment variables to set. If all of the executables for the above binary packages are accessible from the command line, no additional changes are necessary.
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
