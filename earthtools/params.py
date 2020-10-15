#####
# contains the parameters accessed by earthtools functions
#
# c. 2016-2017 Christopher Anderson
#####

import os
import sys
import platform
import threading
import subprocess
import multiprocessing
from psutil import virtual_memory

## get environment variables
environ = platform.os.environ
pathsep = platform.os.path.sep
pathcat = platform.os.path.join
system  = platform.system()

# check if LTBASE environment variable set
#  points to where lastools is installed
try:
    ltbase = environ['LTBASE']
except KeyError:
    ltbase = ''

# check if GDALBASE environment variable set
#  points to custom install of gdal binaries
try:
    gdalbase = environ['GDALBASE']
except KeyError:
    gdalbase = ''

# check if GDALBASE environment variable set
#  points to custom install of gdal binaries
try:
    otbbase = environ['OTB_BASE']
except KeyError:
    otbbase = ''

# check OUTPUT_DIR environment variable
#  points to where default outputs go if not set
try:
    outdir = environ['OUTPUT_DIR']
except KeyError:
    outdir = './'

# check SCRATCH_DIR environment variable
#  points to where intermediate files go if not set
try:
    scratchdir = environ['SCRATCH_DIR']
except KeyError:
    scratchdir = './'

cores = multiprocessing.cpu_count()
mem = virtual_memory().total

# list of gdal data types
gdalTypes = ["AAIGrid",
            "ACE2",
            "ADRG",
            "AIG",
            "AIRSAR",
            "ARG",
            "BLX",
            "BAG",
            "BMP",
            "BPG",
            "BSB",
            "BT",
            "CALS",
            "CEOS",
            "COASP",
            "COSAR",
            "CPG",
            "CTG",
            "DB2",
            "DDS",
            "DIMAP",
            "DIPEx",
            "DODS",
            "DOQ1",
            "DOQ2",
            "DTED",
            "E00GRID",
            "ECRGTOC",
            "ECW",
            "EHdr",
            "EIR",
            "ELAS",
            "ENVI",
            "EPSILON",
            "ERS",
            "ESAT",
            "FAST",
            "FIT",
            "FITS",
            "FujiBAS",
            "GENBIN",
            "GPKG",
            "GEORASTER",
            "GFF",
            "GIF",
            "GRIB",
            "GMT",
            "GRASS",
            "GRASSASCIIGrid",
            "GSAG",
            "GSBG",
            "GS7BG",
            "GSC",
            "GTA",
            "GTiff",
            "GTX",
            "GXF",
            "HDF4",
            "HDF5",
            "HF2",
            "HFA",
            "IDA",
            "ILWIS",
            "INGR",
            "IRIS",
            "ISCE",
            "ISIS2",
            "ISIS3",
            "JAXAPALSAR",
            "JDEM",
            "JPEG",
            "JPEGLS",
            "JPEG2000",
            "JP2ECW",
            "JP2KAK",
            "JP2MrSID",
            "JP2OpenJPEG",
            "JPIPKAK",
            "KEA",
            "KMLSUPEROVERLAY",
            "KRO",
            "L1B",
            "LAN",
            "LCP",
            "Leveller",
            "LOSLAS",
            "MBTiles",
            "MAP",
            "MEM",
            "MFF",
            "MFF2 (HKV)",
            "MG4Lidar",
            "MRF",
            "MrSID",
            "MSG",
            "MSGN",
            "NDF",
            "NGSGEOID",
            "NITF",
            "netCDF",
            "NTv2",
            "NWT_GRC",
            "NWT_GRD",
            "OGDI",
            "OZI",
            "PAux",
            "PCIDSK",
            "PCRaster",
            "PDF",
            "PDS",
            "PLMosaic",
            "PNG",
            "PostGISRaster",
            "PNM",
            "R",
            "RASDAMAN",
            "Rasterlite",
            "RIK",
            "RMF",
            "ROI_PAC",
            "RPFTOC",
            "RS2",
            "RST",
            "SAFE",
            "SENTINEL2",
            "SAGA",
            "SAR_CEOS",
            "SDE",
            "SDTS",
            "SGI",
            "SNODAS",
            "SRP",
            "SRTMHGT",
            "TERRAGEN",
            "TIL",
            "TSX",
            "USGSDEM",
            "VICAR",
            "VRT",
            "WCS",
            "WEBP",
            "WMS",
            "WMTS",
            "XPM",
            "XYZ",
            "ZMap"]

# define ogr types
ogrTypes = ["AeronavFAA",
            "AmigoCloud",
            "ArcObjects",
            "AVCBin",
            "AVCE00",
            "ARCGEN",
            "BNA",
            "DWG",
            "DXF",
            "Carto",
            "Cloudant",
            "CouchDB",
            "CSV",
            "CSW",
            "VFK",
            "DB2ODBC",
            "DODS",
            "EDIGEO",
            "ElasticSearch",
            "FileGDB",
            "PGeo",
            "SDE",
            "ESRI Shapefile",
            "FMEObjects Gateway",
            "GeoJSON",
            "Geoconcept",
            "Geomedia",
            "GPKG",
            "GeoRSS",
            "GFT",
            "GML",
            "GMT",
            "GPSBabel",
            "GPX",
            "GRASS",
            "GPSTrackMaker",
            "HTF",
            "Idrisi",
            "IDB",
            '"Interlis 1" and "Interlis 2"',
            "INGRES",
            "OpenJUMP .jml",
            "KML",
            "LIBKML",
            "MapInfo File",
            "DGN",
            "MDB",
            "Memory",
            "MongoDB",
            "MySQL",
            "NAS",
            "netCDF",
            "OCI",
            "ODBC",
            "MSSQLSpatial",
            "ODS",
            "OGDI",
            "OpenAir",
            "OpenFileGDB",
            "OSM",
            "PCIDSK",
            "PDF",
            "PDS",
            "PLScenes",
            "PGDump",
            "PostgreSQL/PostGIS",
            "REC",
            "S57",
            "SDTS",
            "SEGUKOOA",
            "SEGY",
            "Selafin",
            "SOSI",
            "SQLite",
            "SUA",
            "SVG",
            "SXF",
            "UK. NTF",
            "TIGER",
            "VRT",
            "WFS",
            "XLS",
            "XLSX",
            "XPLANE",
            "VDV",
            "Walk",
            "WAsP"]
