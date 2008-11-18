from distutils.core import setup, Extension

root = "../.."


# Parse FLIBS environment variable, if it exists, otherwise
# default to gfortran.
extralink,fdirs,flibs = [],[],['gfortran']
import os
if 'FLIBS' in os.environ:
    fdirs,flibs = [],[]
    for str in os.environ['FLIBS'].split():
        if str[:2] == '-L':
            fdirs.append(str[2:])
        elif str[:2] == '-l':
            flibs.append(str[2:])
        else:
            extralink.append(str)
            print "Warning: unknown FLIBS component "+str

#cxxlibs = ['supc++','stdc++']
cxxlibs = []
reflmodel = Extension('reflmodel',
                   ['pysetup.cc'],
                   extra_compile_args=[],
                   extra_link_args=extralink,
                   include_dirs=[root+"/model1d/lib"],
                   library_dirs=[root+"/model1d/lib"]+fdirs,
                   libraries=['refl']+flibs+cxxlibs)

setup(name='reflmodel',
            version='1.0',
            py_modules=[],
            ext_modules=[reflmodel],
            )

__id__ = "$Id: setup.py 2 2005-08-19 17:42:06Z pkienzle $"
