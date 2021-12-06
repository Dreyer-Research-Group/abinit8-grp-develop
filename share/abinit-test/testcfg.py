#!/usr/bin/env python

import sys,os
from pymongo import MongoClient

builder = sys.argv[1]

def fix_unicode(data):
    if isinstance(data, unicode):
        return data.encode('utf-8')
    elif isinstance(data, dict):
        data = dict((fix_unicode(k), fix_unicode(data[k])) for k in data)
    elif isinstance(data, list):
        for i in xrange(0, len(data)):
            data[i] = fix_unicode(data[i])
    return data

data={}
uri="mongodb://bbro:bbro@gitlab.pcpm.ucl.ac.be/buildbot"
client = MongoClient(uri)
db_builders = client.buildbot.builders
data =  db_builders.find_one({"name":builder})

if sys.version_info[:2] < (3,0,0):
    data = fix_unicode(data)  # convert unicode to string

with open('testbot.cfg','w') as f:
   f.write("[testbot]\n")
   slavename = data["name"]
   f.write("slavename = %s\n" % slavename )
   type = data["type"]
   if type != []:
       if 'ref' in type:
          f.write("type = ref\n")
   f.write("ncpus = %u\n" % data["ncpus"])
   if data["mpi_prefix"] != None :
       mpi_prefix = ""
       try:
          mpi_prefix = os.environ['MPI_HOME']
       except:
          pass
       if mpi_prefix == "":
          try:
             mpi_prefix = os.environ['MPIHOME']
          except:
             mpi_prefix = data["mpi_prefix"]
       f.write("mpi_prefix = %s\n" % mpi_prefix )
       f.write("mpirun_np = %s%s\n" % (mpi_prefix, data["mpirun_np"]) )
   if data["mpi_flavor"] == 'poe':
       f.write("poe = %s\n" % data["poe"])
       f.write("poe_args = %s\n" % data["poe_args"])
   if data["keywords"] !=  None:
       f.write("keywords = %s\n" % data["keywords"])
   if data["with_tdirs"] !=  None:
       f.write("with_tdirs = %s\n" % data["with_tdirs"])
   if data["without_tdirs"] !=  None:
       f.write("without_tdirs = %s\n" % data["without_tdirs"])
   f.write("omp_num_threads = %u\n" % data["omp_num_threads"])
   f.write("special = %s\n" % data["special"])
   f.write("fallbacks_prefix = %s\n" % data["fallbacks_prefix"])
   f.write("timeout_time = %d\n" % data["timeout_time"])
   f.write("runmode = %s\n" % data["runmode"])
   f.write("cygwin_dir = %s\n" % "")
   f.write("\n")
