#!/usr/bin/env python
# 
# use for all special tests except abirules & buildsys
#
import sys,os
from glob import glob
import shelve
import re
from datetime import *
try:
	import subprocess
except:
	pass

from os.path import join as pj, abspath as absp, exists as pexists, isfile, basename

try:
   from ConfigParser import SafeConfigParser
except:
   from configparser import SafeConfigParser

basedir, x = os.path.split(absp(__file__))
testbot_cfg = pj(basedir, "testbot.cfg")
Results = dict()

##############################################################################
def _str2list(string): return [s.strip() for s in string.split(",") if s]

def ReadListOfTests(f):
	parser = SafeConfigParser()
	parser.read(testbot_cfg)

	t = parser.get("testbot","special")
	return _str2list(t)

def MakeReport(f):
	special = ReadListOfTests(f)
	if special != []:
		for s in special:
		
			try:
				fname="%s.status" % s
				with open(fname,'r') as f:
					tmp=f.readline()
					Results[s] = tmp.rstrip('\n')
			except:
				if s == "distcheck":
					Results[s] =  "SKIP"
					break
				cmd="grep fail %s/tmp*/report; exit 0" % s
				p = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
				if p == "":
					Results[s] = "OK"
				else:
					Results[s] =  "FAIL"

	try:
		fname = "testbot_summary.shelve"
		d = shelve.open(fname,writeback=True)
		d['special'] = Results
		d['rawdate'] = datetime.now()
		d.close()
	except:
		print("no testbot_summary.shelve")
		pass
	print("Results : ",Results)

def GetListFailed(t):
	listfailed=[]
	fname=glob("%s/tmp*/report" % t)[0]
	#print fname
	with open(fname,'r') as f:
		# cp report for reporting
		freport="%s/report.log" % t
		with open(freport,'w') as fr:
			fr.writelines(f.readlines())
	with open(fname,'r') as f:
		# check if failed test
		for line in f:
			if re.search('failed',line):		    	
				m = re.split(r'^Case_(\d*)',line)
				#print(m)
				listfailed.append(m[1])
	#print(listfailed)
	return listfailed

def BuildsysErrorLogs():
        listfailed=GetListFailed('buildsys')
        #print(listfailed)
        fout="buildsys/out.log"
        ferr="buildsys/err.log"
        with open(fout,'w') as fo:
          with open(ferr,'w') as fl:
            for l in listfailed:
              f_out=glob("buildsys/tmp*/t%s.out" % l)[0]
              f_err=glob("buildsys/tmp*/t%s.err" % l)[0]
              with open(f_out,'r') as f:
                fo.write("test : t%s\n\n" % l)
                fo.writelines(f.readlines())
                fo.write("\n")
              with open(f_err,'r') as f:
                fl.write("test : t%s\n\n" % l)
                fl.writelines(f.readlines())
                fl.write("\n")
        return

def AbirulesErrorLogs():
        listfailed=GetListFailed('abirules')
        #print listfailed
        fdiff="abirules/diff.log"
        fout="abirules/out.log"
        flog="abirules/log.log"
        with open(fdiff,'w') as fi:
          with open(fout,'w') as fo:
             with open(flog,'w') as fl:
               for l in listfailed:
                f_out=glob("abirules/tmp*/t%s.out" % l)[0]
                with open(f_out,'r') as f:
                        fo.write("test : t%s\n\n" % l)
                        fo.writelines(f.readlines())
                        fo.write("\n")
                f_diff=glob("abirules/tmp*/diff.t%s" % l)[0]
                with open(f_diff,'r') as f:
                        fi.write("test : t%s\n\n" % l)
                        fi.writelines(f.readlines())
                        fi.write("\n")
                try:
                        f_log=glob("abirules/tmp*/t%s.log" % l)[0]
                        with open(f_log,'r') as f:
                                fl.write("test : t%s\n\n" % l)
                                fl.writelines(f.readlines())
                                fl.write("\n")
                except:
                        pass
        return

#############################################################################################################

if __name__ == "__main__":

	if len(sys.argv) == 1:
		MakeReport(testbot_cfg)
	else:
		test=sys.argv[1:2][0]
		if test == 'abirules':
			AbirulesErrorLogs()
		else:
			BuildsysErrorLogs()

sys.exit()
