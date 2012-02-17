#!/usr/bin/env python

# Andrey Asadchev <andrey@si.fi.ameslab.gov>.
# Python program to submit GAMESS jobs to Blue Gene/L.
# The program can submit to Blue Gene/L via Cobalt queueing system or
# via mpirun.  The program can be extended to use LoadLever as well.
#
# Set SYSTEM variable to either COBALT or MPIRUN and define appropriate
# command; cqsub for COBALT and mpirun for MPIRUN.
# Define exe, scr, envGamessEXTBAS, envGamessERICFMT, and envGamessMCPPATH
# variables.
# The rest of the variables can be left as it is.
# 

# Define SYSTEM, either COBALT or MPIRUN
SYSTEM = "COBALT"
cobalt="cqsub"   # Cobalt submit command
mpirun="mpirun"  # mpirun command

# Main GAMESS parameters
exe="gamess.00.x"                       # GAMESS executable
cwd="~/scr"                             # Scratch directory
envGamessEXTBAS = "~/bin/basis.txt"     # Extended basis set
envGamessERICFMT = "~/bin/ericfmt.dat"  # ERIC formatted data
envGamessMCPPATH = "~/bin/mcpdata"      # MCP data

# Default runtime system parameters
dest = "default"
nodes = ""
processors = ""
ppn = "2"
mode = "co"
time = ""
environ=""

# Default DDI parameters
gddi = False
header = False

# Default verbosity level
verbose = False

# GAMESS environment
envGamess = {"IRCDATA":"irc", "INPUT":"F05", "PUNCH":"dat",
             "AOINTS":"F08", "MOINTS":"F09", "DICTNRY":"F10",
             "DRTFILE":"F11", "CIVECTR":"F12", "CASINTS":"F13",
             "CIINTS":"F14", "WORK15":"F15", "WORK16":"F16",
             "CSFSAVE":"F17", "FOCKDER":"F18", "WORK19":"F19",
             "DASORT":"F20", "DFTINTS":"F21", "DFTGRID":"F22",
             "JKFILE":"F23", "ORDINT":"F24", "EFPIND":"F25",
             "PCMDATA":"F26", "PCMINTS":"F27", "SVPWRK1":"F26",
             "SVPWRK2":"F27", "MLTPL":" F28", "MLTPLT":"F29",
             "DAFL30":"F30", "SOINTX":"F31", "SOINTY":"F32",
             "SOINTZ":"F33", "SORESC":"F34", "SIMEN":"simen",
             "SIMCOR":"simcor", "GCILIST":"F37", "HESSIAN":"F38",
             "SOCCDAT":"F40", "AABB41":"F41", "BBAA42":"F42",
             "BBBB43":"F43", "MCQD50":"F50", "MCQD51":"F51",
             "MCQD52":"F52", "MCQD53":"F53", "MCQD54":"F54",
             "MCQD55":"F55", "MCQD56":"F56", "MCQD57":"F57",
             "MCQD58":"F58", "MCQD59":"F59", "MCQD60":"F60",
             "MCQD61":"F61", "MCQD62":"F62", "MCQD63":"F63",
             "MCQD64":"F64", "NMRINT1":"F61", "NMRINT2":"F62",
             "NMRINT3":"F63", "NMRINT4":"F64", "NMRINT5":"F65",
             "NMRINT6":"F66", "DCPHFH2":"F67", "DCPHF21":"F68" }

# GAMESS GDDI environment
envGamessGDDI = { "IRCDATA":"F04", "OUTPUT":"F06", "PUNCH":"F07" }


import sys
import getopt
import os
import shutil
import re
import math

def main():
    global dest, nodes, processors, ppn, mode, time, environ
    global exe, cwd
    global gddi, header
    global verbose
    
    shortOptions = "d:n:p:m:t:e:C:X:GHvh"
    longOptions = ["dest=", "nodes=", "processors=", "mode=", "time=", "environ", "cwd",
                   "exe=", "gddi", "header", "verbose", "help"]
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], shortOptions, longOptions)
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    # input files have to be specified
    if not args:
        usage()
        sys.exit(1)
        
    # validata system
    validateSystem()

    # initialize environ and output
    environ = environ.strip()
    output = ""
    
    for o, a in opts:
        if o in ("-d", "--dest"): dest = a
        if o in ("-n", "--nodes"): nodes = a
        if o in ("-p", "--processors"): processors = a
        if o in ("-m", "--mode"): mode = a
        if o in ("-t", "--time"): time = a
        if o in ("-e", "--environ"):
            if a.strip():
                if environ: environ = environ + ":" + a
                else: environ = a
        if o in ("-C", "--cwd"): cwd = a
        if o in ("-O", "--output"): output = a
        if o in ("-X", "--exe"): exe = a
        if o in ("-G", "--gddi"): gddi = True
        if o in ("-H", "--header"): header = True
        if o in ("-v", "--verbose"): verbose = True
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)

    # validate global arguments
    validateGlobalArguments()

    # set DDI header if needed.
    if header is True:
        envvar = "DDI_HEADER=RUNTIME"
        if environ: environ = environ + ":" + envvar
        else: environ = envvar

    # resolve paths
    exe = getPath(exe)
    cwd = getPath(cwd)
    output = getPath(output)

    for inputFile in args:
        envJob = ""
        # setup job
        try: (cwdJob, outputJob, envJob) = setupGamessJob(inputFile, cwd, output)
        except Exception,e :
            print e
            continue
        # configure environmental variables
        if not envJob: envJob = environ
        elif environ: envJob = environ + ":" + envJob
        # run job
        jobCmd = getJobCmd(dest, nodes, processors, ppn, mode, time, envJob, exe, cwdJob, outputJob, verbose)
        try:
            if verbose: print jobCmd
            os.system(jobCmd)
        except Exception, e: print e

    sys.exit(0)

# validate SYSTEM
def validateSystem():
    global SYSTEM
    global cobalt
    global mpirun

    # Validate SYSTEM
    if SYSTEM == "COBALT":
        if not cobalt:
            print "Cobalt command must be defined in the" + sys.argv[0]
            sys.exit(1)
    elif SYSTEM == "MPIRUN":
        if not mpirun:
            print "mpirun command must be defined in the" + sys.argv[0]
            sys.exit(1)

# validate arguments
def validateGlobalArguments():
    global nodes, processors, ppn, mode, time, exe

    nodes = nodes.strip()
    processors = processors.strip()
    ppn = ppn.strip()
    mode = mode.strip()
    time = time.strip()
    exe = exe.strip()

    if not isBlankOrInt(nodes):
        print "Number of nodes must be an integer."
        sys.exit(1)
        
    if not isBlankOrInt(processors):
        print "Number of processors must be an integer."
        sys.exit(1)

    if not isBlankOrInt(ppn):
        print "Number of processors per node must be an integer."
        sys.exit(1)

    if not isBlankOrInt(time):
        print "Time must be an integer."
        sys.exit(1)

    if mode != "co" and mode != "vn":
        print "Mode must be defined either as co or vn"
        sys.exit(1)

    if not exe:
        print "Executable must be defined in the" + sys.argv[0]
        print "or as the command line argument to -X"
        sys.exit(1)

# return true if string is blank or an integer
def isBlankOrInt(string):
    if string:
        try: int(string)
        except ValueError, e:
            return False
    return True
        
    
# expand variables and resolve to absolute path
def getPath(path):
    path = path.strip()
    path = os.path.expanduser(path)
    path = os.path.expandvars(path)
    path = os.path.abspath(path)
    return path
    
# set up GAMESS job.  I/O errors should be handled by the caller.
def setupGamessJob(inputFile, cwd, output):
    envFile = ""
    
    # determine job, job directory, job file
    job = getJob(inputFile)
    jobDirectory = getPath(os.path.join(cwd, job))
    jobFile = os.path.join(jobDirectory, job + ".F05")
    
    # determine job output
    if os.path.isdir(output):
        output = os.path.join(output, job)
        
    # erase any previous job directory
    if os.path.isdir(jobDirectory):
        shutil.rmtree(jobDirectory)
        
    # create new job directory, job file, and set GAMESS environment.
    try:
        os.mkdir(jobDirectory)
        shutil.copyfile(inputFile, jobFile)
        envFile = setEnvGamess(jobDirectory, job, isGDDI(jobFile, gddi))
    except Exception, e:
        try:
            shutil.rmtree(jobDirectory)
        except Exception, e:
            pass
        raise
            
    # append necessary environmental variables
    environ = "DDI_SCRATCH=" + jobDirectory + ":" + "DDI_CWD=" + jobDirectory
    if envFile: environ = environ + ":" + "ENVFIL" + "=" + envFile

    return (jobDirectory, output, environ)

# return job name.
def getJob(inputFile):
    job = os.path.basename(inputFile).split(".")
    return ".".join(job[0:len(job)-1])

# determine if the jobFile is GDDI.
def isGDDI(jobFile, gddi):
    pattern = re.compile("^ \$GDDI\s+", re.IGNORECASE)
    jobFileHandle = open(jobFile, "r")

    # Force GDDI
    if gddi is True: return True

    # Determine from the jobFile
    for line in jobFileHandle:
        if pattern.match(line):
            jobFileHandle.close()
            return True
    jobFileHandle.close()
    return False

# set GAMESS environmental variables
def setEnvGamess(jobDirectory, job, isGDDI):
    envFile = os.path.join(jobDirectory, job + ".env")
    envFileHandle = open(envFile, 'w')
    envGamess = getEnvGamess(jobDirectory, job, isGDDI)
    for k in envGamess:
        envFileHandle.write(k + "=" + envGamess[k] + "\n")
    envFileHandle.close()
    return envFile

# return dict of GAMESS environmental variables
def getEnvGamess(jobDirectory, job, isGDDI):
    global envGamess
    global envGamessGDDI
    global envGamessEXTBAS
    global envGamessERICFMT
    global envGamessMCPPATH

    _envGamess = {}
    _envGamessEXTBAS = envGamessEXTBAS.strip()
    _envGamessERICFMT = envGamessERICFMT.strip()
    _envGamessMCPPATH = envGamessMCPPATH.strip()

    # GAMESS environment
    if _envGamessEXTBAS: _envGamess['EXTBAS'] = getPath(_envGamessEXTBAS)
    if _envGamessERICFMT: _envGamess['ERICFMT'] = getPath(_envGamessERICFMT)
    if _envGamessMCPPATH: _envGamess['MCPPATH'] = getPath(_envGamessMCPPATH)
    
    # job environment
    for k in envGamess:
        v = envGamess[k].strip()
        if v: _envGamess[k] = os.path.join(jobDirectory, job + "." + v)

    # GDDI job environment
    if isGDDI:
        for k in envGamessGDDI:
            v = envGamessGDDI[k].strip()
            if v: _envGamess[k] = os.path.join(jobDirectory, job + "." + v)

    return _envGamess

# get job command to run
def getJobCmd(dest, nodes, processors, ppn, mode, time, environ, exe, cwd, output, verbose):
    if SYSTEM == "COBALT":
        return getJobCmdCobalt(dest, nodes, processors, ppn, mode, time, environ, exe, cwd, output, verbose)
    if SYSTEM == "MPIRUN":
        return getJobCmdMpirun(dest, nodes, processors, ppn, mode, time, environ, exe, cwd, output, verbose)

# get job command to run via Cobalt
def getJobCmdCobalt(dest, nodes, processors, ppn, mode, time, environ, exe, cwd, output, verbose):
    global cobalt

    # determine number of nodes
    if not nodes:
        if processors:
            if mode == "vn" and ppn: nodes = str(int(math.ceil(float(processors) / float(ppn))))
            else: nodes = processors
    
    cmd = cobalt
    if dest:       cmd = cmd + " -q " + dest
    if nodes:      cmd = cmd + " -n " + nodes
    if processors: cmd = cmd + " -c " + processors
    if mode:       cmd = cmd + " -m " + mode
    if time:       cmd = cmd + " -t " + time
    if environ:    cmd = cmd + " -e " + environ
    if cwd:        cmd = cmd + " -C " + cwd
    if output:     cmd = cmd + " -O " + output
    cmd = cmd + " " + exe
    
    return cmd

# get job command to run via mpirun
def getJobCmdMpirun(dest, nodes, processors, ppn, mode, time, environ, exe, cwd, output, verbose):
    global mpirun
    
    # determine number of processors
    if not processors:
        if nodes:
            if mode == "vn" and ppn: processors = str(int(nodes) * int(ppn))
            else: processors = nodes

    # stdout and stderr
    if not output:
        output = "mpirun"
    output = output + ".output"
    output = output + ".error"
    
    cmd = mpirun
    if dest:       cmd = cmd + " -partition " + dest
    if processors: cmd = cmd + " -np " + processors
    if mode:       cmd = cmd + " -mode " + mode.capitalize()
    if time:       cmd = cmd + " -timeout " + str(int(time) * 60)
    if environ:    cmd = cmd + " -env " + "\"" + environ.replace(":", " ") + "\""
    if cwd:        cmd = cmd + " -cwd " + cwd
    if verbose:    cmd = cmd + " -v " + "1"
    # run command in background with nohup
    cmd = "nohup " + cmd + " -exe " + exe + " > " + output + " 2> " + error + " &"

    return cmd
    
# print usage
def usage():
    print "Usage: " + sys.argv[0] + " [OPTION]... FILE [FILE]..."
    print "Submit GAMESS input FILEs to Blue Gene/L."
    print ""
    print "-d, --dest=DEST       Submit to queue/parition DEST."
    print "-n, --nodes=NN        Request NN nodes."
    print "-p, --processors=NP   Request NP processors"
    print "-m, --mode=MODE       Blue Gene/L mode."
    print "-t, --time=MIN        Request MIN minutes."
    print "-e, --environ=ENV     Pass column separated environmental"
    print "                      variables list ENV to the job runtime."
    print "-C, --cwd=DIR         Set working directory to DIR."
    print "-X, --gamess=EXE      Use GAMESS executable EXE"
    print "-G, --gddi            Force GDDI."
    print "-H, --header          Print DDI header."
    print "-v, --verbose         Verbose output."
    print "-h, --help            Display this help and exit."
    return

if __name__ == "__main__":
    main()
    
