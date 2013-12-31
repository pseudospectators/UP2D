# a list over the different suffixes that are in use
suffixes=[]         

# a dictionary where the flag/group is the key, and the data is a list
# of the object files in that group
object_files={}

# now, the new source file data structure stored in .mkdep.map, and is
# reloaded when mkdep is executed in order for mkdep to remember
# things about the last run
#
# the keys are full paths to source files, the data is a tuple
#  (type, is_include, compflag, mtime, deps, reparse, deleted)
#   type = free, fixed, fixed132 etc
#   is_include is 0 for source files
#   compflag is ?
#   mtime = time for last parse of this source file
#   deps = list of object files that this file depends on
#   reparse = 1 if this file should be reparsed.
#   deleted = 1 if this file is missing from the new project file list.

clean_source_files2={}         

# the mod time for this file
project_mtime=0                

# a map from module to tuple (file name, [reference symbols], [provided symbols])
module_map = {}     

# a dictionary where the key is the routine name in upper case, and
# the data is a tuple with the file name, a list of other routines
# called and the host routines. also contain all function and program
# names.
routines={}         

# a dictionary like above, but only with functions
functions={}        

# like above, but only containing the programs.
programs={}

incdirs = []

calltree=0
prevcalltree=0

incpathstring=""

# true if we want to prepend "build" to object file locations:
use_build_dir=1
