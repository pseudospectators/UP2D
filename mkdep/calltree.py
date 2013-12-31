import sys,pickle

def usage():
    print """
pretty print the possible call trees over a given routine name
from the routines listed in .mkdep_routines

usage:
   calltree ROUTINE
"""
    sys.exit(2)

def print_names(names, dict, level, stack=[]):
    level += "    "
    
    for name in names:
        file= ""
        if dict.has_key(name):
            file = dict[name][0]            
        print (level+name).ljust(60) + file
        if name in stack:
            print level + "Recursion?"
            return        
        stack.append(name)        
        if dict.has_key(name):
            print_names(dict[name][1], dict, level, stack)
        last = stack.pop()

        
# get the routines dictionary from the restart file
f=open('.mkdep_restart_file','r')
for fivetimes in [1,2,3,4,5]:
    routines = pickle.load(f)
f.close()

if len(sys.argv)<2:
    usage()

print "In the trees of the current mkdep run (.mkdep_restart_file)"
print "the possible tree over "+sys.argv[1]+" looks like the below;"

level=""
stack=[]
print_names([sys.argv[1]], routines, level, stack)
