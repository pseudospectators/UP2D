import sys,pickle


def usage():
    print """
pretty print the possible call stacks for a given routine name
from the routines listed in .mkdep_routines

usage:
   callstack ROUTINE
"""
    sys.exit(2)
    
def write_callees(name, level="", called=True):
    """
    print the given routine name with the given indentation level
    find its parent and call self with parent as arg.
    """
    if called:
        print (level+name).ljust(30)+" "+routines[name][0]
    else:
        print (level+name+"*").ljust(30)+" "+routines[name][0]
    level += "   "
    for key in routines.keys():
        # check if someone called routine "key"
        #if len(routines[key][2]) > 0:
        # if so, check if "name" is called from routine "key"
        if name in routines[key][1]:
            # if so, "key" is a dad, indent and print
            if level=="":
                print ""
            if len(routines[key][2]) > 0:
                write_callees(key, level)
            else:
                write_callees(key, level, False)
   

# get the routines dictionary from the restart file
f=open('.mkdep_restart_file','r')
for fivetimes in [1,2,3,4,5]:
    routines = pickle.load(f)
f.close()

if len(sys.argv)<2:
    usage()

print "In the trees of the current mkdep run (.mkdep_restart_file)"
print "the possible callees of "+sys.argv[1]+" looks like the below;"
print "routines annotated * are not called from the given files.\n"

write_callees(sys.argv[1])
