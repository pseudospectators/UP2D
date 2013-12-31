"""
some utilities for mkdep.py
Helge Avlesen <avle@bccs.uib.no>  para//ab
"""

import sys, string, os, re, pickle, commands

debug=0
insert_group=0
full_classify=0
dump_modmap=0
treetop=""
show_warnings=0
calltree_title = "mkdep calltree"
root_path = "./"
assemblefile=""

import globals

# matcher for a valid and upper case Fortran name
match_fortran_name = re.compile("[A-Z][A-Z0-9_]*")
match_single_quote = re.compile("'")
match_double_quote = re.compile('"')
match_exclam       = re.compile("!")


def classify_files(projfile='', incfile='', rootdir=''):
    """
    go through the list of file names, assign default file types by suffix,
    look for duplicate files and include directories.
    """

    ## the default guesses for file types
    
    filetype = { ".F":"fixed", ".F90":"free", ".f":"fixed", \
                 ".f90":"free", ".c":"C", ".h":"C", ".H":"C",\
                 ".inc":"fixed", ".INC":"fixed", ".SYM":"fixed"}

    ## argument 1 should be a file with a list of Fortran files
    try:
        inputfile = open(projfile,'r')
    except IOError:
        print "Could not open project file called ",projfile
        sys.exit(2)
        
    source_files = inputfile.readlines()
    inputfile.close()

    globals.project_mtime = os.stat(projfile).st_mtime

    load_maps()

    longnamesize=20
    shortnamesize=12
    shortnamelist={}
    dups=0

    rs = rootdir.strip()
    if len(rs)>0:
        if rs[-1] in ["/","\\"]:
            rootdir = rs[:-1]
        else:
            rootdir = rs

    ## edit the list of source files by adding default values to
    ## missing columns for language, file type, optimization.
    ## leave comments or blank lines.

    if insert_group or full_classify:
        newsrcfile=open(projfile,'w')

    if debug:
	classifylist = open('.mkdep_filelist','w')

    if assemblefile <> "":
        assem=open(assemblefile, 'w')
        
    line_number=0

    for line in source_files:
        line_number+=1
        current_line = line.strip()

        ## skip blank or comment lines
        if current_line=="":
            if insert_group or full_classify: newsrcfile.write(line)
            continue
        if current_line[0]=='#':
            if insert_group or full_classify: newsrcfile.write(line)
            continue

        ## if the previous parsing was without parsing, force a full reparse
        reparse=0
	if globals.calltree and globals.prevcalltree==0:
            reparse=1

       
        ## keep comments at end of line
        comment=""
        cpos=current_line.find("#")
        if cpos>-1:
            comment=current_line[cpos:]
            current_line=current_line[:cpos]

        current_line = current_line.split()
        cols=len(current_line)

        if rootdir=="":
            longname = current_line[0]
        else:
            longname = rootdir.strip() + os.path.sep + current_line[0]

        try:
            statobj = os.stat(longname)
        except OSError:
            print "mkdep error: please check your file list, could not find file:",longname
            sys.exit()
        
        if globals.clean_source_files2.has_key(longname):
            # source file is modified after last parse
            if statobj.st_mtime > globals.clean_source_files2[longname][3]:
                reparse = 1
                
        shortname = longname.split("/")[-1].strip()
        dotpos = shortname.rfind(".")
        suffix = shortname[dotpos:]


        ## the default group for object files
        group_name="OPT"
        if cols>1:
            group_name=current_line[1]
        do_opt=group_name

        if cols>2:
            type=current_line[2]
            if type not in ['free','fixed','fixed132','C']:
                print 'line',line_number,'column 3 "'+type+\
		      '" should preferrably be one of fixed,fixed132,C or free:'
                print line,
        else:
            if filetype.has_key(suffix):
                type=filetype[suffix]
            else:
                type="unknown"

        # if the file type has changed, we may need to parse in a different way
        if globals.clean_source_files2.has_key(longname):
            if type <> globals.clean_source_files2[longname][0]:
                reparse=1
            
                
        inc="source"
        is_include=0
        if cols>3:
            inc=current_line[3]
            if inc not in ['include','source']:
                print 'line',line_number,'column 4 "'+inc+'" must be one of include,source:'
                print line,

            if inc=="include": is_include=1
        else:
            if suffix.lower() in [".h",".i",".inc",".com",".sym"]:
                is_include=1
                inc="include"
                
        ## assemble source files into one
        if not is_include and assemblefile <> "":
            tmpfile = open(longname,'r')
            assem.write( tmpfile.read() )
            tmpfile.close()            

            
        # if is_include has changed, and it has changed to 0, reparse
        if globals.clean_source_files2.has_key(longname):
            if is_include <> globals.clean_source_files2[longname][1] and not is_include:
                reparse=1
                
        slashpos = longname.rfind("/")
        dirname = longname[:slashpos]


        ## make a suggestion for a list of include dirs
        if is_include:
            if dirname not in globals.incdirs:
                globals.incdirs.append(dirname)

        if cols>4:
            themodifiedline = line
	    fullmodifiedline = line
        else:
            if len(longname) >= longnamesize: longnamesize = len(longname) + 2
            themodifiedline = longname.ljust(longnamesize) + \
                              group_name.ljust(10) + \
			      comment + "\n"
            fullmodifiedline = longname.ljust(longnamesize) + \
                              group_name.ljust(10) + \
			      type.ljust(9) + \
			      inc.ljust(8) + \
			      comment + "\n"
        if insert_group:
	    newsrcfile.write(themodifiedline)
	elif full_classify:
	    newsrcfile.write(fullmodifiedline)
	if debug:
	    classifylist.write(fullmodifiedline)

        if reparse or not globals.clean_source_files2.has_key(longname):
            deps = []
            reparse = 1
        else:
            deps = globals.clean_source_files2[longname][4]
            
        statobj=os.stat(longname)

        # a file that will be processed, remove the default deleted status:
        deleted=0
        
        globals.clean_source_files2[longname] = (type, is_include, do_opt, statobj.st_mtime, deps, reparse, deleted)
        
        if shortnamelist.has_key(shortname):
            shortnamelist[shortname].append(longname)
            dups=1
        else:
            shortnamelist[shortname]=[longname]

    if assemblefile <> "":
        assem.close()
    
    ## delete files not in project file
    for key in globals.clean_source_files2.keys():
        if globals.clean_source_files2[key][6]:
            print "INFO:", key, " has apparently been removed from the list of source files"
            del globals.clean_source_files2[key]

            for list in globals.object_files.keys():
                try:
                    globals.object_files[list].remove(key)
                except ValueError:
                    pass
                
            # also remove object file from dependency lists
            for key2 in globals.clean_source_files2.keys():
                dotpos = key.rfind(".")
                objf = key[:dotpos] + ".o"
                #print "obj to remove=",objf, "current key=", key2, "deps=", globals.clean_source_files2[key2][4]
                try:
                    globals.clean_source_files2[key2][4].remove(objf)
                except ValueError:
                    pass
           
            
    if insert_group: newsrcfile.close()

    
    ## if none of the source files are marked MAIN, assume the first
    ## file to contain the PROGRAM statement
    
    #opttypes = map(lambda x: x[3], clean_source_files)
    #try:
    #    opttypes.index('MAIN')
    #except ValueError:
    ## since we cannot modify tuples we gotta do this...        
    #    first = clean_source_files[0]
    #    clean_source_files[0] = (first[0],first[1],first[2],'MAIN')


    if dups:
        print "\nIn the given list of source files, some of the file names were"
        print "identical, which indicates a duplicate. See file .mkdep_duplicates"
        print "for a list over files with identical names.\n"
        snf = open('.mkdep_duplicates','w')
        for key in shortnamelist:
            if len(shortnamelist[key])>1:
                if len(key)>shortnamesize: shortnamesize = len(key)+2
                snf.write(key+": " + string.join(shortnamelist[key]) + "\n")
        snf.close()

    ## if argument 2 is present, it is a file with a list of include directories
    ## we also write a suggestion for this list, based on the info from sys.argv[1]

    if len(globals.incdirs)>0:
        ifile = open(".mkdep_suggested_incdirs","w")
        ifile.write("# This is a list over directories where '"+projfile+"'\n")
        ifile.write("# suggest we have an include file:\n\n")
        ifile.write(string.join(globals.incdirs,"\n"))
        ifile.close()


    ifile = open(".mkdep_includes","w")
    ifile.write("INC=")
    includepath=[]

    if incfile<>'':
        inputfile = open(incfile,'r')
        lines=inputfile.readlines()
        lines.reverse()

        for line in lines:
            line=line.strip()
            if line=="":
                continue
            if line[0]=="#":
                continue

            if rootdir=="":
                line = line.strip()
            else:
                line = rootdir.strip() + os.path.sep + line.strip()

            includepath.append(line)
            globals.incpathstring += " -I" + line

        ifile.write(globals.incpathstring)
        print "Include files will be picked up in this order:"
        print globals.incpathstring + "\n"
        
    elif len(globals.incdirs)>0:
        print """
Info: you did not supply a file with a list of include directories to
search. It may be necessary to provide include paths (-Idir1 -Idir2 -I...)
in the Makefile manually. You can use the file .mkdep_suggested_incdirs
as a template for the file to use with the -i option.
"""        
    ifile.write("\n")
    ifile.close()

    return includepath


#def system_command(command):
#    pipe = os.popen(command,'r')
#    text = pipe.read()
#    stat = pipe.close()
#    return stat, text


def getsrc(cppcommand, cppflags, filename):
    """
    fetch the source from "filename" into a string
    """
    if cppcommand=="":
        tmpfile = open(filename,'r')
        return tmpfile.read()
    
    cpp = cppcommand + " " + cppflags + " " + filename

    stat, out = commands.getstatusoutput(cpp)
    
    if len(out) < 3:
        print "cpp output looks too small:", out
        print "cpp command attempted:", cpp
        sys.exit(2)

    return out


def strip_comments(line):
    """
    line is one line of free format Fortran. return the line
    stripped for comments as well as leading and trailing white
    space.
    """
    line = line.lstrip()             # remove leading spaces
    ccol = line.find("!")            # find position of first !
    # TODO: ! could be inside a string...
    if ccol<0:
        return line
    if ccol==0:
        return ""                    # it was a comment line
    if ccol>0:
        return line[:ccol].rstrip()  # remove everything after !

def strip_qfcomments(line):
    """
    line is one line of free format Fortran. return the line
    stripped for comments as well as leading and trailing white
    space.
    """
    ccol = line.find("!")            # find position of first !
    # TODO: ! could be inside a string...
    if ccol<0:
        return line
    if ccol==0:
        return ""                    # it was a comment line
    if ccol>0:
        return line[:ccol].rstrip()  # remove everything after !

def match_module(names):
    """
    line is a stripped Fortran statement. (no comment or
    white space at beginning or end.)  if it contains the beginning of
    a module definition, return the module name in upper case
    """
    if len(names) < 2: return ""
    if names[0] == "MODULE" and names[1] <> "PROCEDURE":
        return names[1]
    return ""

def match_module_end(names):
    """
    """
    if len(names) < 3: return ""
    if names[0]=="END" and names[1]=="MODULE":
        return names[2]
    return ""

def match_end_subprogram(names):
    """
    """
    if len(names) > 2:
        if names[0]=="END" and names[1] in ["SUBROUTINE","FUNCTION","PROGRAM"]:
            return names[2]
    if len(names) == 1:
        if names[0] == "END":
            return "OLDSTYLE_END"
    return ""

def match_end(names):
    """
    """    
    if len(names) == 1:
	if names[0]=="END":
	    return 1
    return 0

def match_subroutine_end(names):
    """
    """
    if len(names) < 3: return ""   
    if names[0]=="END" and names[1]=="SUBROUTINE":
        return names[2]
    return ""

def match_function_end(line):
    """
    """
    if len(names) < 3: return ""
    if names[0]=="END" and names[1]=="FUNCTION":
        return names[2]
    return ""

def match_program_def(names):
    """
    """
    if len(names) < 2: return ""
    if names[0]=="PROGRAM":
        return names[1]
    return ""


def match_and_chop(line, pos, delimiter):
    """
    line has delimiter at pos, find its match, and return line
    without the contained string (also remove the delimiters)
    """
    mpos = line[pos+1:].find(delimiter)
    if mpos>0:
        if pos+1+mpos < len(line)-1:
            if line[pos+mpos+2]==delimiter:            
#                return line[:pos+mpos+1]+line[pos+mpos+3:]
                return line[:pos+mpos+1]+"A_STRING"+line[pos+mpos+3:]                
            else:
#                return line[:pos]+line[pos+1+mpos+1:]
                return line[:pos]+"A_STRING"+line[pos+1+mpos+1:]
        else:
#            return line[:pos]+line[pos+1+mpos+1:]
            return line[:pos]+"A_STRING"+line[pos+1+mpos+1:]
    elif mpos==0:
        return line[:pos] + "A_STRING" + line[pos+2:]
    else:
        print "error libmkdep.empty_string: only found one " + delimiter + " in line:"
        print line
        sys.exit(2)

def strip_strings(line):
    """
    line is a stripped (python wise) fortran statement.  return line, with 
    all fortran strings removed, delimiters included...
    """
    quotes_in_string=1
    while quotes_in_string:
        dqpos=line.find('"')
        sqpos=line.find("'")
        if dqpos<0 and sqpos<0:
            quotes_in_string=0
        else:
            if dqpos<0:
                line = match_and_chop(line, sqpos, "'")
            else:
                if sqpos<0:
                    line = match_and_chop(line, dqpos, '"')
                else:
                    if dqpos<sqpos:
                        line = match_and_chop(line, dqpos, '"')
                    else:
                        line = match_and_chop(line, sqpos, "'")
    return line

def match_function_def(names):
    """
    """
    # double precision function xxx
    if len(names) >=4:
        if names[2]=="FUNCTION":
            return names[3]    
    if len(names) >=3:
        if names[1]=="FUNCTION" and names[0]<>"END":
            return names[2]
    if len(names) >=2:
        if names[0]=="FUNCTION":
            return names[1]
    return ""

def match_subroutine_def(names):
    """
    """
    if len(names) < 2: return ""
    if names[0]=="SUBROUTINE":
        return names[1]
    return ""

def match_subroutine_call(names):
    """
    """
    if len(names) < 1: return ""
    try:
        i=names.index("CALL")
    except ValueError:
        return ""
    return names[i+1]

def match_function_call(line, names, functions):
    """
    """
    functions_called=[]
    for name in names:
        if functions.has_key(name):
            # ok, poss1ible call. but we must check wether it is a call
            # or a local variable. we assume functions are called with
            # an argument list. this is not lways the case.  if this
            # is not enough, we must know if there is a local variable
            # with the function name, which requires parsing...
            pattern = name + " *\("
            if len(re.findall(pattern, line)) >= 1:
                functions_called.append(name)
    return functions_called

def match_use(line, map, filename):
    """
    line is a stripped Fortran statement. If it is a use
    statement, return the filename containing the used module, if
    the module is found in the MODULE->FILENAME map.
    if not, we assume the module will be found elsewhere.
    
    also return a list of symbols made available with the use statement
    """
    mcol = line.find("USE")
    if mcol<>0:
        return "",[]
    if line[mcol+3]<>" " and line[mcol+3]<>",":
        return "",[]
    ccol = line.find(",")
    colcol = line.find(":")
    
# USE MODULE_NAME, ONLY : LOCAL_NAME1 => USE_NAME1, USE_NAME2, ...
	    
    if colcol+ccol >= 0:
        if colcol >= 0:
            # got "only", split the stuff after : on the ,        
            symlist=line[colcol+1:].split(",")

            # now handle renames.
            # rename_matcher = re.compile("[A-Z][A-Z0-9_]* *=> *[A-Z][A-Z0-9_]*")
            # for rename in rename_matcher.findall(line):
            #    from, to = rename.split('=>')
            #    ...
            
            #initially, only use the global name, if rename
            symlist2=[]
            for item in symlist:
                # split on rename, use the global sum (last of list)
                symlist2.append( item.split("=>")[-1] )
            del symlist
            symlist=symlist2

        if ccol >= 0:
            # got , split the stuff after the first , into substrings.
            # (if colcol, we must also have a , after the module name)
            name = line[4:ccol].strip()
            
            #symlist = line[ccol+1:].split(",")
            # later. ignore renames for now. add all provided symb to namespace
            if map.has_key(name):
                symlist = map[name][2]
       
    else:
        # only got a name, add all symbols to name space.
        name = line.rstrip()[4:]    # got the module name
    
    if map.has_key(name):
        return map[name][0],[]
    else:
        print "Warning: "+filename+": module "+name+" source file not given"
        return "",[]


def match_include(line, include_path=["."]):
    """
    match fortran and c include statements on a stripped line,
    return 1,filename for a match, 0,"" otherwise.
    """
    tmp=line
    icol = tmp.upper().find("INCLUDE ")
    if icol<0:
        return ""
    if icol>1 and line[0]=="#":   ## if col 0 not a #, no real include
        return ""

    ## BUG! we strip strings before we parse for include
    ## statements... this increadibly clever trick makes it impossible
    ## to get dependencies for include files right...
    
    ## locate the quotes, that can be ' or ". filename inbetween.
    firstquote = line[icol+8:].find('"')
    if firstquote>-1:
        lastquote = line[icol+8:].rfind('"')
    else:
        firstquote = line[icol+8:].find("'")
        if firstquote>-1:
            lastquote = line[icol+8:].rfind("'")
        else:
            return ""

    name=line[icol+9+firstquote:icol+8+lastquote]

    ## now locate the files...
    found=0
    
    for path in include_path:
        if path <> ".":
            testname = path + "/" + name
        else:
            testname = name
        if os.path.isfile(testname):
            name=testname
            found=1
            
    if not found:
        print "\nInfo: include file",name,"not found in include paths"
        name=""
    return name


def join_and_strip_free(line, joined_line):
    """
    line is a free format fortran line from the scanner.
    if line has an ampersand at the beginning, append line to joined line.
    if line has an ampersand at the end, return the joined line and a flag=1.
    otherwise flag=0.
    """
    ## strip off free form comments
    ## this is safe(r) if we have an even number of " or ' up until the !

    ## fails on
    ##   write(sf,'("!! -*-f90-*-",/,"!!",/,a)', advance='no') "!! Restart file written by BOM5 on "

    expos = line.find("!")


    if expos >= 0:
        # figure out if it is safe to run the brute force
        # strip_comments by counting quotes. if the string before the
        # ! has uneven number of single OR double quotes it is not
        # safe, because ! probably is in a string.
        ns1=count_single_quotes(joined_line)
        ns2=count_single_quotes(line[:expos])
        nd1=count_double_quotes(joined_line)
        nd2=count_double_quotes(line[:expos])

        unsafe = (ns1+ns2)%2 or (nd1+nd2)%2
        
        if not unsafe:            
            line = strip_comments(line)

        
    if line=="":
        return joined_line,1

    ## remove ampersand at beginning   
    line=line.lstrip()
    if line<>"":
	if line[0]=="&":
	    line=line[1:]

    ## then check for ampersand at the end
    amppos = line.rfind("&")

    if amppos>-1:
        ## if there is something on the line after the last amp. it could mean the amp is in a string.
	if line[amppos+1:].strip() == "":	
	    joined_line += line[:amppos].strip() + " "
	    return joined_line, 1
	else:
	    joined_line += line.strip()	    
    else:
        joined_line += line.strip()
    return joined_line, 0


def handle_last(line, long_line, next_line):
    warning_flag=0
    if line.strip()=="":
        if long_line <> "":
            return long_line, 0, "",warning_flag
        if next_line <> "":
            return next_line, 0, "",warning_flag
        warning_flag=1
        return "",0,"",warning_flag
    if line[0] in ["C","c","!","*"]:
       if long_line<>"":
            return long_line, 0, "",warning_flag
       if next_line<>"":
            return next_line, 0, "",warning_flag
    contline=0
    if len(line)>6:
        if line[5]<>" ":
            contline=1
    line=line[6:]
    if contline:
        if long_line <> "":
            long_line += line.strip()
        if next_line <> "":
            long_line = next_line + line.strip()
        return long_line, 0, "",warning_flag
    if not contline:
        if long_line<>"":
            return long_line + ";" + line.rstrip(), 0, "",warning_flag
        if next_line<>"":
            return next_line+";"+line.rstrip(), 0, "",warning_flag
    return line.rstrip(), 0, "", warning_flag


def count_quotes(string):
    """
    return the number of quotes in string
    """
    return len(match_single_quote.findall(string)) + len(match_double_quote.findall(string))

def count_single_quotes(string):
    """
    return the number of quotes in string
    """
    return len(match_single_quote.findall(string))

def count_double_quotes(string):
    """
    return the number of quotes in string
    """
    return len(match_double_quote.findall(string))
    

def strip_trailing_comment(line, n_prev=None):
    """
    line is the beginning of a valid fortran statement.
    we want to return the statement without the comment.
    we count ' and ", if we have an even number of both,
    we can be pretty sure a ! is a comment.

    if n_prev is present, it is the number of quotes from a previous
    string to be prepended
    """
    i=0     # position of current ! in the substring after the previous
    gi=0    # global position of previous !

    # when i is -1 we have tested all !'s
    while i >= 0:
        gi += i+1
        i = line[gi+1:].find("!")
        if i >=0:
            n = count_quotes(line[:gi+1+i])
            
            # if the string is a continuation, n_prev is the number of
            # quotes in the statement we are appending to
            if n_prev: n += n_prev
            
            if n>0:
                # if n is even (n%2==0) we can chop:
                if not n%2:
                    return line[:gi+1+i]
            else:
                return line[:gi+1+i]
    return line


def join_and_strip_fixed(line, long_line, next_line, last, type="fixed"):
    """
    fixed format is a bit more tedious than free form; when the scanner feeds us a
    line we dont know yet if the next line is a continuation line or not,

    we can therefore conclude that we have a complete line only when a new line is found.
    therefore the parsing of fixed form code is done with a one line lag
    w.r.t. the free form (that we started out with), so unfortunately
    this routine became a bit messy...

    line is the thing we get from the scanner.

    if this line is not a continuation line, long_line contains a line
    of fortran we can parse. we then signal that we want to parse
    long_line, and store "line" in "next_line"
    
    if line is a cont line, we must add it to long_line and say we
    need one more line. next_line is therefore set to "".

    return 1 as the second result if the scanner should fetch more lines, 0 if
    we have a complete line ready to be parsed.
    """

    # if 1st bit: tab on line, 2nd bit: possibly continuation over cpp
    warning_flag=0

    if line.find("\t")>-1:
        if debug:
            print "Warning: this line has a TAB:"
            print line
        else:
            warning_flag |= 2**0

    if debug:
        print "1>"+line.replace("\n","")+"<"
        print "2>"+long_line.replace("\n","")+"<"
        print "3>"+next_line.replace("\n","")+"<"

    ## special treatment of the last line in the file
    if last:
        a1,b1,c1,d1 = handle_last(line, long_line, next_line)
        return a1,b1,c1,d1

    ## if line is blank we should check if next_line or long_line has
    ## something to be written: if long_line is <>"" we can send it to
    ## splitting, and set next_line="". next_line and long_line cannot
    ## normally be "" simultaneously since long_line is "" after the
    ## splitter. if next_line is<>"" while line is blank it means we
    ## can send it to the splitter

    if line.strip()=="":
        if long_line <> "":
            return long_line, 0, "",warning_flag
        if next_line <> "":
            return next_line, 0, "",warning_flag
        return "", 1, "",warning_flag

    ## if we are a comment line, we may need to join more

    if line[0] in ["C","c","!","*"]:
        return long_line, 1, next_line,warning_flag

    ## (cpp handling goes here if needed)
    
    ## strip stuff after col 72 for pure fixed form
    if type=="fixed":
        if len(line)>72:
            line=line[:73].rstrip()

    ## decide if "line" is a new line or a continuation line
    contline=0
    if len(line)>6:
        if line[5]<>" ":
            contline=1
            
    ## at the very beginning of files, this may happen:
    if long_line==next_line=="":
        if line.strip()=="":
            return "",1,"",warning_flag
        return line[6:].rstrip(), 1, "",warning_flag

    line=line[6:]

    ## if we are a continuation line: if long_line<>"" add to it.
    ## if next_line<>"" make it long_line and add to that instead.
    ## return with next_line=""

    if contline:
        # try to remove trailing comments with the "' rule
        n = count_quotes(long_line)
        line = strip_trailing_comment(line, n)
        
        if long_line <> "":
            long_line += " " + line.strip()
        if next_line <> "":
            long_line = next_line + " " + line.strip()
        return long_line, 1, "",warning_flag

    ## if we are a new line: if long line<>"", send it to splitter
    ## and return next_line=line
    ## if next_line<>"", set long_line=next_line, next_line=line
    ## and return to splitter

    line = strip_trailing_comment(line)
    
    if long_line<>"":
        next_line=line.rstrip()
        return long_line, 0, next_line,warning_flag
    if next_line<>"":
        long_line=next_line
        next_line=line.rstrip()
        return long_line, 0, next_line,warning_flag
    
    print "join_and_strip_fixed error. you should not see this message."
    sys.exit(1)


def add2dict(dict, name, item):
    """    
    checks whether name is in dict before adding the item.
    file is a string containing a file name,
    """
    if debug and dict.has_key(name):
        print "\nWarning: the key ", name, "already in dict will be overwritten"
	print "old value=",dict[name]
        print "new value=",item    
    dict[name] = item
    return 1, name


def add2tree(tree, currentname, newname):
    if tree.has_key(currentname):
        for name in newname:
            if tree[currentname][1].count(name)==0:
                tree[currentname][1].append(name)
    else:
        print "error libmkdep.add2tree: tree does not have currentname=",currentname
        print tree
        sys.exit(2)



def register_main_program(filename):
    """
    the main program must be put into a separate group
    """
    
    ## set group of current file to MAIN in globals.clean_source_files2

    if globals.clean_source_files2[filename][2] <> "MAIN":
	t = globals.clean_source_files2[filename]
	globals.clean_source_files2[filename] = (t[0],t[1],"MAIN",t[3],t[4],t[5],t[6])

    ## move it from whatever group to MAIN in globals.object_files
    
    for list in globals.object_files.keys():
	try:
	    if list<>'MAIN':
		print "\n",filename,"possible main program\n"
		globals.object_files[list].remove(filename)
	except ValueError:
	    pass
    if globals.object_files.has_key('MAIN'):
	try:
	    i = globals.object_files['MAIN'].index(filename)
	except ValueError:
	    globals.object_files['MAIN'].append(filename)
    else:
	globals.object_files['MAIN'] = [filename]


        
def write_header(file):
    file.write(\
"""
<html>
<head>
   <script language="javascript" src=".mkdep/mktree.js"></script>
   <link rel="stylesheet" href=".mkdep/mktree.css">
</head>
<body>

<table width=100% height=100%>
<tr>
<td width=600 height=100%>
<iframe width="100%" height=100% name="code" src="" frameborder=0>no iframe support in your browser</iframe>
</td>
<td width=300 style="vertical-align: top;">
<h1>
"""
+ calltree_title + \
"""
</h1>
mkdep was invoked like this:
<pre>
"""
+ string.join(sys.argv) +
"""
</pre>
* If a routine is called several times in a name space, only the first call is recorded.<br>
* Functions that are not declared will not show up in the call tree.<br>

<input type='button' onclick="expandToItem('tree1', document.getElementById('myText').value); return false;" value='Search for ROUTINE' />
<input type='text' id='myText' /><br>
<input type='button' onclick="expandTree('tree1'); return false;" value="expand all"/>
<input type='button' onclick="collapseTree('tree1'); return false;" value="collapse all"/>

<p>Format: ROUTINE (PATH TO CONTAINING FILE)

<ul class="mktree" id="tree1">
""")


def print_names(treefile, html, names, dict, level, stack=[]):
    """
    recursive function to write call trees.
    treefile is a file handler,
    html is the file handler for the html tree file,
    names is a list of names in the dictionary dict that we want to print trees for,
    dict contains tuples of (filename, [a list of calls], hosts),
    level is the current indent level.
    """

    htmllevel = "" # or level.
    
    if level == "":
        write_header(html)
    else:
        if len(names) > 0:
            html.write((htmllevel+'<ul>\n'))

    level += "    "
    
    #if root_path <> "":
    #    if root_path[-1] not in ["/","\\"]:
    #        root_path += "/"
    
    for name in names:
        file= ""
        if dict.has_key(name):
            file = dict[name][0]
        treefile.write((level + name + "   (" + file + ")\n" ))

        if file:
            html.write((htmllevel+'<li id="'+name+'">'+name+' <a target="code" href="'+ root_path + file+'">('+file+')</a> \n'))
        else:
            html.write((htmllevel + '<li id="'+name+'">' + name + '\n'))
        
        # detect recursion
        if name in stack:
            treefile.write((level + "Recursion?" + "\n"))
            html.write(('\n' + htmllevel + '</ul>\n'))
            return
        
        stack.append(name)
        
        if dict.has_key(name):
            print_names(treefile, html, dict[name][1], dict, level, stack)

        last = stack.pop()
        
    if len(names) > 0:
        html.write(('\n' + htmllevel + '</ul>\n'))

        
def write_call_trees():
    """
    Write the call tree('s if several main programs) to a separate file
    """
    if globals.calltree:
        trefile = open('.mkdep_trees', 'w')
        htmlfile = open((calltree_title + '.html'),'w')

        level=""   # the indentation level that increases as we go deeper into a tree
        stack=[]   # a call stack used to detect recursion

        # by default print call tree for each of the main programs
        if len(globals.programs.keys()) > 0 and treetop=="":
            print_names(trefile, htmlfile, globals.programs.keys(), globals.routines, level, stack)
        else:
            # if no programs are present, do it for each routine
            if treetop=="":
                print_names(trefile, htmlfile, globals.routines.keys(), globals.routines, level, stack)
            #unless the top of the tree is given on the command line
            else:
                print_names(trefile, htmlfile, [treetop], globals.routines, level, stack)
        trefile.close()

        htmlfile.write("</ul>\n</td>\n</tr>\n</table>\n</body>\n</html>\n")
        htmlfile.close()

        
def write_makefile(suffixes, use_mklib, f95, cc):
    """
    conditionally construct a makefile as a string and dump it.
    """
    
    makefile = '## Use GNU Make (gmake) to build\n\n\
PROJECT = proj1.exe\n\n\
include .mkdep_includes\n\
include .mkdep_objects\n\n\
F95=' + f95 + '\n\
BASEFLAGS=-g\n'

    if 'f90' in suffixes:  makefile += 'FREEFLAGS=$(INC) $(BASEFLAGS)\n'
    if 'f'   in suffixes:  makefile += 'FIXEDFLAGS=$(INC) $(BASEFLAGS)\n'
    if 'F90' in suffixes:  makefile += 'PFREEFLAGS=$(INC) $(BASEFLAGS)\n'
    if 'F'   in suffixes:  makefile += 'PFIXEDFLAGS=$(INC) $(BASEFLAGS)\n'
    if 'c'   in suffixes:  makefile += 'CFLAGS=$(INC) $(BASEFLAGS)\n'

    makefile += 'LFLAGS=\n'

    if 'c' in suffixes:
        makefile += 'CC=' + cc + '\n\
CFLAGS=\n'

    if use_mklib:
        makefile += 'AR=ar\n\
ARADDFLAG=r\n\
ARDELFLAG=d\n\n\
MAINLIB = lib$(PROJECT).a\n'

    makefile += '\n$(PROJECT):$(OBJS)\n'

    if use_mklib:
        makefile += '\tmklib -n 200 $(MAINLIB) .mkdep_objects\n\
\t$(AR) $(ARDELFLAG) $(MAINLIB) $(MAIN)\n\
\t$(F95) -o $(PROJECT) $(MAIN) $(MAINLIB) $(LFLAGS)\n\n'

    else:
        makefile += '\t$(F95) -o $(PROJECT) $(OBJS) $(LFLAGS)\n\n'

    makefile += 'include .mkdep_dependencies\n\n'

    if 'f90' in suffixes: makefile += '%.o : %.f90 ; $(F95) -c $(FREEFLAGS)   -o $@ $<\n'
    if 'f'   in suffixes: makefile += '%.o : %.f   ; $(F95) -c $(FIXEDFLAGS)  -o $@ $<\n'
    if 'F90' in suffixes: makefile += '%.o : %.F90 ; $(F95) -c $(PFREEFLAGS)  -o $@ $<\n'
    if 'F'   in suffixes: makefile += '%.o : %.F   ; $(F95) -c $(PFIXEDFLAGS) -o $@ $<\n'
    if 'c'   in suffixes: makefile += '%.o : %.c   ; $(CC)  -c $(CFLAGS)      -o $@ $<\n'

    makefile += '\ndep:\n\
\t' + string.join(sys.argv) + '\n\
clean:\n\
\trm -f *.mod *~ $(PROJECT) $(MAINLIB)\n\
\trmobjs .mkdep_objects\n'


    if os.access('Makefile', os.F_OK):
        #f=open('Makefile.tmp','w')
        #this can be annoying.
        #print "\nA template Makefile will be written to Makefile.tmp"
        print "found old Makefile, will not overwrite"
    else:
        f=open('Makefile','w')
        print "\nA template Makefile has been created"
        f.write(makefile)
        f.close()



def write_object_lists(object_files):
    """
    given a dictionary of object files
    replace the source suffixes with a .o, and write the
    objects to a new list to be included by our simple Makefile
    """
    
    print "\nIn the Makefile, the full list of object files is available as $(OBJS)."
    if len(object_files) > 1:
        print "This list of objects is split into subgroups, accessible as the variables;"

    ofile=open('.mkdep_objects','w')
    ostr="\n\nOBJS="

    for key in object_files:
        if len(object_files)>1:        
            print "$("+key+"),",
        ostr += "$("+key+") "
        ofile.write(('\n\n'+key+'='))
        for file in object_files[key]:
            dotpos = file.rfind(".")
            seppos = file.rfind("/")
            if globals.use_build_dir:
                objf = "build/" + file[seppos+1:dotpos]+".o"
            else:
                objf = file[:dotpos]+".o"
            ofile.write((" \\\n"+objf))

    if len(object_files)>1:
        print "which in the Makefile\ncan be utilized to compile groups of files with separate flags."

    ofile.write((ostr+"\n"))
    ofile.close()


def write_modmap():
    """
    writes the list of modules and their source file to .mkdep_modmap
    """
    mfile = open('.mkdep_modmap','w')
    mfile.write("This is a list over all modules defined in the given list \
of source files\nas well as the full path to the file with the \
definition used\n\n")

    for key in globals.module_map.keys():
	mfile.write((key.ljust(33)+" " + globals.module_map[key][0].strip()+"\n"))
    mfile.close()


def write_dependencies():
    """
    pretty print dependencies to file
    """
    dfile=open('.mkdep_dependencies','w')
    dfile.write('# This is the file with dependencies to be included in a Makefile\n')
    dfile.write('# If there is a change in a file on the right hand side of a colon, it will\n')
    dfile.write('# cause the file on the left hand side of the colon to be rebuilt.\n\n')


    # create a list of object files from the tuple object_files
    src=[]
    objs=[]
    for key in globals.object_files.keys():
        src += globals.object_files[key]
    for i in range(len(src)):
        objs.append(string.replace(src[i],".f90",".o"))
        
    for item in globals.clean_source_files2.items():
	(filename, tuple) = item
	deps = tuple[4]
	is_include=tuple[1]

	if len(deps)>0:
	    objf = filename
            seppos = filename.rfind("/")
	    dotpos = filename.rfind(".")
            
            ## first write the object file name on the left side of :
	    if is_include==0:
                if globals.use_build_dir:
                    objf = "build/" + filename[seppos+1:dotpos]+".o"
                else:
                    objf = filename[:dotpos]+".o"                    
	    str = objf + ":"

            ## then add dependencies. break lines if they are longer than 80 chars
	    for ind in range(len(deps)):
                ## if we use a build dir, all .o files are there.
                if globals.use_build_dir:
                    if deps[ind] in objs:
                        seppos = deps[ind].rfind("/")
                        str += " " + "build/" + deps[ind][seppos+1:]
                    else:
                        str += " " + deps[ind]
                else:
                    str += " " + deps[ind]
		if len(str)>80 and ind<len(deps)-1:
		    str += " \\\n  "
		    dfile.write(str)
		    str=""
		    continue
	    dfile.write((str+"\n"))

    dfile.close()


def save_maps():
    """
    dump the maps to file
    """    
    mapfile=open('.mkdep_restart_file','w')

    # mark all files for deletion
    for (key,(type, is_include, do_opt, mtime, deps, reparse, deleted)) in globals.clean_source_files2.items():
        globals.clean_source_files2[key]=(type, is_include, do_opt, mtime, deps, reparse, 1)
    
    pickle.dump(globals.clean_source_files2, mapfile)
    pickle.dump(globals.object_files, mapfile)
    pickle.dump(globals.incdirs, mapfile)
    pickle.dump(globals.module_map, mapfile)
    pickle.dump(globals.routines, mapfile)
    pickle.dump(globals.functions, mapfile)
    pickle.dump(globals.programs, mapfile)
    pickle.dump(globals.suffixes, mapfile)
    pickle.dump(globals.calltree, mapfile)
    mapfile.close()
       

def load_maps():
    if os.access('.mkdep_restart_file', os.F_OK):
        mapfile=open('.mkdep_restart_file','r')
        try:
            globals.clean_source_files2 = pickle.load(mapfile)
            globals.object_files= pickle.load(mapfile)
            globals.incdirs = pickle.load(mapfile)
            globals.module_map = pickle.load(mapfile)
            globals.routines = pickle.load(mapfile)
            globals.functions = pickle.load(mapfile)
            globals.programs = pickle.load(mapfile)
	    globals.suffixes = pickle.load(mapfile)
            globals.prevcalltree = pickle.load(mapfile)
        except EOFError:
            print "pickle file too short, full rescan necessary"
            globals.object_files={}
            globals.clean_source_files2={}
            globals.module_map = {}
            globals.routines={}
            globals.functions={}
            globals.programs={}
            globals.incdirs = []
            globals.suffixes = []
        mapfile.close()

        
def install_files():
    """
    copy css and icons for the call tree to the .mkdep directory
    """
    if not os.access('.mkdep', os.F_OK):
        mkdeppath = string.join(sys.argv[0].split(os.sep)[:-1], os.sep)
        os.mkdir('.mkdep')
        import shutil
        for f in ['bullet.gif','minus.gif','mktree.css','mktree.js','plus.gif']:
            shutil.copy((mkdeppath+os.sep+'.mkdep'+os.sep+f), '.mkdep')


def find_uncalled():
    """
    go through globals.routines and look for names that are never referenced
    """
    # for each routine
    for key in globals.routines.keys():
        # loop over all the routines that are called
        for name in globals.routines[key][1]:
            # and add the callee to the called-from list in the dict
            if globals.routines.has_key(name):
                globals.routines[name][2].append(key)
            else:
                # probably no source file for this routine
                pass

    uncalled=[]
    for key in globals.routines.keys():
        if len(globals.routines[key][2]) == 0:
            uncalled.append((key, globals.routines[key][0] ))

    if len(uncalled)>0:        
        print "\n\nSeveral routines are never called by the set of files in the project file"
        print "(see .mkdep_uncalled). If the set of files in the project is a library,"
        print "you may remove only those files that are not a part of the interface"
        
        uncfile=open('.mkdep_uncalled','w')
        uncfile.write("These routines are never called by the set of files in the project file\n\n")
        
        for key,file in uncalled:
            uncfile.write( (key+" ("+file+")\n"))
        uncfile.close()
            
