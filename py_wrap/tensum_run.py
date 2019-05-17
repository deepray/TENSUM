#!/usr/bin/env python
# Python wrapper for tensum
# Author: Deep Ray, TIFR CAM, Bangalore
# Date  : 14 December, 2015

import os
import numpy
from sys import argv, exit
import re
import os, errno

def clean(string):
    sa, sb = '', ''
    a = string.split(' ')
    for s in a:
        sa = sa + s
    if "\t" in sa:
        sb = sa.strip("\t")
        for s in sb:
            sa = sa + s
    return sa

def extract_info(s):
    if len(s)==2:
        key = s[0]
        value = s[1]
        return str(clean(key)), str(value)
    else:
        return "error", "error"

def rs(string):
    ns = ''
    for i in string:
        if i not in [' ']:
            ns = ns + i
    return ns


def read_parameter(input):
    in_file = open(input, 'r')
    d = {}
    for line in in_file:
        s =  line[:-1].split("=")
        key, value = extract_info(s)
        if key != "error":
            d[key] = value
    in_file.close()
    return d

def gen_and_part_mesh(d):
    fname = rs(d['mesh_file_name'])
    nparts = int(rs(d['mesh_parts']))
    dim = rs(d['mesh_dimension'])
    part_dir_loc = rs(d['partition_dir_loc'])
    grid_gen_log = "grid_gen.log"
    file_head = fname[:-4]
    file_ext  = fname[-4:]
    if(file_ext == '.geo'):
        mesh_file = file_head+'.msh'
    elif(file_ext == '.msh'):
        mesh_file = fname
    else:
        exit(' Error: Unknown format of mesh-file "'+fname+'". Only ".geo" or ".msh"'+\
              ' formats allowed.\n')        
    if(nparts > 1):
        comline = "gmsh -"+ dim  + " -part " + str(nparts) + " " + fname + " > " + grid_gen_log
    else:
        comline = "gmsh -"+ dim + " " + fname + " > " + grid_gen_log    
    print ("  Generating mesh using gmsh from ", fname)
    print ("  --- log output available in", grid_gen_log)
    os.system(comline)
    
    comline = "grid_part -D " + dim + " -I " + mesh_file+ " -P " + nparts+\
              " -L " + part_dir_loc 
    #print "  --- Partitioning mesh into ", nparts, " parts\n\n"
    os.system(comline)

def init_solver(d):
    pname = rs(d['solver_param_file'])
    nparts = rs(d['mesh_parts'])
    nprocs = rs(d['nprocs'])
    dim = rs(d['mesh_dimension'])
    part_dir_loc = rs(d['partition_dir_loc'])
    if(dim == '2'):
       solver_exe = "tensum"
    elif(dim == '3'):  
       exit('  Error: 3D solver is still under developement.') 
       solver_exe = ""   
    else:
       exit('  Error: unknown dim = '+dim)
          
    p_flag = ""
    r_flag = ""
    d_flag = ""
    b_flag = ""
    v_flag = ""
    if(rs(d['use_solver']) == 'no'):
       p_flag = "-p"
    if(rs(d['read_from_restart']) == 'yes'):
       r_flag = "-r"
    if(rs(d['print_cells']) == 'yes'):
       d_flag = "-d"
    if(rs(d['print_bounds']) == 'yes'):
       b_flag = "-b"
    if(rs(d['verbose']) == 'yes'):
       v_flag = "-v"
    
    #print "  Starting simulator  "
    comline = "mpirun -n " + nprocs + " " + solver_exe \
              + " -npart " + nparts + " " \
              + " -L " + part_dir_loc + " " \
              + p_flag + " "\
              + r_flag + " "\
              + b_flag + " "\
              + v_flag + " "\
              + d_flag \
              + " -i " + pname
    #print(comline)
    os.system(comline)

# wrapper begins here:
script, input = argv
param_map = read_parameter(input)

print ("---------------------------------------------------------------------------\n")
#print " Initiating requested tasks: \n"

# execution of main operations
if rs(param_map['gen_mesh_and_part']) == 'yes':
    gen_and_part_mesh(param_map)
if rs(param_map['initiate_solver']) == 'yes':
    init_solver(param_map)

