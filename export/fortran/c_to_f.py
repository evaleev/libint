#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Convert selected C structs to Fortran derived types.

Usage:
C input to stdin, Fortran output to stdout.
Specify names of C structs to be converted as command line arguments.
'''

import re, sys

c_struct_names = sys.argv[1:]
c_in = sys.stdin
f_out = sys.stdout

def extract_c_struct_fields(c_code, struct_name):
    struct_re = re.compile("^\s*typedef\s*struct\s*\{([^\}]*)\}\s*"+struct_name+"\s*;",re.MULTILINE | re.DOTALL)
    m = struct_re.search(c_code)
    if not m: return []
    fields = m.group(1).splitlines()
    fields = filter(lambda k: (not k.startswith('#')) and k.strip(), fields)
    fields = [_.strip() for _ in ' '.join(fields).split(';') if _]
    return fields

def parse_c_struct_fields(fields):
    dtype = []
    size = []
    name = []
    for f in fields:
        words = f.split()
        var = words.pop()
        m = re.match("(\w*)(?:\[(\w+)\])?", var)
        name.append(m.group(1))
        size.append(m.group(2))
        dtype.append(' '.join(words))

    return (dtype, size, name)

def convert_name_c_to_f(name):
    if name.startswith('_'):
        name = 'f' + name
    return name

def convert_type_c_to_f(dtype):
    dtype = re.sub("(const|volatile|mutable)", "", dtype)
    dtype = dtype.strip()
    if dtype.endswith("*"):
        return "type(c_ptr)"
    elif dtype == "float":
        return "real(c_float)"
    elif dtype == "double":
        return "real(c_double)"
    elif dtype == "int":
        return "integer(c_int)"

instring = c_in.read()

for name in c_struct_names:
    fields = extract_c_struct_fields(instring, name)
    if not fields: continue
    dtype, size, variable = parse_c_struct_fields(fields)
    variable = [convert_name_c_to_f(_) for _ in variable]
    dtype = [convert_type_c_to_f(_) for _ in dtype]

    f_out.write("type, bind(c) :: {}\n".format(name))
    for t, s, v in zip(dtype,size,variable):
        if s:
            f_out.write("   {}, DIMENSION({}) :: {}\n".format(t, s, v))
        else:
            f_out.write("   {} :: {}\n".format(t, v))
    f_out.write("end type\n")

