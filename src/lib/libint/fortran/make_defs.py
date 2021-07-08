import re, sys
f_in = open(sys.argv[1], 'r')
f_out = open(sys.argv[2], 'w')
for line in f_in:
    if re.match('^#', line) and not re.match('#include', line):
        f_out.write(line)
f_in.close()
f_out.close()
