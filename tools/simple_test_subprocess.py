#

import sys
import os
import re
import glob
import subprocess

################

x=[]
cmdstr = ['ls /Users/codar/', 'ls -l /Users/Codar', 'ls -la /Users/Codar']
x.extend(cmdstr) # extending a list 
x


# careful your extending a string and not a list anymore
x=[]
x.extend(cmdstr[2]) # extending a string type
x



x=[]
x.append(cmdstr[2])
x
# Out[10]: ['ls -la /Users/Codar']

