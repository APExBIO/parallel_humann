#!/home/shaojun/Tools/miniconda3/envs/humann3/bin/python

# -*- coding: utf-8 -*-
import re
import sys

from humann.humann import main

if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw?|\.exe)?$', '', sys.argv[0])
    sys.exit(main())