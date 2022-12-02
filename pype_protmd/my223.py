import os
import sys
from lib2to3.main import main
if __name__=='__main__':
    sys.path.append(os.path.dirname(os.path.realpath(__file__)))
    sys.exit(main('py223_fixers'))
    
