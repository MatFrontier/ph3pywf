#!/usr/bin/python

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

import sys
import time
from ph3pywf.utils.mission_control import check_progress_and_rerun

def main():
    tag = sys.argv[1]
    while True: # need to add some conditions
        check_progress_and_rerun(tag)
        time.sleep(3600)
    

if __name__ == "__main__":
    main()