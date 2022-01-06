#!/usr/bin/python

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

import sys
# import time
from datetime import datetime
# import subprocess
from ph3pywf.utils.mission_control import check_progress_and_rerun

def main():
    tag = sys.argv[1] # TODO: change this to use fw_id in the future
    user = "jerrylai"
    # if len(sys.argv) == 3:
    #     user = sys.argv[2]
    
    check_progress_and_rerun(tag)
    

if __name__ == "__main__":
    main()