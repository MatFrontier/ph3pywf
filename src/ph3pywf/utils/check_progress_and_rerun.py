#!/usr/bin/python

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

import sys
import time
import subprocess
from ph3pywf.utils.mission_control import check_progress_and_rerun

def main():
    tag = sys.argv[1]
    user = "jerrylai"
    # if len(sys.argv) == 3:
    #     user = sys.argv[2]
        
    while True:
        sq = subprocess.run(["squeue","-u",user], capture_output=True, text=True).stdout
        n_jobs = sq.count("\n") - 1
        print(f"{n_jobs} jobs running")
        if n_jobs == 0:
            print("no slurm job running")
            break
        
        check_progress_and_rerun(tag)
        time.sleep(30)
    

if __name__ == "__main__":
    main()