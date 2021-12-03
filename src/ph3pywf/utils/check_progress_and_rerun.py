#!/usr/bin/python

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

import sys
import time
from datetime import datetime
import subprocess
from ph3pywf.utils.mission_control import check_progress_and_rerun

def main():
    tag = sys.argv[1] # TODO: change this to use fw_id in the future
    user = "jerrylai"
    # if len(sys.argv) == 3:
    #     user = sys.argv[2]
    
    sleep_time = 1200
    
    while True:
        print(datetime.utcnow())
        
        sq = subprocess.run(["squeue","-u",user], capture_output=True, text=True).stdout
        n_jobs = sq.count("\n") - 1
        print(f"{n_jobs} jobs running")
        if n_jobs == 0:
            print("no slurm job running, exiting")
            break
        
        check_progress_and_rerun(tag)
        
        print(f"sleep for {sleep_time} seconds")
        time.sleep(sleep_time)
    

if __name__ == "__main__":
    main()