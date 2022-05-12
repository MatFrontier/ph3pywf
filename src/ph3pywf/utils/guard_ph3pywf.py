#!/usr/bin/python

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

from ast import Raise
from logging import raiseExceptions
import sys
import time
from datetime import datetime

# import subprocess
from ph3pywf.utils.mission_control import check_progress_and_rerun, get_tag_from_fw_id

SLEEP_TIME = 600  # TODO: allow user to specify this via argv


def main():
    if len(sys.argv) == 1:
        print("tag or fw_id not specified")
        raise
    elif sys.argv[1].isnumeric():
        tag = get_tag_from_fw_id(int(sys.argv[1]))
    else:
        tag = sys.argv[1]

    user = "jerrylai"
    # if len(sys.argv) == 3:
    #     user = sys.argv[2]

    while True:
        print(datetime.utcnow())
        print("tag = {}".format(tag))
        check_progress_and_rerun(tag)
        print(f"sleep for {SLEEP_TIME} seconds")
        time.sleep(SLEEP_TIME)


if __name__ == "__main__":
    main()