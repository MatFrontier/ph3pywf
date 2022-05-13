#!/usr/bin/python

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

from ast import Raise
from logging import raiseExceptions
import sys
import time
from datetime import datetime
import argparse

# import subprocess
from ph3pywf.utils.mission_control import check_progress_and_rerun, get_tag_from_fw_id

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--fw_ids",
        nargs="*",
        type=int,
        default=[],
    )
    parser.add_argument(
        "--tags",
        nargs="*",
        type=str,
        default=[],
    )
    parser.add_argument(
        "--sleep",
        nargs=1,
        type=int,
        default=600,
    )

    # if len(sys.argv) == 1:
    #     print("tag or fw_id not specified")
    #     raise
    # elif sys.argv[1].isnumeric():
    #     tag = get_tag_from_fw_id(int(sys.argv[1]))
    # else:
    #     tag = sys.argv[1]

    args = parser.parse_args()
    SLEEP_TIME = args.sleep
    tags = args.tags
    
    for fw_id in args.fw_ids:
        tag = get_tag_from_fw_id(fw_id)
        if tag not in tags:
            tags.append(tag)

    while True:
        print(datetime.utcnow())
        for tag in tags:
            print("tag = {}".format(tag))
            check_progress_and_rerun(tag)
        print(f"sleep for {SLEEP_TIME} seconds")
        time.sleep(SLEEP_TIME)


if __name__ == "__main__":
    main()