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
        help="fw_id(s) of any FW in workflow(s)",
        type=int,
        default=[],
    )
    parser.add_argument(
        "--tags",
        nargs="*",
        help="tag(s) or task_label(s) of workflow(s)",
        type=str,
        default=[],
    )
    parser.add_argument(
        "-s",
        "--sleep",
        nargs="?",
        help="sleep time in seconds",
        type=int,
        default=600,
    )

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