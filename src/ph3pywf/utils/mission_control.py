# coding: utf-8

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

from fireworks import LaunchPad
from fireworks.fw_config import LAUNCHPAD_LOC
from atomate.vasp.database import VaspCalcDb
import os

def check_progress_and_rerun(tag):
    
    launchpad = LaunchPad.auto_load()
    
    # Check progress 
    RUNNING_ids = launchpad.get_fw_ids(
        {"state": "RUNNING", "name": {"$regex": tag}})
    COMPLETED_ids = launchpad.get_fw_ids(
        {"state": "COMPLETED", "name": {"$regex": tag}})
    FIZZLED_ids = launchpad.get_fw_ids(
        {"state": "FIZZLED", "name": {"$regex": tag}})
    ALL_ids = launchpad.get_fw_ids({"name": {"$regex": tag}})

    print(f"total: {len(ALL_ids)}\nRUNNING: {len(RUNNING_ids)}\nFIZZLED: {len(FIZZLED_ids)}\nCOMPLETED: {len(COMPLETED_ids)}")

    # Rerun FIZZLED
    for fw_id in FIZZLED_ids:
        print(f"INFO: rerunning FIZZLED fw: {fw_id}")
        launchpad.rerun_fw(fw_id)

    # Rerun lost runs
    LOST_ids = launchpad.detect_lostruns(expiration_secs=1*60*60,
                                         query={"name": {"$regex": tag}})[1]
    for fw_id in LOST_ids:
        print(f"INFO: rerunning lost run: {fw_id}")
        launchpad.rerun_fw(fw_id)

    print("DONE")

def check_time_of_each_run(tag):
    db_file = os.path.join(os.path.dirname(LAUNCHPAD_LOC), "db.json")
    
    mmdb = VaspCalcDb.from_db_file(db_file)
    docs_p = mmdb.collection.find(
        {
            "task_label": {"$regex": f"{tag}"},
        }
    )

    docs = []
    for p in docs_p:
        docs.append(p)
    # print(type(docs[0]))
    for d in docs:
        print("{{task_label:\"{}\"}}".format(d["task_label"]))
        if "nsites" in d:
            print("nsites = {}".format(d["nsites"]))
        if "input" in d:
            print("EDIFF = {}".format(d["input"]["incar"]["EDIFF"]))
        if "run_stats" in d:
            if "standard" in d["run_stats"]:
                print("Maximum memory used (kb) = {}".format(d["run_stats"]["standard"]["Maximum memory used (kb)"]))
            print("Elapsed time (sec) = {}".format(d["run_stats"]["overall"]["Elapsed time (sec)"]))
        print("")