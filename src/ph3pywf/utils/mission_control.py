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
    RUNNING_ids = launchpad.get_fw_ids({"state": "RUNNING", "name": {"$regex": tag}})
    COMPLETED_ids = launchpad.get_fw_ids(
        {"state": "COMPLETED", "name": {"$regex": tag}}
    )
    FIZZLED_ids = launchpad.get_fw_ids({"state": "FIZZLED", "name": {"$regex": tag}})
    ALL_ids = launchpad.get_fw_ids({"name": {"$regex": tag}})

    print(
        f"total: {len(ALL_ids)}\nRUNNING: {len(RUNNING_ids)}\nFIZZLED: {len(FIZZLED_ids)}\nCOMPLETED: {len(COMPLETED_ids)}"
    )

    # Rerun FIZZLED
    for fw_id in FIZZLED_ids:
        d = launchpad.get_fw_dict_by_id(fw_id)
        n_archived_launches = len(d["archived_launches"])

        # print a warning if a FW has failed more than 1 time
        if n_archived_launches > 1:
            print(f"WARNING: fw: {fw_id}, has FIZZLED {n_archived_launches} times")

        print(f"INFO: rerunning FIZZLED fw: {fw_id}")
        launchpad.rerun_fw(fw_id)

    # Rerun lost runs
    LOST_ids = launchpad.detect_lostruns(
        expiration_secs=1 * 60 * 60, query={"name": {"$regex": tag}}
    )[1]
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
        print('{{task_label:"{}"}}'.format(d["task_label"]))
        if "nsites" in d:
            print("nsites = {}".format(d["nsites"]))
        if "input" in d:
            print("EDIFF = {}".format(d["input"]["incar"]["EDIFF"]))
        if "run_stats" in d:
            if "standard" in d["run_stats"]:
                print(
                    "Maximum memory used (kb) = {}".format(
                        d["run_stats"]["standard"]["Maximum memory used (kb)"]
                    )
                )
            print(
                "Elapsed time (sec) = {}".format(
                    d["run_stats"]["overall"]["Elapsed time (sec)"]
                )
            )
        print("")


def get_tag_from_fw_id(fw_id):
    launchpad = LaunchPad.auto_load()
    wf_dict = launchpad.workflows.find_one({"nodes": fw_id})

    return wf_dict["metadata"]["label"]


def get_resource_report(tag):
    db_file = os.path.join(os.path.dirname(LAUNCHPAD_LOC), "db.json")

    mmdb = VaspCalcDb.from_db_file(db_file)
    docs_p = mmdb.collection.find(
        {
            "task_label": {"$regex": f"{tag}"},
        }
    )

    docs = []
    print("Fetching tasks...")
    for p in docs_p:
        docs.append(p)
    # print(type(docs[0]))

    n_tasks = 0
    total_time = 0
    total_mem_used = 0
    cores_used = {}
    for d in docs:
        if "run_stats" in d:
            if "standard" in d["run_stats"]:
                total_mem_used += d["run_stats"]["standard"]["Maximum memory used (kb)"]
                cores = d["run_stats"]["standard"]["cores"]
                cores_used[cores] = cores_used[cores] + 1 if cores in cores_used else 0
            total_time += d["run_stats"]["overall"]["Elapsed time (sec)"]
        n_tasks += 1

    print("Total number of tasks = {}".format(n_tasks))
    print("Total computation time (h) = {}".format(total_time / 3600))
    print("Average time per task (s) = {}".format(total_time / n_tasks))
    print("Average memory used (kb) = {}".format(total_mem_used / n_tasks))
    for cores, count in cores_used.items():
        print(f"{count} tasks used {cores} cores")