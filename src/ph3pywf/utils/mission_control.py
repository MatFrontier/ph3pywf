# coding: utf-8

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

from fireworks import LaunchPad

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