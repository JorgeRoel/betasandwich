import argparse
import datetime
import os
from logger import LoggingManager
from loop_closure import close_hairpins, close_arches
from protocols import sandwich
from sequence_design import prepare_xml, design_sequence
from util import tidy_pdb, split

from pyrosetta import *

pyrosetta.init()


def decompose_code(code):

    characters = [char for char in code]

    return int(characters[0]), int(characters[1]), int(characters[2]), int(characters[3])


def parse_command_line():

    parser = argparse.ArgumentParser(prog="beta_sandwich")

    parser.add_argument(
        "--code", 
        "-code",
        "-c",
        help="Code for generating sandwich topology", 
        dest="code"
    )

    parser.add_argument(
        "--rotation", 
        "-rotation",
        "-r",
        help="Options are + or - for either positive (left) or negative (right) rotational sampling", 
        dest="rotation"
    )

    parser.add_argument(
        "--suffix", 
        "-suffix",
        "-s",
        help="Suffix for output files", 
        dest="suffix"
    )

    return parser.parse_args()


if __name__ == '__main__':

    args = parse_command_line()

    LIBRARY_PATH = os.path.join(os.path.realpath(__file__).replace("beta_sandwich.py",""), "library")

    code = int(args.code)
    rotation = args.rotation
    suffix = args.suffix

    nres1, nres2, ori1, ori2 = decompose_code(args.code)
    beta1 = f"beta_sheet_{nres1}_{ori1}.pdb"
    beta2 = f"beta_sheet_{nres2}_{ori2}.pdb"

    log_file = f"sandwich_{suffix}.log"
    pdb_bb = f"backbone_{suffix}.pdb"
    pdb_cls = f"sandwich_{suffix}.pdb"
    xml_file = f"sandwich_{suffix}.xml"

    log = LoggingManager.get_logger("beta_sandwich", file_name=log_file)

    log.info(f"STARTING BUILDER: {pdb_cls} ({datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')})")

    sheet1 = pose_from_pdb(os.path.join(LIBRARY_PATH, "4strands", beta1))
    sheet2 = pose_from_pdb(os.path.join(LIBRARY_PATH, "3strands", beta2))

    backbone = sandwich(sheet1, sheet2, nres1, nres2, ori1, ori2, rotation, log_file)
    
    backbone.dump_pdb(pdb_bb)
    
    new_numbering = tidy_pdb(pdb_bb, nres1, nres2, ori1, ori2)

    with open(pdb_bb, "w") as renum_pdb:
        for line in new_numbering:
            renum_pdb.write(line)
    
    pose, bp_backbone, hairpins = close_hairpins(pdb_bb, nres1, nres2, code, suffix)

    try:
        pose, nres_total, arches = close_arches(pose, nres1, nres2, bp_backbone, suffix)

        if pose.size() == nres_total:
            pose.dump_pdb(pdb_cls)
            os.remove(pdb_bb)
            log.info(f"COMPLETED BUILDER: {pdb_cls} ({datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')})")
            log.info(f"STARTING DESIGNER: {pdb_cls} ({datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')})")
            prepare_xml(nres1, nres2, hairpins, arches, nres_total, suffix)
            
            for sequence in range(5):
                pose_dsg = pose.clone()
                pdb_seq = f"sandwich_{suffix}_{sequence+1}_seq.pdb"
                report_file = f"sandwich_{suffix}_{sequence+1}.out"
                design_pose = design_sequence(pose_dsg, report_file, xml_file, nres_total)
                design_pose.dump_pdb(pdb_seq)
                log.info(f"COMPLETED DESIGNER: {pdb_seq} ({datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')})")

    except:
        log.info(f"FAILED: {pdb_cls} ({datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')})")
