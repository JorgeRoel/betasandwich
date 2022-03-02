import os
import random
from blueprint import *
from chiralities import codes, chiralities
from filters import check_short_hairpins, check_long_hairpins
from logger import LoggingManager
from protocols import blueprint_hairpins

from pyrosetta import *

pyrosetta.init()


def close_hairpins(initial_backbone, nres1, nres2, code, suffix):

    ref_pose = pose_from_pdb(initial_backbone)
    bp_hairpins = f"hp_{suffix}.bp"
    bp_backbone = f"bb_{suffix}.bp"
    log_file = f"sandwich_{suffix}.log"
    hairpins = []

    log = LoggingManager.get_logger("close_hairpins", file_name=log_file)

    log.info("HAIRPINS")

    chiral = chiralities[codes.index(code)]
    chiral1 = chiral[0]
    chiral2 = chiral[1]
    chiral3 = chiral[2]

    if chiral1 == "R":
        lh = 5
        abego_string = "BAAGB"
        pyrosetta.init("-vall filtered.vall.hp5")
    elif chiral1 == "L":
        lh = 2
        abego_string = None
        pyrosetta.init("-vall filtered.vall.hp2")
    
    blueprint_hairpins(nres1, nres2, lh, bp_hairpins, abego_string)

    bdr = pyrosetta.rosetta.protocols.fldsgn.BluePrintBDR()
    censcorefxn = pyrosetta.create_score_function("fldsgn_cen.wts")
    bdr.scorefunction(censcorefxn)
    bdr.ss_from_blueprint(True)

    bp = Blueprint(bp_hairpins)
    
    if lh == 5:
        bp.bp_data[nres1+2][1] = "G"
        bp.bp_data[nres1+3][1] = "G"
        bp.bp_data[nres1+4][1] = "G"
    
    bp.dump_blueprint(bp_hairpins)
    bdr.set_blueprint(bp_hairpins)
    
    for i in range(10):
        
        log.info(f"HAIRPIN_1: TRYING LENGTH {lh}")
        idx = nres2
        pose1 = ref_pose.clone()
        bdr.apply(pose1)

        try:
            if lh == 5:
                check, abego = check_long_hairpins(pose1, idx)
            elif lh == 2:
                check, abego = check_short_hairpins(pose1, idx)
        except:
            continue

        if check == True:
            idx = idx + lh
            nres_total = nres1*4 + nres2*3 + lh
            log.info(f"HAIRPIN_1: {i+1} ATTEMPTS WITH ABEGO {abego}")
            hairpins.append(lh)
            break

    if pose1.size() == nres_total:

        if chiral2 == "R":
            lh = 5
            abego_string = "BAAGB"
            pyrosetta.init("-vall filtered.vall.hp5")
        elif chiral2 == "L":
            lh = 2
            abego_string = None
            pyrosetta.init("-vall filtered.vall.hp2")

        idx = idx + nres2 + nres1
        bdr = pyrosetta.rosetta.protocols.fldsgn.BluePrintBDR()
        bdr.scorefunction(censcorefxn)
        bdr.ss_from_blueprint(True)
        
        bp = Blueprint(bp_hairpins)
        bp.add_loop_between_segments(idx=idx+1, length=lh, abego_string=abego_string)

        if lh == 5:
            bp.bp_data[idx+2][1] = "G"
            bp.bp_data[idx+3][1] = "G"
            bp.bp_data[idx+4][1] = "G"

        bp.dump_blueprint(bp_hairpins)
        bdr.set_blueprint(bp_hairpins)

        for i in range(10):

            log.info(f"HAIRPIN_2: TRYING LENGTH {lh}")
            pose2 = pose1.clone()
            bdr.apply(pose2)

            try:
                if lh == 5:
                    check, abego = check_long_hairpins(pose2, idx)
                elif lh == 2:
                    check, abego = check_short_hairpins(pose2, idx)
            except:
                continue

            if check == True:
                idx = idx + lh
                nres_total = nres_total + lh
                log.info(f"HAIRPIN_2: {i+1} ATTEMPTS WITH ABEGO {abego}")
                hairpins.append(lh)
                break

    if pose2.size() == nres_total:

        if chiral3 == "R":
            lh = 5
            abego_string = "BAAGB"
            pyrosetta.init("-vall filtered.vall.hp5")
        elif chiral3 == "L":
            lh = 2
            abego_string = None
            pyrosetta.init("-vall filtered.vall.hp2")

        idx = idx + nres1*2 + nres2
        bdr = pyrosetta.rosetta.protocols.fldsgn.BluePrintBDR()
        bdr.scorefunction(censcorefxn)
        bdr.ss_from_blueprint(True)
    
        bp = Blueprint(bp_hairpins)
        bp.add_loop_between_segments(idx=idx+1, length=lh, abego_string=abego_string)

        if lh == 5:
            bp.bp_data[idx+2][1] = "G"
            bp.bp_data[idx+3][1] = "G"
            bp.bp_data[idx+4][1] = "G"

        bp.dump_blueprint(bp_hairpins)
        bdr.set_blueprint(bp_hairpins)

        for i in range(10):

            log.info(f"HAIRPIN_3: TRYING LENGTH {lh}")
            pose3 = pose2.clone()
            bdr.apply(pose3)

            try:
                if lh == 5:
                    check, abego = check_long_hairpins(pose3, idx)
                elif lh == 2:
                    check, abego = check_short_hairpins(pose3, idx)
            except:
                continue

            if check == True:
                bp.reindex_blueprint(0)
                bp.freeze_all()
                bp.bp_data[0][-2] = "R"
                bp.bp_data[-1][-2] = "R"
                bp.dump_blueprint(bp_backbone)
                os.remove(bp_hairpins)
                log.info(f"HAIRPIN_3: {i+1} ATTEMPTS WITH ABEGO {abego}")
                log.info("HAIRPINS COMPLETED")
                hairpins.append(lh)
                return pose3, bp_backbone, hairpins
                break


def close_arches(pose, nres1, nres2, bp_backbone, suffix):

    pyrosetta.init()

    nres_total_ref = pose.size()
    bp_arch1 = f"ac1_{suffix}.bp"
    bp_arch2 = f"ac2_{suffix}.bp"
    bp_sandwich = f"sandwich_{suffix}.bp"
    log_file = f"sandwich_{suffix}.log"
    arches = []

    cbreak = pyrosetta.rosetta.protocols.simple_filters.ChainBreak()

    log = LoggingManager.get_logger("close_arches", file_name=log_file)

    log.info("ARCHES")

    for i in range(10):

        pose1 = pose.clone()

        # ARCH 5-6
        bp = Blueprint(bp_backbone)

        idx = bp.segment_dict["E3"].bp_data[0][0]+nres1+nres2
        l = random.sample([3, 4, 5, 6], 1)[0]

        log.info(f"ARCH_5-6: TRYING LENGTH {l}")

        bp.add_loop_between_segments(idx=idx, length=l)
        bp.dump_blueprint(bp_arch1)
        
        bdr = pyrosetta.rosetta.protocols.fldsgn.BluePrintBDR()
        censcorefxn = pyrosetta.create_score_function("fldsgn_cen.wts")
        bdr.scorefunction(censcorefxn)
        bdr.ss_from_blueprint(True)
        bdr.set_blueprint(bp_arch1)
        bdr.apply(pose1)

        if pose1.size() == nres_total_ref + l:
            
            nres_total = pose1.size()
            abego = list(rosetta.core.sequence.get_abego(pose1))[idx:idx+l]
            bp.reindex_blueprint(0)
            bp.freeze_all()
            bp.bp_data[0][-2] = "R"
            bp.bp_data[-1][-2] = "R"
            bp.dump_blueprint(bp_arch1)
            log.info(f"ARCH_5-6: {i+1} ATTEMPTS WITH LENGTH {l} AND ABEGO {abego}")
            arches.append(l)
            break

        elif i == 9:
            os.remove(bp_backbone)
            os.remove(bp_arch1)
            os.remove(bp_arch2)
            log.info("ARCHE5_6: FAILED")

    if nres_total == nres_total_ref + l:
    
        for i in range(10):

            pose2 = pose1.clone()

            # MODEL BOTH ARCHES (2-3 / 4-5) AT THE SAME TIME
            bp = Blueprint(bp_arch1)

            idx = bp.segment_dict["E2"].bp_data[0][0]+nres2
            other_idx = bp.segment_dict["E3"].bp_data[0][0]+nres1+l
            
            l = random.sample([3, 4, 5, 6], 1)[0]
            other_l = random.sample([3, 4, 5, 6], 1)[0]

            log.info(f"ARCH_2-3: TRYING LENGTH {l}")
            log.info(f"ARCH_4-5: TRYING LENGTH {other_l}")

            bp.add_loop_between_segments(idx=idx, length=l, other_idx=other_idx, other_length=other_l)
            bp.dump_blueprint(bp_arch2)

            bdr = pyrosetta.rosetta.protocols.fldsgn.BluePrintBDR()
            censcorefxn = pyrosetta.create_score_function("fldsgn_cen.wts")
            bdr.scorefunction(censcorefxn)
            bdr.ss_from_blueprint(True)
            bdr.set_blueprint(bp_arch2)
            bdr.apply(pose2)

            # Apply chain break filter
            break_filter = cbreak.report_sm(pose2)

            if pose2.size() == nres_total+l+other_l and break_filter == 0.0:

                abego = list(rosetta.core.sequence.get_abego(pose2))[idx:idx+l]
                other_abego = list(rosetta.core.sequence.get_abego(pose2))[other_idx:other_idx+other_l]
                bp.reindex_blueprint(0)
                bp.freeze_all()
                bp.bp_data[0][-2] = "R"
                bp.bp_data[-1][-2] = "R"
                bp.dump_blueprint(bp_sandwich)
                os.remove(bp_backbone)
                os.remove(bp_arch1)
                os.remove(bp_arch2)
                log.info(f"ARCH_2-3: {i+1} ATTEMPTS WITH LENGTH {l} AND ABEGO {abego}")
                log.info(f"ARCH_4-5: {i+1} ATTEMPTS WITH LENGTH {other_l} AND ABEGO {other_abego}")
                log.info("ARCHES: COMPLETED")
                arches.append(l)
                arches.append(other_l)
                return pose2, nres_total+l+other_l, arches
                break

            elif i == 9:
                os.remove(bp_backbone)
                os.remove(bp_arch1)
                os.remove(bp_arch2)
                log.info("ARCHES2_3/4_5: FAILED")
