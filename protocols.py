import numpy as np
import os
import random
from blueprint import *
from functions import rotation_matrix
from logger import LoggingManager
from util import beta_strand, antiparallel_strand

from pyrosetta import *

pyrosetta.init()


def sandwich(pose1, pose2, nres1, nres2, ori1, ori2, angle, log_file):

    log = LoggingManager.get_logger("sandwich", file_name=log_file)

    if angle == "+":
        theta = random.uniform(0.0, 60.0)
    elif angle == "-":
        theta = random.uniform(-60.0, 0.0)

    distance = random.uniform(10.0, 12.0)
    pertu_x = random.uniform(-5.0, 5.0)
    pertu_y = random.uniform(-5.0, 5.0)

    sandwich_pose = pose1.clone()

    log.info(f"ANGLE: {theta}")
    log.info(f"TRANSLATION: {distance}")
    log.info(f"PERTURBATION_X: {pertu_x}")
    log.info(f"PERTURBATION_Y: {pertu_y}")

    # Align centers of mass
    com1 = pyrosetta.rosetta.core.pose.get_center_of_mass(pose1)
    com2 = pyrosetta.rosetta.core.pose.get_center_of_mass(pose2)
    trans = np.array(com1 - com2)
    pyrosetta.bindings.pose.translate(pose2, trans)

    # Translate sheet 3 strands with respect 4
    res1 = nres1 + 1 #1st res 2nd strand
    res2 = nres1 + 2 #2nd res 2nd strand
    res3 = nres1 * 2 - 1 #N-1 res 2nd strand
    res4 = nres1 * 2 #N res 2nd strand
    
    ca_1 = pose1.residue(res1).xyz("CA")
    ca_2 = pose1.residue(res2).xyz("CA")
    ca_3 = pose1.residue(res3).xyz("CA")
    ca_4 = pose1.residue(res4).xyz("CA")

    com_1 = [(ca_2[0] + ca_1[0])/2., (ca_2[1] + ca_1[1])/2., (ca_2[2] + ca_1[2])/2.]
    com_2 = [(ca_4[0] + ca_3[0])/2., (ca_4[1] + ca_3[1])/2., (ca_4[2] + ca_3[2])/2.]
    v = [x1 - x2 for (x1, x2) in zip(com_1, com_2)]

    res2 = res4 #N res 2nd strand
    res3 = nres1 * 2 + 1 #1st res 3rd strand
    res4 = nres1 * 3 #N res 3rd strand

    com_3 = pyrosetta.rosetta.protocols.geometry.center_of_mass(pose1, res1, res2)
    com_4 = pyrosetta.rosetta.protocols.geometry.center_of_mass(pose1, res3, res4)
    u = np.array(com_4 - com_3)

    trans = np.cross(u, v)
    trans = trans / np.linalg.norm(trans)
    pyrosetta.bindings.pose.translate(pose2, trans * distance)

    # Rotate and return to "aligned" position
    com_ori = pyrosetta.rosetta.core.pose.get_center_of_mass(pose2)
    rotation = rotation_matrix(trans, theta=theta)
    pyrosetta.bindings.pose.rotate(pose2, rotation)
    com_rot = pyrosetta.rosetta.core.pose.get_center_of_mass(pose2)
    back = np.array(com_ori - com_rot)
    pyrosetta.bindings.pose.translate(pose2, back)

    # Apply pseudo-x/y perturbations to 3 strand-sheet
    v = v / np.linalg.norm(v)
    u = u / np.linalg.norm(u)

    pyrosetta.bindings.pose.translate(pose2, v * pertu_x)
    pyrosetta.bindings.pose.translate(pose2, u * pertu_y)

    pyrosetta.rosetta.core.pose.append_pose_to_pose(sandwich_pose, pose2)

    log.info("BACKBONE COMPLETED")
    
    return sandwich_pose


def blueprint_hairpins(nres1, nres2, lh, outfile, abego_string=None):

    blueprint = Blueprint()
    blueprint.add_segment_to_segment(anchor="L1", ss="E", length=nres2, append=True)
    blueprint.add_segment_to_segment(anchor="E1", ss="L", length=lh, abego=abego_string, append=True)
    blueprint.add_segment_to_segment(anchor="L2", ss="E", length=nres2, append=True)
    blueprint.add_segment_to_segment(anchor="E2", ss="E", length=nres1, append=True)
    blueprint.add_segment_to_segment(anchor="L2", ss="E", length=nres1, append=True)
    blueprint.add_segment_to_segment(anchor="E2", ss="E", length=nres2, append=True)
    blueprint.add_segment_to_segment(anchor="E2", ss="E", length=nres1, append=True)
    blueprint.add_segment_to_segment(anchor="E2", ss="E", length=nres1, append=True)
    blueprint.add_termini()
    blueprint.freeze_all()

    pos = 1
    for res in blueprint.bp_data:
        if res[2][0] == "E":
            res[0] = pos
            pos += 1

    blueprint.remodel_segment(id="L1")
    blueprint.remodel_segment(id="L2")
    blueprint.remodel_segment(id="L3")
    
    blueprint.dump_blueprint(outfile)

