import numpy as np
from constants import DISTANCE_STRANDS, RESIDUE_LENGTH, CO_BOND, CC_BOND
from functions import caca_vector, co_vector, orientation_vector, rotation_matrix
from pyrosetta import *


pyrosetta.init(silent=True)


def sheet_numbering(parameters):

	res_beta = list(list(zip(*parameters))[2])
	count = 0
	strands = []
    
	for p in parameters:
		count += int(p[2])
	residues = list(range(1, count+1))

	for i in res_beta:
		beta = []
		for j in range(int(i)):
			beta.append(residues[0])
			residues.pop(0)
		strands.append(beta)

	return strands


def beta_strand(strand, phi=-140.0, psi=130.0):

	pose = pose_from_sequence(strand)

	for uid,i in enumerate(strand):
		res_id = uid + 1
        # phi = random.uniform(-150.0, -130.0)
        # psi = random.uniform(120.0, 140.0)
		pose.set_phi(res_id, phi)
		pose.set_psi(res_id, psi)

	return pose


def parallel_strand(pose, strand, nbeta, nres):

	ref_pose = pose.clone()
	add_pose = beta_strand(strand)
	trans_vector = co_vector(ref_pose)

	pyrosetta.bindings.pose.translate(add_pose, trans_vector * (DISTANCE_STRANDS*(nbeta-1)))

	pyrosetta.rosetta.core.pose.append_pose_to_pose(pose, add_pose)

	return pose


def antiparallel_strand(pose, strand, nbeta, nres, quadruple = True):

	length = RESIDUE_LENGTH * nres

	ref_pose = beta_strand(strand)
	trans_vector = co_vector(ref_pose)
	back_vector = caca_vector(ref_pose)
	orien_vector = orientation_vector(ref_pose)
	rotation = rotation_matrix(trans_vector)
	add_pose = beta_strand(strand)

	if (nres % 2) == 0:
		orien_rotation = rotation_matrix(orien_vector, theta=90)
		trans_vector = co_vector(ref_pose)
	elif (nres % 2) != 0:
		orien_rotation = rotation_matrix(orien_vector, theta=-90)
	
	pyrosetta.bindings.pose.rotate(add_pose, orien_rotation)
	pyrosetta.bindings.pose.rotate(add_pose, rotation)
	pyrosetta.bindings.pose.translate(add_pose, np.array(back_vector) * length)
	pyrosetta.bindings.pose.translate(add_pose, np.array(trans_vector) * (DISTANCE_STRANDS*(nbeta-1)))

	pyrosetta.rosetta.core.pose.append_pose_to_pose(pose, add_pose)

	return pose


def tidy_pdb(in_pdb, nres1, nres2, ori1, ori2):

    chid = "A"

    strand1 = list(range(nres2+1))[1:]
    strand2 = list(range(strand1[-1]+1, strand1[-1]+1+nres2))
    strand3 = list(range(strand2[-1]+1, strand2[-1]+1+nres1))
    strand4 = list(range(strand3[-1]+1, strand3[-1]+1+nres1))
    strand5 = list(range(strand4[-1]+1, strand4[-1]+1+nres2))
    strand6 = list(range(strand5[-1]+1, strand5[-1]+1+nres1))
    strand7 = list(range(strand6[-1]+1, strand6[-1]+1+nres1))

    topology = [strand7, strand6, strand3, strand4, strand1, strand2, strand5]
    topology = [item for sublist in topology for item in sublist]

    ini_res = 1
    end_res = strand7[-1]
    res_lis = list(range(ini_res, end_res+1))
    renumbered_lis = []
    sorted_lis = []
    
    with open(in_pdb, "r") as inpdb:
    	for line in inpdb:
    		for i, resid in enumerate(res_lis):
    			if line.startswith("ATOM"):
    				if line[24:26].strip() == str(resid):
    					renumbered_lis.append(line[:21] + chid + str(topology[i]).rjust(4) + line[26:])

    for line in sorted(renumbered_lis, key=lambda line: line[22:26]):
    	sorted_lis.append(line)

    return sorted_lis


def sort_pdb(in_pdb, out_pdb):

	with open(in_pdb, "r") as lines:
		with open(out_pdb, "w") as outfile:
			for line in sorted(lines, key=lambda line: line[22:26]):
				outfile.write(line)


def arch_distance(sandwich_pdb):

    pose = pose_from_pdb(sandwich_pdb)

    xyz_a = pose.residue(18).xyz("C")
    xyz_b = pose.residue(19).xyz("N")
    d2_3 = math.sqrt(((xyz_b[0]-xyz_a[0])**2 + (xyz_b[1]-xyz_a[1])**2 + (xyz_b[2]-xyz_a[2])**2))

    xyz_a = pose.residue(36).xyz("C")
    xyz_b = pose.residue(37).xyz("N")
    d4_5 = math.sqrt(((xyz_b[0]-xyz_a[0])**2 + (xyz_b[1]-xyz_a[1])**2 + (xyz_b[2]-xyz_a[2])**2))

    xyz_a = pose.residue(45).xyz("C")
    xyz_b = pose.residue(46).xyz("N")
    d5_6 = math.sqrt(((xyz_b[0]-xyz_a[0])**2 + (xyz_b[1]-xyz_a[1])**2 + (xyz_b[2]-xyz_a[2])**2))

    return str(round(d2_3,1)), str(round(d4_5,1)), str(round(d5_6,1))


def split(word):
    return [char for char in word]


















	
