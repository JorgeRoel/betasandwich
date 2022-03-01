import math
import numpy as np
from pyrosetta import *

pyrosetta.init(silent=True)


def caca_vector(pose):

	caca = np.array(pose.residue(3).xyz("CA") - pose.residue(1).xyz("CA"))

	return caca / np.linalg.norm(caca)


def co_vector(pose, reverse=False):

	return np.array(pose.residue(1).xyz("O") - pose.residue(1).xyz("C"))


def orientation_vector(pose):

	com1 = pyrosetta.rosetta.protocols.geometry.center_of_mass(pose, 1, 2)
	com2 = pyrosetta.rosetta.protocols.geometry.center_of_mass(pose, 3, 4)

	return np.array(com2 - com1)


def rotation_matrix(axis, theta=180):

    axis = np.asarray(axis)
    theta = np.asarray(theta * math.pi / 180.0 )
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])




