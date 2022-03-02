import numpy as np
from pyrosetta import *

pyrosetta.init()


def check_chirality(pose, nres1, nres2, reverse=False):

	"Check chiralities of hairpins as defined in Lin. PNAS. 2015"

	chiralities = {

					"h1": "",
					"h2": "",
					"h3": ""
	}

	# Check h1
	com1 = pyrosetta.rosetta.protocols.geometry.center_of_mass(pose, 1, 2)
	com2 = pyrosetta.rosetta.protocols.geometry.center_of_mass(pose, nres2-1, nres2)
	u = np.array(com2 - com1)
	u = u / np.linalg.norm(u)

	com1 = pyrosetta.rosetta.protocols.geometry.center_of_mass(pose, 1, nres2)
	com2 = pyrosetta.rosetta.protocols.geometry.center_of_mass(pose, nres2+1, nres2*2)
	v = np.array(com2 - com1)
	v = v / np.linalg.norm(v)

	uv = np.cross(u, v)
	uv = uv / np.linalg.norm(uv)

	ca = pose.residue(1).xyz("CA")
	cb = pose.residue(1).xyz("CB")
	w = np.array(cb - ca)
	w = w / np.linalg.norm(w)

	if np.dot(uv,w) > 0:
		chiralities["h1"] = "R"
	else:
		chiralities["h1"] = "L"

	# Check h2
	com1 = pyrosetta.rosetta.protocols.geometry.center_of_mass(pose, nres2*2+1, nres2*2+2)
	com2 = pyrosetta.rosetta.protocols.geometry.center_of_mass(pose, nres2*2+nres1-1, nres2*2+nres1)
	u = np.array(com2 - com1)
	u = u / np.linalg.norm(u)

	com1 = pyrosetta.rosetta.protocols.geometry.center_of_mass(pose, nres2*2+1, nres2*2+nres1)
	com2 = pyrosetta.rosetta.protocols.geometry.center_of_mass(pose, nres2*2+nres1+1, nres2*2+nres1*2)
	v = np.array(com2 - com1)
	v = v / np.linalg.norm(v)

	uv = np.cross(u, v)
	uv = uv / np.linalg.norm(uv)

	ca = pose.residue(nres2*2+nres1).xyz("CA")
	cb = pose.residue(nres2*2+nres1).xyz("CB")
	w = np.array(cb - ca)
	w = w / np.linalg.norm(w)

	if np.dot(uv,w) > 0:
		chiralities["h2"] = "R"
	else:
		chiralities["h2"] = "L"

	# Check h3
	com1 = pyrosetta.rosetta.protocols.geometry.center_of_mass(pose, nres2*3+nres1*2+1, nres2*3+nres1*2+2)
	com2 = pyrosetta.rosetta.protocols.geometry.center_of_mass(pose, nres2*3+nres1*3-1, nres2*3+nres1*3)
	u = np.array(com2 - com1)
	u = u / np.linalg.norm(u)

	com1 = pyrosetta.rosetta.protocols.geometry.center_of_mass(pose, nres2*3+nres1*2+1, nres2*3+nres1*3)
	com2 = pyrosetta.rosetta.protocols.geometry.center_of_mass(pose, nres2*3+nres1*3+1, nres2*3+nres1*4)
	v = np.array(com2 - com1)
	v = v / np.linalg.norm(v)

	uv = np.cross(u, v)
	uv = uv / np.linalg.norm(uv)

	ca = pose.residue(nres2*3+nres1*3).xyz("CA")
	cb = pose.residue(nres2*3+nres1*3).xyz("CB")
	w = np.array(cb - ca)
	w = w / np.linalg.norm(w)

	if np.dot(uv,w) > 0:
		chiralities["h3"] = "R"
	else:
		chiralities["h3"] = "L"

	return chiralities


def check_short_hairpins(pose, idx):

	lh = 2
	abego = list(rosetta.core.sequence.get_abego(pose))

	aa = ["A", "A"]
	bg = ["B", "G"]
	ea = ["E", "A"]
	gg = ["G", "G"]

	h = abego[idx:idx+lh]

	if h != aa and h != bg and h != ea and h != gg:
		return False
	else:
		return True, h


def check_long_hairpins(pose, idx):

	lh = 5
	abego = list(rosetta.core.sequence.get_abego(pose))

	aa = ["B", "A", "A", "G", "B"]

	h = abego[idx:idx+lh]

	if h != aa:
		return False
	else:
		return True, h
