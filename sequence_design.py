import os
from pyrosetta import *

pyrosetta.init()


def prepare_xml(nres1, nres2, hairpins, arches, nres_total, suffix):

	LIBRARY_PATH = os.path.join(os.path.realpath(__file__).replace("sequence_design.py",""), "library")
	output = f"sandwich_{suffix}.xml"

	beta1 = f"1-{nres2}"
	hairpin1 = f"{nres2+1}-{nres2+hairpins[0]}"
	beta2 = f"{nres2+hairpins[0]+1}-{2*nres2+hairpins[0]}"
	arch1 = f"{2*nres2+hairpins[0]+1}-{2*nres2+hairpins[0]+arches[1]}"
	beta3 = f"{2*nres2+hairpins[0]+arches[1]+1}-{2*nres2+hairpins[0]+arches[1]+nres1}"
	hairpin2 = f"{2*nres2+hairpins[0]+arches[1]+nres1+1}-{2*nres2+hairpins[0]+arches[1]+nres1+hairpins[1]}"
	beta4 = f"{2*nres2+hairpins[0]+arches[1]+nres1+hairpins[1]+1}-{2*nres2+hairpins[0]+arches[1]+2*nres1+hairpins[1]}"
	arch2 = f"{2*nres2+hairpins[0]+arches[1]+2*nres1+hairpins[1]+arches[2]}-{2*nres2+hairpins[0]+arches[1]+2*nres1+hairpins[1]+arches[2]}"
	beta5 = f"{2*nres2+hairpins[0]+arches[1]+2*nres1+hairpins[1]+arches[2]+1}-{3*nres2+hairpins[0]+arches[1]+2*nres1+hairpins[1]+arches[2]}"
	arch3 = f"{3*nres2+hairpins[0]+arches[1]+2*nres1+hairpins[1]+arches[2]+1}-{3*nres2+hairpins[0]+arches[1]+2*nres1+hairpins[1]+arches[2]+arches[0]}"
	beta6 = f"{3*nres2+hairpins[0]+arches[1]+2*nres1+hairpins[1]+arches[2]+arches[0]+1}-{3*nres2+hairpins[0]+arches[1]+3*nres1+hairpins[1]+arches[2]+arches[0]}"
	hairpin3 = f"{3*nres2+hairpins[0]+arches[1]+3*nres1+hairpins[1]+arches[2]+arches[0]+1}-{3*nres2+hairpins[0]+arches[1]+3*nres1+hairpins[1]+arches[2]+arches[0]+hairpins[2]}"
	beta7 = f"{3*nres2+hairpins[0]+arches[1]+3*nres1+hairpins[1]+arches[2]+arches[0]+hairpins[2]+1}-{nres_total}"

	strands = ",".join([beta1,beta2,beta3,beta4,beta5,beta6,beta7])
	loops = ",".join([hairpin1,hairpin2,hairpin3,arch1,arch2,arch3])

	with open(os.path.join(LIBRARY_PATH, "xml.xml"), "r") as template:
		lines = [line for line in template]
		for i,line in enumerate(lines):
			if "STRANDS" in line:
				s = line.replace("STRANDS", strands)
				lines[i] = s
			if "LOOPS" in line:
				l = line.replace("LOOPS", loops)
				lines[i] = l

	with open(output, "w") as outxml:
		for line in lines:
			outxml.write(line)


def design_sequence(pose, report_file, xml_file, nres_total):

	xml_object = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects()
	protocol = xml_object.create_from_file(xml_file)
	fastdesign = protocol.get_mover("fast_design")
	score = protocol.get_filter("total_score")
	hpatch = protocol.get_filter("hpatch_score")
	hydrophobics = protocol.get_filter("exposed_hydrophobics")
	total_bsa = protocol.get_filter("total_bsa")

	fastdesign.apply(pose)

	sc = score.score(pose)
	sc_norm = sc / nres_total
	hp = hpatch.score(pose)
	hy = hydrophobics.score(pose)
	tb = total_bsa.score(pose)

	metrics = [report_file.replace("out","pdb"), sc, sc_norm, hp, hy, tb]

	with open(report_file, "w") as report:
		report.write("    ".join(str(m) for m in metrics))

	return pose
