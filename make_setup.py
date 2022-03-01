import sys
import os


number_models = int(sys.argv[1])
code = sys.argv[2]
rotation = sys.argv[3]

with open("job.txt", "w") as outfile:
	for i in range(number_models):
		suffix = '{0:04}'.format(i+1)
		cmd = "python3 /data/jrtcri/PROJECTS/design/betasandwich/beta_sandwich.py -code {} -rotation {} -suffix {} > /tmp/trash 2> /tmp/trash;".format(code, rotation, suffix)
		outfile.write(cmd)
		outfile.write(os.linesep)
