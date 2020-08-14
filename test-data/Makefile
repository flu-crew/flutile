# This makefile will test the expected against observed output.
# If the make command fails, then there is a problem

all:
	flutile represent --max-day-sep=60 --min-pident-sep=100 1B.2.1.fna.aln > .obs-reps
	diff .obs-reps .exp-reps

	flutile represent --min-pident-sep=100 --print-groups 1B.2.1.fna.aln > .obs-groups
	grep -c 'A/' .obs-groups > .A
	grep -c 'A/' 1B.2.1.fna.aln > .B
	diff .A .B
	rm .A .B .obs-groups
