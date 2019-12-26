# This makefile will test the expected against observed output.
# If the make command fails, then there is a problem

all:
	flutile represent --max-day-sep=60 --min-pident-sep=100 1B.2.1.fna.aln > .obs-reps
	diff .obs-reps .exp-reps
