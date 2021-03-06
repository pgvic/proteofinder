import matplotlib
import matplotlib.pyplot as plt
from Bio import Phylo
import io
from Bio.Align.Applications import MuscleCommandline
import subprocess as sp



def alignment(input, output):
	'''
	Realiza un proceso de alineamiento usando la herramienta MuscleCommandline
	'''
	muscle_cline= MuscleCommandline(
									input=input,
									out=output
									)
	stdout, stderr = muscle_cline()
	return

def tree(input, output):
	'''
	Crea un arbol filogenético por método Neigbour-joining usando la herramienta
	Muscle.
	'''
	sp.check_call(['muscle',
				'-maketree',
				'-in', input,
				'-out', output,
				'-cluster',
				'neighborjoining'
				])
	return

def tree_drawer(input, output):
	'''
	Crea una imagen png a partir de un arbol filogenético en formato Newick (.nw)
	'''
	tree = Phylo.read(input, "newick")
	matplotlib.rc('font', size=6)
	fig = plt.figure(figsize=(14, 6), dpi=100)
	axes = fig.add_subplot(1,1,1)
	Phylo.draw(tree, axes=axes, do_show=False)
	plt.savefig(output, dpi=300)
	return
