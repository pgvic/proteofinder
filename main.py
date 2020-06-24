#!/usr/bin/env python
#Victor Paton Gonzalez Programación para Bioinformática

########### MÓDULOS ############
import os
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
import subprocess as sp

import blast as bls
import align as al
import prosite as pst
import folder as fl

#GUI
import tkinter.filedialog
import tkinter as tk
from tkinter import ttk
from PIL import ImageTk, Image



#FUNCIONES

def default():
	'''
	Establece los parámetros de pident y qcovs por defecto
	'''
	pident = 25
	qcovs = 50
	return(pident, qcovs)

def genbank_to_multifasta (genbank_files):
	'''
	Convierte un arcivo Genbank en multifasta. Recoge el locus tag, el nombre y
	la secuencia proteica (CDS).
	'''
	output_handle = open("multifasta.txt", "w")
	for i in genbank_files:
		input_handle  = open(i, "r")
		for seq_record in SeqIO.parse(input_handle, "genbank"):
		    for seq_feature in seq_record.features :
		        if seq_feature.type=="CDS" :
		            try:
			            output_handle.write(">%s@%s\n%s\n" % (
			                   seq_feature.qualifiers['locus_tag'][0],
			                   seq_record.annotations['organism'].replace(" ","_"),
			                   seq_feature.qualifiers['translation'][0]))
		            except:
		            	pass
		input_handle.close()
	output_handle.close()
	return

def blaster(input, pident, qcovs):
	'''
	Iterando sobre un archivo fasta, realiza un blastp de cada query contra un
	archivo multifasta que actúa como subject. Además, realiza un árbol filogenético
	(formatos .nw y .png) y busca dominios proteicos en cada una de las secuencias
	resultantes del blast (tanto query como subject)
	'''
	dir_results = fl.main_folder()
	for record in SeqIO.parse(input, "fasta"):
		dir_query=fl.query_folder(dir_results, record.id)

		#Blast y filtrado
		bls.blast(record, qcovs, pident, dir_query, "multifasta.txt")
		blast_aligned = record.id + "aligned.fasta"
		blast_fasta = record.id + "blast_fasta.fa"

		#Alineamiento y creación de árboles, si estos contienen más de una secuencia
		try:
			al.alignment(blast_fasta, blast_aligned)
			nw_tree = dir_query + record.id + ".nw"
			al.tree(blast_aligned, nw_tree)
			tree_img = dir_query + record.id + ".png"
			if save_pngtree.get()==1:
				al.tree_drawer(nw_tree, tree_img)
		except:
			print("Número de secuencias insuficiente para crear un árbol")

		#Búsqueda de dominios proteicos si hay una db válida
		try:
			domain_file = dir_query + record.id + "_domains.txt"
			pst.domain_scanner(blast_fasta, domain_file, prosite_db)
			print("Base de datos de PROSITE encontrada")
		except:
			print('No se ha introducido ninguna base de datos Prosite, no se buscarán dominios.\n')

		#Eliminar archivos temporales
		os.remove(blast_aligned)
		os.remove(blast_fasta)
	return

def blast_parser_to_fasta (ids, output, query_fasta, multifasta):
	'''
	Crea un archivo con las secuencias completas de una lista de IDs procedentes
	de un Blast. Además, añade la secuencia query que se ha usado para hacer Blast.
	'''
	with open(query_fasta, "r") as query, open(multifasta, "r") as multi, open(output, "w") as output:
		for record in SeqIO.parse(multi, "fasta"):
			for id in ids:
				if id == record.id:
					output.write(">{}\n{}\n\n".format(record.id, record.seq))
				else:
					pass
		for line in query:
			output.write(line)
	return


#Selección de archivos
def genbank_window():
	'''
	Abre una ventana de selección de archivos para seleccionar los archivos subject.
	Filtra por extensiones de formato GenBank y .txt para evitar formatos incompatibles
	'''
	directorio = os.getcwd()
	global genbank_files
	genbank_files = tk.filedialog.askopenfilenames(
												parent=root,
												initialdir = directorio,
												title = "Seleccione archivos Genbank",
												filetypes = (("Genbank files","*.gbff *.txt"),("Todos los archivos","*.*"))
												)
	g_files = str(len(genbank_files))
	g_added = tk.Label(root, text=g_files + " files added", font =("Times New Roman", 10))
	g_added.place(x=150, y=453)
	return

def query_window():
	'''
	Abre una ventana de selección de archivos para seleccionar los archivos query.
	Filtra por extensiones de formato fasta y .txt para evitar formatos incompatibles
	'''
	global query_files
	query_files = tk.filedialog.askopenfilenames(
												parent=root,
												initialdir = directorio,
												title = "Seleccione archivos fasta",
												filetypes = (("Fasta files","*.fasta *.fa *.mpfa *.fna *.fsa *.fas *.faa *.txt"),("Todos los archivos","*.*"))
												)
	q_files = str(len(query_files))
	q_added = tk.Label(root, text=q_files + " files added", font =("Times New Roman", 10))
	q_added.place(x=150, y=533)
	return

def prosite_window():
	'''
	Abre una ventana de selección de archivos para seleccionar los archivos
	que se van a usar como base de datos de dominios.
	Filtra por extensiones de formato .dat y .txt para evitar formatos incompatibles
	'''
	global prosite_db
	prosite_db = tk.filedialog.askopenfilename(
												parent=root,
												initialdir = directorio,
												title = "Seleccione la base de datos de Prosite",
												filetypes = (("Prosite database","*.dat .txt"),("Todos los archivos","*.*"))
												)
	p_files = str(len(prosite_db))
	p_added = tk.Label(root, text="Database added", font =("Times New Roman", 10))
	p_added.place(x=150, y=613)
	return


def gui_msg(texto, col="black"):
	'''
	Mensaje de texto que aparece en la interfaz gráfica
	'''
	msg=tk.Label(root, text=texto, font =("Times New Roman", 10), fg = col, bg = 'gold2')
	msg.place(x=20, y=350)
	return

#Función principal
def launch_function():
	'''
	Lanza el programa. Avisa cuando éste ha terminado, y si lo ha hecho de
	forma exitosa o debido a un error.
	'''
	try:
		pident = float(pident_txtbox.get())
		qcovs = float(qcovs_txtbox.get())
	except:
		pident, qcovs= default()
		print("Usando valores por defecto, no se han introducido valores válidos")
	print("Usando valores de pident "+str(pident)+" y qcovs "+str(qcovs))
	try:
		genbank_to_multifasta(genbank_files)
		for i in query_files:
			blaster(i, pident, qcovs)
		os.remove("multifasta.txt")
		texto="¡Análisis completado!"
		gui_msg(texto)
	except:
		texto="Ha ocurrido un error"
		gui_msg(texto, "red")
	return

########GUI#########
#Ventana principal
directorio = os.getcwd()
root = tk.Tk()
root.geometry('1600x900')
root.resizable(height = False, width = False)
root.title("Proteofinder")

#Icono de ventana
icon = ImageTk.PhotoImage(file=directorio+"/GUI/dna.ico")
root.iconphoto(True, icon)
#Iconos diseñados por <a href="https://www.flaticon.es/autores/freepik"
#title="Freepik">Freepik</a> from <a href="https://www.flaticon.es/"
#title="Flaticon"> www.flaticon.es</a>
root.resizable(False, False)

logo = ImageTk.PhotoImage(Image.open(directorio+"/GUI/dna.png"))
logo_panel=tk.Label(root, image=logo)
logo_panel.place(x=0,y=0)

help = ImageTk.PhotoImage(Image.open(directorio+"/GUI/ayuda.png"))
help_panel=tk.Label(root, image=help)
help_panel.place(x=300,y=0)

#genbank
genbank = tk.Button(root, text="GenBank", command=genbank_window)
genbank.place(x=30, y=450)

#query
query = tk.Button(root, text="Query",command=query_window)
query.place(x=30, y=530)

#pident
pident_txtbox=tk.Entry(root, width = 5)
pident_txtbox.place(x=50, y=703)

#qcovs
qcovs_txtbox=tk.Entry(root, width = 5)
qcovs_txtbox.place(x=50, y=753)

#Guardar arboles png
save_pngtree= tk.IntVar()
check_save = tk.Checkbutton(root, variable = save_pngtree, onvalue = 1, offvalue = 0)
check_save.place(x=40, y= 800)

#Base de datos de prosite
prosite = tk.Button(root, text="Prosite DB",command=prosite_window)
prosite.place(x=30, y=610)

#Launch
launch=tk.Button(root, text="LAUNCH", command=launch_function)
launch.place(x=100, y=840)


root.mainloop()
