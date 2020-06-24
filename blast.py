import os
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO

#blastp
def blast(record, qcovs, pident, dir_query, subject):
    '''
	Iterando sobre un archivo fasta, realiza un blastp de cada query contra un
	archivo multifasta que actúa como subject. Filtra por valores de pident,
    qcovs y evalue.
    '''
    query_file = record.id + "query_file.fasta"
    with open(query_file, "w") as query:
        query.write(">{}\n{}\n".format(record.id, record.seq))
    blast_nofiltered=record.id + "_blast_nf.tsv"
    blast_filtered = dir_query + record.id + "_blast.tsv"
    blast_cline = NcbiblastpCommandline(
                    cmd = 'blastp',
                    query=query_file,
                    subject=subject,
                    evalue=1e-5,
                    qcov_hsp_perc = qcovs,
                    outfmt = '6 sseqid pident qcovs evalue sseq',
                    out = blast_nofiltered
                    )
    stdout,stderr = blast_cline()

    subject_ids=pident_filter(blast_nofiltered, pident, blast_filtered)
    os.remove(blast_nofiltered)
    blast_fasta = record.id + "blast_fasta.fa"
    blast_parser_to_fasta(subject_ids, blast_fasta, query_file, "multifasta.txt")
    os.remove(query_file)
    return

def pident_filter (input, pident, output):
	ids=[]
	with open (input, "r") as nofilter, open(output, "w") as filter:
		filter.write("ID\tpident\tqcovs\tevalue\tsequence".format())
		for line in nofilter:
			parts = line.split("\t")
			if float(parts[1]) > pident:
				filter.write("{}\t{}\t{}\t{}\t{}".format(
														parts[0],
														parts[1],
														parts[2],
														parts[3],
														parts[4]))
				ids.append(parts[0])
	return(ids)

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

def folder(id):
	'''
	Genera un directorio de resultados. Si este ya existe, añade un distintivo
	(1), (2)...
	'''
	try:
		dir_results= os.getcwd()+"/results"
	except:
		pass
	counter=""
	while True:
		try:
			dir_name=dir_results+"/"+id+"_results" +counter+"/"
			os.makedirs(dir_name)
			break
		except:
			if counter:
				counter = '('+str(int(counter[1:-1])+1)+')'
			else:
				counter = '(1)'
			pass
	return(dir_name)
