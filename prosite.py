from Bio.ExPASy import Prosite,Prodoc
import re
from Bio import SeqIO

def pattern_replace(input):
	'''
	Convierte las regex de Prosite a regex identificables por el módulo re
	'''
	pattern = input
	pattern = pattern.replace('-', '')
	pattern = pattern.replace('{', '[^')
	pattern = pattern.replace('}', ']')
	pattern = pattern.replace('(', '{')
	pattern = pattern.replace(')', '}')
	pattern = pattern.replace('x', '.')
	pattern = pattern.replace('>', '$')
	pattern = pattern.replace('<', '^')
	return pattern

def prosite_data(prosite_db):
	'''
	Extrae del archivo prosite.dat (base de datos de Prosite) el accession, la
	regex que identifica a ese dominio y una breve descripción del mismo
	'''
	ids = []
	patterns = []
	descriptions = []
	with open(prosite_db, "r") as raw:
		records = Prosite.parse(raw)
		for record in records:
			pattern = pattern_replace(record.pattern)
			if pattern != "":
				ids.append(record.accession)
				patterns.append(pattern)
				descriptions.append(record.description)
	return(ids, patterns, descriptions)

def domain_scanner(input, output, prosite_db):
	'''
	Itera sobre las proteínas y sobre ellas intenta encontrar patrones de
	dominios de Prosite. Devuelve un archivo con los dominios encontrados
	'''
	pattern_id=prosite_data(prosite_db)[0]
	pattern_regex=prosite_data(prosite_db)[1]
	description = prosite_data(prosite_db)[2]
	with open(input, "r") as input_handle, open(output, "w") as output_handle:
		for record in SeqIO.parse(input_handle, "fasta"):
			seq = str(record.seq)
			for i in range(0,len(pattern_id)):
				matches = re.finditer(pattern_regex[i], seq)
				for m in matches:
					start = m.start()
					end = m.end()
					seq = m.group()
					output_handle.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
															record.id,
															pattern_id[i],
															start,
															end,
															seq,
															description[i],
															))
	return
