#MODULOS
import os

def main_folder():
	'''
	Genera un directorio de resultados. Si este ya existe, añade un distintivo
	(1), (2)...
	'''
	counter=""
	while True:
		try:
			dir_results= os.getcwd()+"/results"+counter
			os.mkdir(dir_results)
			break
		except:
			if counter:
					counter = '('+str(int(counter[1:-1])+1)+')'
			else:
				counter = '(1)'
				pass

	return(dir_results)

def query_folder(dir_results, id):
	dir_name=dir_results+"/"+id+"_results/"
	counter=""
	while True:
		try:
			dir_name=dir_results+"/"+id+"_results/"
			os.makedirs(dir_name)
			break
		except:
			if counter:
					counter = '('+str(int(counter[1:-1])+1)+')'
			else:
				counter = '(1)'
				pass

	return(dir_name)
