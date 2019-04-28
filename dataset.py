from Bio import SeqIO

def load_data(filename):
	"""
	Load patterns and texts from a FASTA dataset file
	:params: filename, a string
	:return: patterns and texts in the file, both list of strings
	"""
	patterns = []
	texts = []
	fasta_sequences = SeqIO.parse(open(filename),'fasta')
	for record in fasta_sequences:
		record_id, sequence = record.id, str(record.seq)
		if (record_id[0] == 'P'):
			patterns.append(sequence)

		if (record_id[0] == 'T'):
			texts.append(sequence)

	return patterns, texts

if __name__ == '__main__':
	patterns, texts = load_data('example.fasta')
	print(patterns)
	print(texts)