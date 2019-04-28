import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

chars = 'ATCG'

def get_seq(length):
	"""
	Generate a sequence of length
	:return: the generated sequence
	"""
	return ''.join(random.choice(chars) for _ in range(length))

def get_text(text, patterns, num_insert):
	"""
	Randomly insert serveral patterns into the base text,
	to generate the text for the dataset
	:param text: the base text, should be a string
	:param patterns: the patterns to insert into base text, should be list of strings
	:param num_insert: the total times of insertions, integer
	:return: the generated text
	"""
	insert_idxs = random.sample(range(len(text)), num_insert)
	insert_idxs.sort()

	inserted_len = 0
	for idx in insert_idxs:
		insert_pattern = random.choice(patterns)
		insert_idx = idx + inserted_len
		text = text[:insert_idx] + insert_pattern + text[insert_idx:]
		inserted_len += len(insert_pattern)

	return text

def save_dataset(patterns, texts, filename):
	"""
	Save patterns and texts to a FASTA file,
	pattern sequence id is initiated with P,
	text sequence id is initiated with T
	"""
	sequences = []
	for idx, pattern in enumerate(patterns):
		sequences.append(SeqRecord(Seq(pattern), id="P"+str(idx), description="PATTERN"))

	for idx, text in enumerate(texts):
		sequences.append(SeqRecord(Seq(text), id="T"+str(idx), description="TEXT"))

	with open(filename, "w") as output_handle:
		SeqIO.write(sequences, output_handle, "fasta")

def get_dataset(num_pattern, len_pattern, num_text, len_base_text, num_insert):
	"""
	Generate a dataset consisting of sets of patterns and sets of texts,
	and save it to current path as a FASTA format file,
	the texts are created by inserting several patterns into base text
	:param num_pattern: total number of patterns to generate, integer
	:param len_pattern: length of each pattern, integer
	:param num_text: total number of texts to generate, integer
	:param len_base_text: length of each base text, integer
	:param num_insert: number of insertions of patterns to each text, integer
	:return none
	"""
	patterns = []
	for _ in range(num_pattern):
		patterns.append(get_seq(len_pattern))

	texts = []
	for _ in range(num_text):
		base_text = get_seq(len_base_text)
		texts.append(get_text(base_text, patterns, num_insert))
	filename = "dataset_P{}L{}_T{}L{}I{}.fasta"
	filename = filename.format(num_pattern, len_pattern, num_text, len_base_text, num_insert)
	save_dataset(patterns, texts, filename)

get_dataset(100, 10, 1, 100000, 200)
get_dataset(100, 50, 1, 100000, 20)
get_dataset(100, 50, 20, 10000, 10)


