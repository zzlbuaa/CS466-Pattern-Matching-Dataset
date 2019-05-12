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
	pattern_len_sum = 0
	max_pattern_len = 0
	for _ in range(num_pattern):
		cur_pattern_len = len_pattern + random.randint(-len_pattern // 2, len_pattern // 2)
		pattern_len_sum += cur_pattern_len
		patterns.append(get_seq(cur_pattern_len))
		max_pattern_len = max(max_pattern_len, cur_pattern_len)
	avg_pattern_len = pattern_len_sum // len(patterns)

	texts = []
	text_len_sum = 0
	num_total_queries = 0
	max_text_len = 0
	for _ in range(num_text):
		cur_base_text_len = len_base_text + random.randint(-len_base_text // 2, len_base_text // 2)
		base_text = get_seq(cur_base_text_len)
		cur_num_insert = num_insert + random.randint(-num_insert // 2, num_insert // 2)
		cur_text = get_text(base_text, patterns, cur_num_insert)
		num_total_queries += cur_num_insert
		texts.append(cur_text)
		text_len_sum += len(cur_text)
		max_text_len = max(max_text_len, len(cur_text))
	avg_text_len = text_len_sum // len(texts)

	filename = "P{}L{}M{}_T{}L{}M{}_Q{}.fasta"
	filename = filename.format(num_pattern, avg_pattern_len, max_pattern_len, num_text, avg_text_len, max_text_len, num_total_queries)
	save_dataset(patterns, texts, filename)

# get_dataset(100, 10, 1, 100000, 200)
get_dataset(30, 30, 50, 10000, 12)
get_dataset(10, 50, 60, 8000, 5)
get_dataset(60, 60, 20, 30000, 18)
# get_dataset(100, 50, 20, 10000, 10)
# get_dataset(50, 20, 1, 20000, 50)


