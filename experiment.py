from ahocorasick import AhoKeywordTree
import time
start_time_1 = time.time()
from dataset import load_data
from skt import SuffixKeywordTree


patterns_1,text_1 = load_data("2_P50L19M30_T1L31304M31304_Q70.fasta")
# pattern_1 = ['he', 'hers', 'his', 'she']




skt = SuffixKeywordTree(text_1)
for kwd in patterns_1:
        res = skt.search_kwd(kwd)
def Merge(dict1, dict2): 
    return(dict2.update(dict1))
# print(text_1)

# kwt_indices_1 = {"a":0}
# for pattern in patterns_1:
kwt = AhoKeywordTree(patterns_1)
curr_kwt_indices_1 = kwt.find_keywords(text_1[0])
print(curr_kwt_indices_1)
	# kwt_indices_1 =  Merge(kwt_indices_1, curr_kwt_indices_1)



# add for loops for either text or pattern

print("first test")
# print(kwt_indices_1)
print("--- %s seconds ---" % (time.time() - start_time_1))

#### part 2 



pattern_2,texts_2 = load_data("6_P60L58M90_T20L29345M45582_Q339.fasta")
# print(texts_2)

###
start_time_2 = time.time()
pattern_2 = ['SMEAR', 'EARN', 'ARNOLD', 'OLDER', 'DERIVE', 'VENO']
text_2 = 'BSMEARNOLDERIVENO'
kwt = AhoKeywordTree(pattern_2)
kw_indices_2 = kwt.find_keywords(text_2)
kwt = AhoKeywordTree(pattern_2)
kwt_indices_1 = {'a':0}
print(kwt)
for text in texts_2:
	curr_kwt_indices_1 = kwt.find_keywords(text)
	# print(curr_kwt_indices_1)
	# kwt_indices_1 =  Merge(kwt_indices_1, curr_kwt_indices_1)


	##part 4 
# for text in texts_2:
# 	skt = SuffixKeywordTree(text)
# 	print(skt)

# 	for kwd in pattern_2:
#         	res = skt.search_kwd(kwd)
#         	# print(kwd)
        	# print(res)
# print(text_1)	
print("second test")
# print(kw_indices_2)
print("--- %s seconds ---" % (time.time() - start_time_2))
