import urllib
import urllib.parse
import urllib.request
from urllib.request import urlopen
import sys, string, os
import re
from extract_pathways_and_gos import Get_Pages

protein_list = []
search_key = ''
ensembl_query = ''
refseq_query = ''
uniprot_query = ''
ensembl_key =''
refseq_key = ''
uniprot_key = ''
""" change initial parameters with uniprot keywords """
##print(len(sys.argv))
### python uniprot_acc_id_finder.py --ensembl 'ENSG00000141510'
try:
  for i in range(len(sys.argv)):
     if sys.argv[i] == "--ensembl":
       ensembl_query = sys.argv[i+1].split(',')
     if sys.argv[i] == "--refseq":
       refseq_query = sys.argv[i+1].split(',')
     if sys.argv[i] == "--uniprot_id":
       uniprot_query = sys.argv[i+1].split(',')
except:
    print("Please write the command in a structure like\npython uniprot_acc_id_finder.py --ensembl 'ENSG00000167110'")
dic = {}
dic["ENSEMBL_ID"] = ensembl_query
dic["P_REFSEQ_AC"] = refseq_query
dic["ID"] = uniprot_query
class Find_Protein_ACC():

     def __init__(self, protein_id, protein_list):

         self.url = 'http://www.uniprot.org/uploadlists/'
         self.uniprot_acc = []
         for key,value in dic.items():
             self.connect_to_uniprot(key,value)
             if value != '':
                print(key + " " + "=======>" + " " + str(value))
             elif value == '':
                print(key+ " " + "=======>" + " " + "No" + " " + key)
         self.get_page_urls = [] 
         self.uniprot_acc = list(set(self.uniprot_acc))
         
         for url in self.uniprot_acc:
             self.get_page_urls.append('http://www.uniprot.org/uniprot/' + str(url) + '.txt/')
             #print('http://www.uniprot.org/uniprot/' + str(url) + '.txt/')
         get_pages = Get_Pages(self.get_page_urls)

     def connect_to_uniprot(self,search_key, arr):
         query = " ".join(arr)
         self.params = {
               'from':search_key,
               'to':'ACC',
               'format':'fasta',
               'query':query
         }
         self.page = ''
         data = urllib.parse.urlencode(self.params)
         request = urllib.request.Request(self.url, data.encode())
         request.add_header('User-Agent','Python')
         response = urlopen(request)
         self.page = response.read(200000)
         self.page = self.page.decode().split('\t')
         self.page = ' '.join(self.page)
         self.get_acc_information()

     def get_acc_information(self):
         pattern = '>(sp|tr)' + '(.*)' + '\s.*'       
         for every_line in self.page.split('\n'):    
             try: 
                info = re.search(pattern, every_line)
                info1 = info.group(2)
                acc_information = info1.split("|")[1]
                self.uniprot_acc.append(acc_information)
                k += 1
             except:
                pass

Find_Protein_ACC(search_key, protein_list)   
