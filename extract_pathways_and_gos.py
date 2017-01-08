import urllib
import urllib.request
import re
## python uniprot_acc_id_finder.py --ensembl 'ENSG00000206394','ENSG00000159197','ENSG00000110881' --refseq 'NP_001029186.1'

## python uniprot_acc_id_finder.py --refseq 'NP_000955.1','NP_001169999.1'
## python uniprot_acc_id_finder.py --refseq 'NP_000955.1','NP_001169999.1','NP_001278850.1'
## python uniprot_acc_id_finder.py --refseq 'NP_037305.2','NP_001139633.1'
## python uniprot_acc_id_finder.py --refseq 'NP_000608.1','NP_001139633.1','NP_037305.2'
## python uniprot_acc_id_finder.py --ensembl 'ENSMUSG00000023030'
## python uniprot_acc_id_finder.py --refseq 'NP_388317.1','WP_005082954.1','NP_416893.1'
## python uniprot_acc_id_finder.py --refseq 'NP_388317.1','WP_005082954.1','NP_416893.1','WP_003898649.1'

commons = open('commons_for_proteins.txt', 'w')

fr = open('frequencies.txt', 'w')
class Get_Pages:
        
    def __init__(self, urls): 
       self.urls = urls
       #self.protein_information = {}
       self.protein_features = {}

       self.protein = ""
       self.all_keywords = ['go_biological', 'go_cellular','go_molecular','reactome','biocyc', 'signor', 'signalink', 'unipathway', 'enzyme', 'brenda', 'sabio']
       # create regular exceptions for go annotations and reactome pathways
       self.pattern_cellular = "DR\s*GO;\s*(.*); C:.*\;"
       self.pattern_biological = "DR\s*GO;\s*(.*); P:.*\;"
       self.pattern_molecular = "DR\s*GO;\s*(.*); F:.*\;"
       self.pattern_reactome = "DR\s*Reactome;\s*(.*);.*"
       self.pattern_biocyc = "DR\s*BioCyc;\s*(.*);.*"  
       self.pattern_signor = "DR\s*SIGNOR;\s*(.*);.*"  
       self.pattern_signalink = "DR\s*SignaLink;\s*(.*);.*"  
       self.pattern_unipathway = "DR\s*UniPathway;\s*(.*);.*" 
       self.pattern_enzyme = "DR\s*ENZYME;\s*(.*);.*"  
       self.pattern_brenda = "DR\s*BRENDA;\s*(.*);.*"  
       self.pattern_sabio = "DR\s*SABIO-RK;\s*(.*);.*"   
       for url in self.urls:
          page = urllib.request.urlopen(url).read()
          # decode page to utf-8 from byte 
          page = (page).decode('utf-8')
          # split page with new line to get it line by line           
          info = page.split('\n') 
          #self.protein_information[(url.split('/')[-2])] = page
          go_on = False
          for line in info: 
             if go_on == False:
                self.protein = line.split(" ")
                self.protein = [protein for protein in self.protein if protein != '']
                self.protein = self.protein[1]
           
                self.create_dictionary()
                go_on = True

             else:
                self.get_reactome(line)
                self.get_biocyc(line)
                self.get_signor(line)
                self.get_signalink(line)
                self.get_unipathway(line)
                self.get_enzyme(line)
                self.get_brenda(line)
                self.get_sabio(line)
                self.get_go_annotations(line)
                      
       
       for i in self.all_keywords:
            self.frequency(i)
            self.common_func(i)
       self.close_files()

    def close_files(self):
          fr.close()
          commons.close()
 
    
    def frequency(self, keyword):
        k = []
        for i in range(len(self.protein_features)):
            
            for k_ in list(self.protein_features.items())[i][1][keyword]:
                 k.append(k_)
       
        fr.write("-------" + keyword + "-----------\n\n")
        for ki in list(set(k)):
            count = 0
            for kin in k:
                if ki== kin:
                     count += 1
            
            fr.write(ki + " occuries " + str(count) + " " + "out of " + str(len(self.protein_features)) + " proteins." + "\n")
            count = "%.4f" % (count/len(self.protein_features))
            fr.write("So frequency of" + " " + ki +  " is  " + str(count) + "\n")
        fr.write("\n\n")
  
    def common_func(self, keyword):
        k = []
       
        protein_names = []
        for i in range(len(self.protein_features)):
            protein_names.append(list(self.protein_features.items())[i][0])
           
            for k_ in list(self.protein_features.items())[i][1][keyword]:
                 k.append(k_)
       
        commons.write("-------" + keyword + "-----------\n\n")
        
        count = 0
        for ki in list(set(k)):
           
           proteins = []
           for i in range(len(self.protein_features)):
             
             for k_ in list(self.protein_features.items())[i][1][keyword]:
                  if k_ == ki:
                       proteins.append(list(self.protein_features.items())[i][0])
             
           if len(proteins) == len(self.protein_features):
                count = 1
                commons.write(ki + " occuries in the all proteins " + ",".join(proteins) + "\n")
           
        commons.write("\n")  

    def create_dictionary(self):
        self.protein_features[str(self.protein)] = {}
        self.protein_features[str(self.protein)]['reactome'] = []
        self.protein_features[str(self.protein)]['biocyc'] = []
        self.protein_features[str(self.protein)]['signor'] = []
        self.protein_features[str(self.protein)]['signalink'] = []
        self.protein_features[str(self.protein)]['unipathway'] = []
        self.protein_features[str(self.protein)]['enzyme'] = []
        self.protein_features[str(self.protein)]['brenda'] = []
        self.protein_features[str(self.protein)]['sabio'] = []
        self.protein_features[str(self.protein)]['go_cellular'] = []
        self.protein_features[str(self.protein)]['go_biological'] = []
        self.protein_features[str(self.protein)]['go_molecular'] = []
            
    def get_reactome(self, text):
         reactome = re.search(self.pattern_reactome, text)
         if reactome != None:
            self.protein_features[str(self.protein)]['reactome'].append(reactome.group(1))

    def get_biocyc(self, text):
         biocyc = re.search(self.pattern_biocyc, text)
         if biocyc != None:
            self.protein_features[str(self.protein)]['biocyc'].append(biocyc.group(1))

    def get_signor(self, text):
        signor = re.search(self.pattern_signor, text)
        if signor != None:
            self.protein_features[str(self.protein)]['signor'].append(signor.group(1))
        
    def get_signalink(self, text):
        signalink = re.search(self.pattern_signalink, text)
        if signalink != None:
            self.protein_features[str(self.protein)]['signalink'].append(signalink.group(1))

    def get_unipathway(self, text):
        unipathway = re.search(self.pattern_unipathway, text)
        if unipathway != None:
            self.protein_features[str(self.protein)]['unipathway'].append(unipathway.group(1))

    def get_enzyme(self, text):
        enzyme = re.search(self.pattern_enzyme, text)
        if enzyme != None:
            self.protein_features[str(self.protein)]['enzyme'].append(enzyme.group(1))

    def get_brenda(self, text):
        brenda = re.search(self.pattern_brenda, text)
        if brenda != None:
            self.protein_features[str(self.protein)]['brenda'].append(brenda.group(1))

    def get_sabio(self, text):
        sabio = re.search(self.pattern_sabio, text)
        if sabio != None:
            self.protein_features[str(self.protein)]['sabio'].append(sabio.group(1))

    def get_go_annotations(self, text):
        goC = re.search(self.pattern_cellular, text)
        goB = re.search(self.pattern_biological, text)
        goM = re.search(self.pattern_molecular, text)
        if goC != None:
            self.protein_features[str(self.protein)]['go_cellular'].append(goC.group(1))
        if goB != None:
            self.protein_features[str(self.protein)]['go_biological'].append(goB.group(1))
        if goM != None:
            self.protein_features[str(self.protein)]['go_molecular'].append(goM.group(1))
