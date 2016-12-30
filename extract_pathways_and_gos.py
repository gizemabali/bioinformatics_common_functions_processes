import urllib
import urllib.request
import re


# page = urllib.request.urlopen('http://www.uniprot.org/uniprot/P10316.txt/').read()
# page = page.split(b'\n')
# for i in x:
#   print(i.decode('utf-8'))
protein_cellular = open('Output/common_cellular_gos.txt', 'w')
protein_biological = open('Output/common_biological_gos.txt', 'w')
protein_molecular = open('Output/common_molecular_gos.txt', 'w')
protein_reactome = open('Output/common_reactome.txt', 'w')
protein_biocyc = open('Output/common_biocyc.txt', 'w')
protein_signor = open('Output/common_signor.txt', 'w')
protein_signalink = open('Output/common_signalink.txt', 'w')
protein_unipathway = open('Output/common_unipathway.txt', 'w')
protein_enzyme = open('Output/common_enzyme.txt', 'w')
protein_brenda = open('Output/common_brenda.txt', 'w')
protein_sabio = open('Output/common_sabio.txt', 'w')
class Get_Pages:
        
    def __init__(self, urls): 
       self.urls = urls
       #self.protein_information = {}
       self.protein_features = {}
       self.sets = []
       self.protein = ""
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
                #print(self.protein)
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
                      
    
       self.common_functions()	
       print(self.sets)
            

    def get_commons(self, index1, index2, keyword, file_name):
            p1 = list(self.protein_features.items())[index1][0]
            p2 = list(self.protein_features.items())[index2][0]
            p1_go = list(self.protein_features.items())[index1][1][keyword]
            p2_go = list(self.protein_features.items())[index2][1][keyword]
            p1_len = len(p1_go)
            p2_len = len(p2_go)
            commons = set(p1_go) & set(p2_go)
            common_count = len(commons) 
            el = []
            el.append(str(p1))
            el.append(str(p2))
            
            if keyword=='go_cellular':
                self.sets.append({'sets':el, 'size':common_count})
            if len(p1_go) == 0:
                p1_go = "There is no " + keyword + " for " + str(p1)
                p1_len = 0
            if len(p2_go) == 0:
                p2_go = "There is no " + keyword + " for " + str(p2)
                p2_len = 0
            if len(commons) == 0:
                commons = "There is no common " + keyword + " for " + str(p1) + " between " + str(p2)
                common_count = 0
            
            file_name.write('-----------------' + p1 + ' & ' + p2 + '------------------------'+ '\n\n')
            file_name.write(p1 +'\t' + 'cellular go count:' + '\t' + str(p1_len) +'\t' + str(p1_go) + '\n')
            file_name.write(p2 +'\t' + 'cellular go count:' + '\t' + str(p2_len) + '\t' + str(p2_go) + '\n')
            file_name.write(p1 + ' & ' + p2 +'\t' + 'common cellular go count:' + '\t' + str(common_count) + '\n')
            file_name.write(p1 + ' & ' + p2 +'\t' + 'common cellular gos:' + '\t' + str(commons) + '\n')
            file_name.write('\n\n')

    def common_functions(self):
       commons = []
       i = 0
       arr = []
  
       for i in range(len(self.protein_features)-1):
         #print(key, self.protein_features[key]['go_cellular'])
         self.sets.append({'sets':[list(self.protein_features.items())[i][0]], 
                    'size':len(list(self.protein_features.items())[i][1]['go_cellular'])})
         for k in range(i+1, len(self.protein_features)):
            self.sets.append({'sets':[list(self.protein_features.items())[k][0]], 
                    'size':len(list(self.protein_features.items())[k][1]['go_cellular'])})
            self.get_commons(i, k, 'go_biological', protein_biological)
            self.get_commons(i, k, 'go_cellular', protein_cellular)
            self.get_commons(i, k, 'go_molecular', protein_molecular)
            self.get_commons(i, k, 'reactome', protein_reactome)
            self.get_commons(i, k, 'biocyc', protein_biocyc)
            self.get_commons(i, k, 'signor', protein_signor)
            self.get_commons(i, k, 'signalink', protein_signalink)
            self.get_commons(i, k, 'unipathway', protein_unipathway)
            self.get_commons(i, k, 'enzyme', protein_enzyme)
            self.get_commons(i, k, 'brenda', protein_brenda)
            self.get_commons(i, k, 'sabio', protein_sabio)


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
