#!/usr/bin/env python

import re
import sys

filename = str(sys.argv[1])

f = open(filename,"r")
nci_t = "NA"
umls_cui = "NA"
label = "NA"
for line in f:
   if re.search(r'<owl:Class rdf',line):
      nci_t = re.sub('^\s+|<owl:Class rdf:about="http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#|">','',line.rstrip())
   if re.search(r'<rdfs:label>',line):
      label = re.sub(r'^\s+|<rdfs:label>|</rdfs:label>','',line.rstrip())
      if nci_t != "NA":
         print(str(nci_t) + '\t' + str(label) + '\t' + str(umls_cui))
      umls_cui = "NA"
      label = "NA"
      nci_t = "NA"
   if re.search(r'<P207',line):
      umls_cui = re.sub(r'<P207>|</P207>','',line.rstrip())
   
