#!/usr/bin/env python

import re
import sys
import csv
from chembl_webresource_client.new_client import new_client
molecule = new_client.molecule

filename = str(sys.argv[1])

with open(filename, 'r') as f:
  reader = csv.DictReader(f, delimiter="\t")
  for row in reader:
    print(row['nci_concept_display_name'])
    res = molecule.search(row['nci_concept_display_name'])
    print(res)
