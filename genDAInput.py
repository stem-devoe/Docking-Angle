# webscraping code adapted from:
# https://www.thepythoncode.com/article/extracting-and-submitting-web-page-forms-in-python
# Abdou Rockikz
# accessed Nov 11, 2020

from pymol import cmd
from bs4 import BeautifulSoup
from requests_html import HTMLSession
from urllib.parse import urljoin
import urllib
import urllib.request
from html.parser import HTMLParser
import requests
import re

cmd.reinitialize()
PDB = str(input("Enter PDB ID: "))
contact_dist = input("Enter desired contact distance in Angstroms: ")
contact_dist = float(contact_dist)
cmd.fetch(PDB)
cmd.split_chains()

##########################################################
#### Obtain chain letter for TRA, TRB, MHCa, MHCb/B2M ####
##########################################################
url = "http://rcsb.org/fasta/entry/" + PDB + "/display"
res = requests.get(url)
seq = res.text
seq = str(seq)

MHCtype = None
if "class II" in seq:
    MHCtype = 2
elif "CLASS II" in seq:
    MHCtype = 2
elif "class I" in seq:
    MHCtype = 1
elif "CLASS I" in seq:
    MHCtype = 1

if MHCtype == None:
    print("\nUnable to class of MHC/HLA from FASTA data")
    print("\n", seq)
    MHCtype = input("Manually input class (1 for class I/CD1/MR1 or 2 for class II): ")
    MHCtype = int(MHCtype)

##if MHCtype == 1:
##    print("Class I")
##elif MHCtype == 2:
##    print("Class II")

if "Chains " in seq:
    print("\nThis pMHC:TCR structure is a multimer.\nPlease enter Chain Identifiers manually.\n")
    print(seq)
    print("\nOnly enter Chain Identifiers for one of the pMHC:TCR structures:\n")
    TRAchain = input("Enter TCR alpha chain identifier: ")
    TRBchain = input("Enter TCR beta chain identifier: ")
    MHCa = input("Enter MHC/HLA alpha chain identifier: ")
    MHCb = input("Enter MHC/HLA beta chain or B2M identifier: ")
    print("Enter peptide chain identifier \nif the peptide sequence is not its own chain,")
    peptide = input("enter \"None\" if part of the MHC or \"Error\" if part of the TCR: ")
else:
    print("This pMHC:TCR structure is NOT a multimer.\nChain Identifiers will be detected automatically.")
    chain = seq.split("Chain ")
    TCRlist = []
    index = 0
    for i in chain:
        #print(i)
        if " TCR " in i or "TCR " in i or " TCR" in i:
            TCRlist.append(index)
            #print(1)
        elif "T-cell" in i or "T-Cell" in i:
            TCRlist.append(index)
            #print(2)
        elif "T-CELL" in i:
            TCRlist.append(index)
            #print(3)
        elif "T CELL" in i:
            TCRlist.append(index)
            #print(4)
        elif "T cell" in i:
            TCRlist.append(index)
            #print(5)
        elif "Tcell" in i:
            TCRlist.append(index)
            #print(6)
        elif "TRBV" in i or "TRAV" in i or "Trav" in i or "Trbv" in i:
            TCRlist.append(index)
            #print(7)
        index = index + 1
    # print(TCRlist)
    TRAchain = None
    TRBchain = None
    for i in TCRlist:
        if "alpha" in chain[i]:
            TRAchain = i
        elif "ALPHA" in chain[i]:
            TRAchain = i
        elif "Alpha" in chain[i]:
            TRAchain = i
        elif "A chain" in chain[i]:
            TRAchain = i
        elif "Beta" in chain[i]:
            TRBchain = i
        elif "beta" in chain[i]:
            TRBchain = i
        elif "BETA" in chain[i]:
            TRBchain = i
        elif "B chain" in chain[i]:
            TRBchain = i
    # print(TRAchain)
    MHClist = []
    index = 0
    for i in chain:
        if "histocomp" in i:
            MHClist.append(index)
        elif "MHC" in i:
            MHClist.append(index)
        elif "HLA" in i:
            MHClist.append(index)
        elif "microglobulin" in i or "MICROGLOBULIN" in i:
            MHClist.append(index)
        elif "HISTOCOMP" in i:
            MHClist.append(index)
        elif "H2" in i:
            MHClist.append(index)
        elif "H-2" in i:
            MHClist.append(index)
        index = index + 1
    # print(MHClist)
    MHCa = None
    MHCb = None
    for i in MHClist:
        if MHCtype == 1:
            if "alpha" in chain[i] or "Alpha" in chain[i] or "ALPHA" in chain[i]:
                MHCa = i
            elif "microglobulin" in chain[i] or "MICROGLOBULIN" in chain[i] or "Microglobulin" in chain[i]:
                MHCb = i
        elif MHCtype == 2:
            # print(chain[i])
            if "alpha" in chain[i]:
                MHCa = i
            elif "Alpha" in chain[i]:
                MHCa = i
            elif "beta" in chain[i]:
                MHCb = i
            elif "Beta" in chain[i]:
                MHCb = i
            elif "ALPHA" in chain[i]:
                MHCa = i
            elif "BETA" in chain[i]:
                MHCb = i
            elif "DRa" in chain[i]:
                MHCa = i
            elif "DQa" in chain[i]:
                MHCa = i
            elif "DPa" in chain[i]:
                MHCa = i
            elif "DRB" in chain[i]:
                MHCb = i
            elif "DQB" in chain[i]:
                MHCb = i
            elif "DPB" in chain[i]:
                MHCb = i
            elif "H2-Ab" in chain[i]:
                MHCb = i
            elif "H2-Eb" in chain[i]:
                MHCb = i
            elif "B chain" in chain[i]:
                MHCb = i
    # print(MHCa)
    # print(MHCb)
    # print(chain[MHClist[1]])
    ## print(chain)
    ## print(len(chain))
    if TRAchain != None:
        TRAchain = chain[TRAchain][0]
    else:
        print("\nUnable to detect TCR alpha chain from FASTA data")
        print("\n", seq)
        TRAchain = input("Manually input the chain identifier for the TCR alpha sequence: ")
        TRAchain = str(TRAchain)
    if TRBchain != None:
        TRBchain = chain[TRBchain][0]
    else:
        print("\nUnable to detect TCR beta chain from FASTA data")
        print("\n", seq)
        TRBchain = input("Manually input the chain identifier for the TCR beta sequence: ")
        TRBchain = str(TRBchain)
    ##print("TCR alpha chain: ", TRAchain)
    ##print("TCR beta chain: ", TRBchain)
    if MHCa != None:
        MHCa = chain[MHCa][0]
    else:
        print("\nUnable to detect MHC/HLA alpha chain from FASTA data")
        print("\n", seq)
        MHCa = input("Manually input the chain identifier for the MHC/HLA alpha sequence: ")
        MHCa = str(MHCa)
    if MHCb != None:
        MHCb = chain[MHCb][0]
    else:
        print("\nUnable to detect MHC beta chain or B2M from FASTA data")
        print("\n", seq)
        MHCb = input("Manually input the chain identifier for the MHC beta chain or B2M sequence: ")
        MHCb = str(MHCb)
    peptide = None
    index = 0
    for i in chain:
        if "peptide" in i:
            peptide = index
        elif "Peptide" in i:
            peptide = index
        elif "PEPTIDE" in i:
            peptide = index
        elif "mimotope" in i:
            peptide = index
        elif "Mimotope" in i:
            peptide = index
        elif "MIMOTOPE" in i:
            peptide = index
        index = index + 1
    if peptide == None:
        print("\nUnable to detect peptide chain identifier from FASTA data")
        print("\n", seq)
        peptide = input("Manually input the chain identifier for the peptide chain\nIf not specified in FASTA descriptors enter \"None\": ")
        peptide = str(peptide)
    else:
        peptide = chain[peptide][0]

if peptide in [TRAchain, TRBchain]:
    peptide = "Error"
    print("Peptide is not a separate chain within PDB file")
    print("This may cause errors with caclulation of contact residues")
elif peptide in [MHCa, MHCb]:
    peptide = None

##print("MHCa chain: ", MHCa)
##if MHCtype == 1:
##    print("B2M chain : ", MHCb)
##elif MHCtype == 2:
##    print("MHCb chain: ", MHCb)

###########################################
#### Obtain TRX FASTAs from Pymol file ####
###########################################
        
TRA_FASTA = cmd.get_fastastr(selection=PDB + '_' + TRAchain, state=-1, quiet=1)
TRB_FASTA = cmd.get_fastastr(selection=PDB + '_' + TRBchain, state=-1, quiet=1)
# print(TRA_FASTA)
# print(TRB_FASTA)

####################################################
####             Identify Species               ####
####################################################

if "mouse" in seq or "mus mus" in seq or "Mus mus" in seq or "MUS MUS" in seq:
    species = "mouse"
elif "homo sapien" in seq or "human" in seq or "Homo" in seq or "HOMO" in seq:
    species = "human"
else:
    print("\nUnable to identify species")
    print("\n", seq)
    species = input("Please enter species (human or mouse): ")
    species = str(species)

########################################################################
#### Start of webscraping to obtain CDR sequences and TRX sequences ####
########################################################################

# initialize an HTTP session
session = HTMLSession()

def get_all_forms(url):
    """Returns all form tags found on a web page's `url` """
    # GET request
    res = session.get(url)
    # for javascript driven website
    # res.html.render()
    soup = BeautifulSoup(res.html.html, "html.parser")
    return soup.find_all("form")


url = "http://www.imgt.org/3Dstructure-DB/cgi/DomainGapAlign.cgi"

# get all form tags
forms = get_all_forms(url)
# print(forms)
# get just the first form
first_form = get_all_forms(url)[0]
action = first_form.attrs.get("action").lower()

if species == "mouse":
    species = "Mus musculus (house mouse)"
elif species == "human":
    species = "Homo sapiens (human)"
data = {"displayed": 1,
        "Sequence": TRA_FASTA,
        "domtype": "V",
        "species": species
        }
# data = urllib.parse.urlencode(data).encode("utf-8")

# join the url with the action (form request URL)
url = urljoin(url, action)
res = session.post(url, data=data)
#print(res.text)

text = res.text
soup2 = BeautifulSoup(text, 'html.parser')

# Extract CDR loop sequences
loops = soup2.find_all("span", class_="loop", limit=3)
loops = loops[0:3]
# print(loops)
a1 = loops[0]
a2 = loops[1]
a3 = loops[2]


def stripLoop(CDR):
    CDR = str(CDR)
    CDR = CDR.lstrip("<span class=\"loop\">")
    if "glycos" in CDR:
        CDR = CDR.replace("<span class=\"N-glycosylation\">", "")
    CDR = CDR.replace("</span>", "")
    CDR = CDR.replace(".", "")
    return (CDR)


a1 = stripLoop(a1)
a2 = stripLoop(a2)
a3 = stripLoop(a3)
##print("Alpha CDR sequences: ", a1, a2, a3)

# Extract TRA sequence
TRAV = soup2.find_all("span", title="V-REGION")
# print(TRAV)
if "glycos" in TRAV:
    TRAV = TRAV.replace("<span class=\"N-glycosylation\">", "")
TRAJ = soup2.find_all("span", title="J-REGION")
if "glycos" in TRAJ:
    TRAJ = TRAJ.replace("<span class=\"N-glycosylation\">", "")
TRAV = str(TRAV)
TRAV = TRAV.split("V-REGION\">")[1].split("</span>")[0]
TRAV = TRAV.replace("\n", "")
TRAJ = str(TRAJ)
TRAJ = TRAJ.split("J-REGION\">")[1].split("</span>")[0]
TRAJ = TRAJ.replace("\n", "")

if "(N-D)-REGION\">" in text:
    print("N additions in TRA")
    TRAnAdd = soup2.find_all("span", title="(N-D)-REGION")
    if "glycos" in TRAnAdd:
        TRAnAdd = TRAnAdd.replace("<span class=\"N-glycosylation\">", "")
    TRAnAdd = str(TRAnAdd)
    TRAnAdd = TRAnAdd.split("(N-D)-REGION\">")[1].split("</span>")[0]
    TRAnAdd = TRAnAdd.replace("\n", "")
    TRA = TRAV + TRAnAdd + TRAJ
else:
    print("No N additions in TRA")
    TRA = TRAV + TRAJ
##print("TRA: ", TRA)

#####################################################################
###### Now form upload to IMGT Domain Gap Align for Beta Chain ######
#####################################################################
data = {"displayed": 1,
        "Sequence": TRB_FASTA,
        "domtype": "V",
        "species": species
        }
# join the url with the action (form request URL)
url = urljoin(url, action)
res = session.post(url, data=data)
text = res.text
# print(text)
soup2 = BeautifulSoup(text, 'html.parser')
loops = soup2.find_all("span", class_="loop", limit=3)
loops = loops[0:3]
# print(loops)
b1 = loops[0]
b2 = loops[1]
b3 = loops[2]
b1 = stripLoop(b1)
b2 = stripLoop(b2)
b3 = stripLoop(b3)
##print("Beta CDR sequences: ", b1, b2, b3)

# Extract TRB sequence
TRBV = soup2.find_all("span", title="V-REGION")
if "glycos" in TRBV:
    TRBV = TRBV.replace("<span class=\"N-glycosylation\">", "")

TRBJ = soup2.find_all("span", title="J-REGION")
if "glycos" in TRBJ:
    TRBJ = TRBJ.replace("<span class=\"N-glycosylation\">", "")
TRBV = str(TRBV)
TRBV = TRBV.split("V-REGION\">")[1].split("</span>")[0]
TRBV = TRBV.replace("\n", "")
TRBJ = str(TRBJ)
TRBJ = TRBJ.split("J-REGION\">")[1].split("</span>")[0]
TRBJ = TRBJ.replace("\n", "")

if "(N-D)-REGION\">" in text:
    print("N-D region detected in TRB")
    TRBD = soup2.find_all("span", title="(N-D)-REGION")
    if "glycos" in TRBD:
        TRBD = TRBD.replace("<span class=\"N-glycosylation\">", "")
    TRBD = str(TRBD)
    TRBD = TRBD.split("(N-D)-REGION\">")[1].split("</span>")[0]
    TRBD = TRBD.replace("\n", "")
    TRB = TRBV + TRBD + TRBJ
else:
    print("No N-D region detected in TRB")
    TRB = TRBV + TRBJ
    
##print("TRB: ", TRB)

print("Copy and Paste the following information into the inputs section of the da.py file:")
print("PDB = \"%s\"" % PDB)
print("MHCtype = %d" % MHCtype)
print("a1 = \"%s\"" % a1)
print("a2 = \"%s\"" % a2)
print("a3 = \"%s\"" % a3)
print("b1 = \"%s\"" % b1)
print("b2 = \"%s\"" % b2)
print("b3 = \"%s\"" % b3)
print("TRAchain = \"%s\"" % TRAchain)
print("TRBchain = \"%s\"" % TRBchain)
print("MHCa = \"%s\"" % MHCa)
print("MHCb = \"%s\"" % MHCb)
print("peptide = \"%s\"" % peptide)
print("TRA = \"%s\"" % TRA)
print("TRB = \"%s\"" % TRB)
print("contact_dist = %d" % contact_dist)
