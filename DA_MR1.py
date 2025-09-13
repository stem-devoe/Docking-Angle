#####################################################################################
####
#### DISCLAIMER: IF THE PEPTIDE SEQUENCE IS COVALENTLY LINKED TO THE TCR BETA    ####
#### CHAIN, THE CODE WILL NOT PROPERLY FUNCTION. I think I managed to catch most ####
#### of the instances where this occurs in the genDAinput script, but if you know####
#### this is the case for your structure, double check it please!                ####
####
#####################################################################################



######## Paste Output from genDAinput.py Here: #############


PDB = "4LCC"
MHCtype = 1
a1 = "TSGFNG"
a2 = "NVLDGL"
a3 = "AVKDSNYQLI"
b1 = "MNHNS"
b2 = "SASEGT"
b3 = "ASSVWTGEGSGELF"
TRAchain = "A"
TRBchain = "B"
MHCa = "C"
peptide = "None"
TRA = "GQNIDQPTEMTATEGAIVQINCTYQTSGFNGLFWYQQHAGEAPTFLSYNVLDGLEEKGRFSSFLSRSKGYSYLLLKELQMKDSASYLCAVKDSNYQLIWGAGTKLIIKP"
TRB = "NAGVTQTPKFQVLKTGQSMTLQCAQDMNHNSMYWYRQDPGMGLRLIYYSASEGTTDKGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASSVWTGEGSGELFFGEGSRLTVL"


######## End data input ##############

contact_dist = [3.5, 4.5]

#####################################################################################
#####           Do Not change code below this point!                            #####
#####################################################################################
import pymol
from pymol import cmd
from Bio import pairwise2
from Bio.Align import substitution_matrices
import numpy
from skspatial.objects import Points, Line
from itertools import groupby, count

######### plane functions ########
'''
Described at PyMOL wiki:
https://pymolwiki.org/index.php/Plane_Wizard

Authors : Troels Schwarz-Linnet
Date    : Dec 2016
Modified: From previous contributors.
'''


from pymol.wizard import Wizard
from chempy import cpv
from pymol.cgo import COLOR, SPHERE, CYLINDER, BEGIN, TRIANGLE_STRIP, NORMAL, VERTEX, END, ALPHA

def makePrimitive(cgo, name):
    az = cmd.get('auto_zoom', quiet=1)
    cmd.set('auto_zoom', 0, quiet=1)
    cmd.load_cgo(cgo, name)
    cmd.set('auto_zoom', az, quiet=1)

def point(p):
    x, y, z = p
    return [COLOR, 1, 1, 1, SPHERE, float(x), float(y), float(z), 0.5]


def line(p1, p2):
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    return [CYLINDER, float(x1), float(y1), float(z1), float(x2), float(y2), float(z2), 0.25, 1, 1, 1, 1, 1, 1]


def plane(corner1, corner2, corner3, corner4, normal, settings):
    planeObj = []
    planeObj.extend(point(corner1))
    planeObj.extend(point(corner2))
    planeObj.extend(point(corner3))
    planeObj.extend(point(corner4))
    planeObj.extend(line(corner1, corner2))
    planeObj.extend(line(corner2, corner3))
    planeObj.extend(line(corner3, corner4))
    planeObj.extend(line(corner4, corner1))

    # Make settings
    if 'ALPHA' in settings:
        planeObj.extend([ALPHA, settings['ALPHA']])

    if 'COLOR' in settings:
        planeObj.extend([COLOR, settings['COLOR'][0], settings['COLOR'][1], settings['COLOR'][2]])
    else:
        planeObj.extend([COLOR, 0.8, 0.8, 0.8]) # greyish

    planeObj.extend([BEGIN, TRIANGLE_STRIP])
    planeObj.append(NORMAL)

    if 'INVERT' in settings:
        if settings['INVERT']==True:
            planeObj.extend(cpv.negate(normal))
        else:
            planeObj.extend(normal)
    else:
        planeObj.extend(normal)


    for corner in [corner1, corner2, corner3, corner4, corner1]:
        planeObj.append(VERTEX)
        planeObj.extend(corner)
    planeObj.append(END)
    return planeObj


def planeFromPoints(p1, p2, p3, vm1=None, vm2=None, center=True, settings={}):
    v1 = cpv.sub(p1, p2)
    v2 = cpv.sub(p3, p2)
    normal = cpv.cross_product(v1, v2)

    if 'translate' in settings:
        vtran = cpv.scale(cpv.normalize(normal), settings['translate'])
        p1_t = cpv.sub(p1, vtran)
        p2_t = cpv.sub(p2, vtran)
        p3_t = cpv.sub(p3, vtran)
        print("New coordinates are:")
        print_info("New", p1_t, p2_t, p3_t)
        print("New coordinates are for normalized plane:")
        v1_t = cpv.normalize(cpv.sub(p1_t, p2_t))
        v2_t = cpv.normalize(cpv.sub(p3_t, p2_t))
        normal_t = cpv.normalize(cpv.cross_product(v1_t, v2_t))
        v2_t = cpv.normalize(cpv.cross_product(normal_t, v1_t))
        p1_t2 = cpv.add(v1_t, p2_t)
        p3_t2 = cpv.add(v2_t, p2_t)
        print_info("Newnormal", p1_t2, p2_t, p3_t2)

    if vm1!=None:
        v1 = cpv.scale(cpv.normalize(v1), vm1)
    if vm2!=None:
        v2 = cpv.scale(cpv.normalize(v2), vm2)

    centrum = p2
    if center:
        corner1 = cpv.add(cpv.add(centrum, v1), v2)
        corner2 = cpv.sub(cpv.add(centrum, v1), v2)
        corner3 = cpv.sub(cpv.sub(centrum, v1), v2)
        corner4 = cpv.add(cpv.sub(centrum, v1), v2)
    else:
        corner1 = cpv.add(cpv.add(centrum, v1), v2)
        corner2 = cpv.add(centrum, v1)
        corner3 = centrum
        corner4 = cpv.add(centrum, v2)

    return plane(corner1, corner2, corner3, corner4, normal, settings)


def print_info(name, coor1, coor2, coor3):
    cs1 = (map(float, [ '%.2f' % elem for elem in coor1 ]) )
    cs2 = (map(float, [ '%.2f' % elem for elem in coor2 ]) )
    cs3 = (map(float, [ '%.2f' % elem for elem in coor3 ]) )
    print("You can also use the function calls with these coordinates")
    print("plane.make_plane_points(name='%s', l1=%s, l2=%s, l3=%s)"%(name, cs1, cs2, cs3))


def make_plane(name,a1='(pk1)',a2='(pk2)',a3='(pk3)', vm1=None, vm2=None, center=True, makepseudo=True, settings={}):
    """
    DESCRIPTION
    Create a CGO plane from three atomic coordinates

    USAGE
    make_plane name, a1, a2, a3

    where each atom is a standard PyMOL selection (defaults to pk1,pk2 and pk3)
    """
    # get coordinates from atom selections
    coor1 = cmd.get_model(a1).get_coord_list()[0]
    coor2 = cmd.get_model(a2).get_coord_list()[0]
    coor3 = cmd.get_model(a3).get_coord_list()[0]

    # Help with alternative
    print_info(name, coor1, coor2, coor3)

    # Get the plane
    plane = planeFromPoints(p1=coor1, p2=coor2, p3=coor3, vm1=vm1, vm2=vm2, center=center, settings=settings)
    makePrimitive(plane, name)
    #cmd.show("cgo", "plane*")

    if makepseudo:
        cmd.pseudoatom("%s_%s"%(name, "l1"), color="tv_blue", pos=coor1)
        cmd.pseudoatom("%s_%s"%(name, "l2"), color="tv_green", pos=coor2)
        cmd.pseudoatom("%s_%s"%(name, "l3"), color="red", pos=coor3)

# Extend function to be called inside pymol
#cmd.extend("make_plane", make_plane)

def make_plane_points(name,l1=None,l2=None,l3=None, vm1=None, vm2=None, center=True, makepseudo=True, settings={}):
    """
    DESCRIPTION
    Create a CGO plane from three atomic coordinates

    USAGE
    make_plane name, l1, l2, l3

    where each xys is a list with floats of x,y,z coordinates
    """
    if l1==None or l2==None or l3==None:
        print("Please provide a list of xyz floats for each 3 positions")
        return
    if type(l1) is not list or type(l2) is not list or type(l3) is not list:
        print(type(l1),type(l2),type(l3))
        print("Please provide 3 list of xyz floats for each 3 positions")
        return

    plane = planeFromPoints(p1=l1, p2=l2, p3=l3, vm1=vm1, vm2=vm2, center=center, settings=settings)
    makePrimitive(plane, name)

    if makepseudo:
        cmd.pseudoatom("%s_%s"%(name, "l1"), color="tv_blue", pos=l1)
        cmd.pseudoatom("%s_%s"%(name, "l2"), color="tv_green", pos=l2)
        cmd.pseudoatom("%s_%s"%(name, "l3"), color="red", pos=l3)

# Extend function to be called inside pymol
#cmd.extend("make_plane_points", make_plane_points)

class PlaneWizard(Wizard):
    def __init__(self):
        Wizard.__init__(self)

        # some attributes to do with picking
        self.pick_count = 0
        self.object_count = 0
        self.object_prefix = "pw"

        self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
        cmd.set("mouse_selection_mode",0) # set selection mode to atomic
        cmd.deselect()

    def reset(self):
        cmd.delete(self.object_prefix + "*")
        cmd.delete("sele*")
        cmd.delete("_indicate*")
        cmd.unpick()
        cmd.refresh_wizard()

    def delete_all(self):
        cmd.delete("plane*")

    def cleanup(self):
        cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode
        self.reset()
        self.delete_all()

    def get_prompt(self):
        self.prompt = None
        if self.pick_count == 0:
            self.prompt = [ 'Please click on the first atom...']
        elif self.pick_count == 1:
            self.prompt = [ 'Please click on the second atom...' ]
        elif self.pick_count == 2:
            self.prompt = [ 'Please click on the third atom...' ]
        return self.prompt

    def do_select(self, name):
        # "edit" only this atom, and not others with the object prefix
        try:
            cmd.edit("%s and not %s*" % (name, self.object_prefix))
            self.do_pick(0)
        except(pymol.CmdException, pmce):
            print(pmce)

    def pickNextAtom(self, atom_name):
        # transfer the click selection to a named selection
        cmd.select(atom_name, "(pk1)")

        # delete the click selection
        cmd.unpick()

        # using the magic of indicate, highlight stuff
        indicate_selection = "_indicate" + self.object_prefix
        cmd.select(indicate_selection, atom_name)
        cmd.enable(indicate_selection)

        self.pick_count += 1
        self.error = None

        # necessary to force update of the prompt
        cmd.refresh_wizard()

    def do_pick(self, picked_bond):

        # this shouldn't actually happen if going through the "do_select"
        if picked_bond:
            self.error = "Error: please select bonds, not atoms"
            print(self.error)
            return

        atom_name = self.object_prefix + str(self.pick_count)
        if self.pick_count < 2:
            self.pickNextAtom(atom_name)
        else:
            self.pickNextAtom(atom_name)

            point1 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "0"))
            point2 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "1"))
            point3 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "2"))
            plane = planeFromPoints(point1, point2, point3)

            planeName = "plane-%02d" % self.object_count

            print_info(planeName, point1, point2, point3)

            self.object_count += 1
            makePrimitive(plane, planeName)
            cmd.show("cgo", "plane*")

            self.pick_count = 0
            self.reset()

    def get_panel(self):
        return [
            [ 1, 'Plane Wizard',''],
            [ 2, 'Reset','cmd.get_wizard().reset()'],
            [ 2, 'Delete All Planes' , 'cmd.get_wizard().delete_all()'],
            [ 2, 'Done','cmd.set_wizard()'],
        ]

######## end plane functions ########

################## Start of Stephanie's Functions ###############
def getNorm(com):
    Ax = com[0][0] - com[1][0]
    Ay = com[0][1] - com[1][1]
    Az = com[0][2] - com[1][2]
    Bx = com[0][0] - com[2][0]
    By = com[0][1] - com[2][1]
    Bz = com[0][2] - com[2][2]
    Norm = numpy.cross([Ax,Ay,Az], [Bx,By,Bz])
    return(Norm)

def findAlphaPoint(tneg, tpos, ref):
    tneg1 = numpy.array(tneg)
    tpos1 = numpy.array(tpos)
    d_tneg = numpy.linalg.norm(tneg1 - ref)
    d_tpos = numpy.linalg.norm(tpos1 - ref)
    if(d_tneg < d_tpos):
        ab = [tneg, tpos]
    elif(d_tneg == d_tpos):
        print("Error. points are equidistant to TRA COM")
        ab = []
    else:
        ab = [tpos, tneg]
    #print(ab)
    return(ab)

def getDirPlane(tneg, tpos, anchor, ref, pointsBool):
    ab = findAlphaPoint(tneg, tpos, ref)
    l = [anchor, ab[0], ab[1]]
    norm = getNorm(l)
    if(pointsBool):
        cmd.pseudoatom('adirpoint', pos = ab[0])
        #print(ab[0])
        cmd.pseudoatom('bdirpoint', pos=ab[1])
    return(norm)

def getAngle(A, B, Model):
    d_prod = numpy.dot(A, B)
    magA = getMag(A)
    magB = getMag(B)
    if Model == "PlanePlane" or Model == "VectorVector":
        Angle = numpy.arccos(d_prod / (magA * magB)) * 180 / numpy.pi
    elif Model == "PlaneVector":
        Angle = numpy.arcsin(d_prod / (magA * magB)) * 180 / numpy.pi
    return(Angle)

def getMag(Vector):
    Mag = numpy.power(Vector[0], 2) + numpy.power(Vector[1], 2) + numpy.power(Vector[2], 2)
    Mag = numpy.power(Mag, 0.5)
    return(Mag)


def getPlane(Norm, PlanePoint):
    d = Norm[0] * PlanePoint[0] + Norm[1] * PlanePoint[1] + Norm[2] * PlanePoint[2]
    Plane = [Norm[0], Norm[1], Norm[2], d]
    return (Plane)


def solveIntersect(Plane, VecDir, VecPoint):
    t = (Plane[3] - Plane[0] * VecPoint[0] - Plane[1] * VecPoint[1] - Plane[2] * VecPoint[2]) / (
                Plane[0] * VecDir[0] + Plane[1] * VecDir[1] + Plane[2] * VecDir[2])
    inter = [VecDir[0] * t + VecPoint[0], VecDir[1] * t + VecPoint[1], VecDir[2] * t + VecPoint[2]]
    return (inter)


def getDist(Intersection, Reference):
    x = numpy.array(Intersection)
    y = numpy.array(Reference)
    dist = numpy.linalg.norm(x - y)
    return (dist)

####### end of stephanie's functions #######


cmd.reinitialize()

# load PDB structure
cmd.fetch(PDB)

# split by chain and color by chain
cmd.split_chains(PDB)
cmd.color('cyan', PDB + '_' + MHCa)
#cmd.color('green', PDB + '_' + MHCb)
cmd.color('yellow', PDB + '_' + TRAchain)
cmd.color('orange', PDB + '_' + TRBchain)

cmd.delete(PDB)
cmd.hide('all')

# obtain points of contact
for i in contact_dist:
    if peptide == "Error":
        print("The peptide is part of another chain within the PDB file.\nContact Residues cannot be determined automatically")
    elif peptide != "None":
        #print("testcheckpoint")
        cmd.select('AlphaContact' + "_" + str(i), '%s near_to %f of (%s or %s)' % (PDB + '_' + TRAchain, i, PDB + '_' + MHCa, PDB + '_' + peptide))
        cmd.select('BetaContact'+ "_" + str(i), '%s near_to %f of (%s or %s)' % (PDB + '_' + TRBchain, i, PDB + '_' + MHCa, PDB + '_' + peptide))
        cmd.select('AlphaContactCa' + "_" + str(i), 'bca. %s' % 'AlphaContact'+ "_" + str(i))
        cmd.select('BetaContactCa'+ "_" + str(i), 'bca. %s' % 'BetaContact'+ "_" + str(i))
    else:
        cmd.select('AlphaContact' + "_" + str(i), '%s near_to %f of ((%s) and not %s)' % (PDB + '_' + TRAchain, i, PDB + '_' + MHCa, PDB + '_' + TRBchain))
        cmd.select('BetaContact' + "_" + str(i), '%s near_to %f of ((%s) and not %s)' % (PDB + '_' + TRBchain, i, PDB + '_' + MHCa, PDB + '_' + TRAchain))
        cmd.select('AlphaContactCa' + "_" + str(i), 'bca. %s' % 'AlphaContact' + "_" + str(i))
        cmd.select('BetaContactCa' + "_" + str(i), 'bca. %s' % 'BetaContact' + "_" + str(i))

# obtain centerofmass and pseudo atom markers for CDRs
CDR_list = [a1, a2, a3, b1, b2, b3]
selNames = ["CDR1a", "CDR2a", "CDR3a", "CDR1b", "CDR2b", "CDR3b"]
CDR_COM_list = [0,0,0,0,0,0]
find_list = [0,0,0,0,0,0]
index = 0
for i in CDR_list:
    print(i)
    find_list[index] = cmd.select(selNames[index], 'pepseq %s in (%s or %s)' % (CDR_list[index], PDB + "_" + TRAchain, PDB + "_" + TRBchain))
    CDR_COM_list[index] = cmd.centerofmass(selection = selNames[index])
    index = index + 1


# obtain centerofmass and pseudo atom markers for TRX
cmd.select('TRA', 'pepseq %s in %s' % (TRA, PDB + "_" + TRAchain))
TRA_COM = cmd.centerofmass('TRA')
###cmd.pseudoatom('TRA_COM', pos = TRA_COM)
cmd.select('TRB', 'pepseq %s in %s' % (TRB, PDB + "_" + TRBchain))
TRB_COM = cmd.centerofmass('TRB')
###cmd.pseudoatom('TRB_COM', pos = TRB_COM)


# merged CDR and TRX planes
cmd.select('allCDR', '((pepseq %s pepseq %s pepseq %s pepseq %s pepseq %s pepseq %s) in (%s or %s))' % (a1, a2, a3, b1, b2, b3, PDB + "_" + TRAchain, PDB + "_" + TRBchain))
allCDR_COM = cmd.centerofmass('allCDR')
###cmd.pseudoatom('allCDR_COM', pos = allCDR_COM)
####make_plane_points('plane_allCDR_TRXs', allCDR_COM, TRA_COM, TRB_COM)

cmd.select('allTRX', '((pepseq %s pepseq %s) in (%s or %s))' % (TRA, TRB, PDB + "_" + TRAchain, PDB + "_" + TRBchain))
allTRX_COM = cmd.centerofmass('allTRX')
###cmd.pseudoatom('allTRX_COM', pos = allTRX_COM)



############### Binding Groove Fitting Using Ca #############
# get atom coordinates
# MR1 a2 helix references
hMR1_2 = {'sequence': 'GSHTYQRMIGCELLEDGSTTGFLQYAYDGQDFLIFNKDTLSWLAVDNVAHTIKQAWEANQHELLYQKNWLEEECIAWLKRFLEYGKDTLQRT',
          'BG': 'NVAHTIKQAWEANQHELLYQKNWLEEECIAWLKRFLEYGKDTLQ',
          # IMGT helix defn NVAHTIKQAWEANQHELLYQKNWLEEECIAWLKRFLEYGKDTLQRT
          'score': 0, 'structSeq': ''}
mMR1_2 = {'sequence': 'GLHTYQRMIGCELLEDGSTTGFLQYAYDGQDFIIFNKDTLSWLAMDYVAHITKQAWEANLHELQYQKNWLEEECIAWLKRFLEYGRDTLERT',
          'BG': 'YVAHITKQAWEANLHELQYQKNWLEEECIAWLKRFLEYGRDTLE',
          # IMGT helix defn YVAHITKQAWEANLHELQYQKNWLEEECIAWLKRFLEYGRDTLERT
          'score': 0, 'structSeq': ''}
Beta = [hMR1_2, mMR1_2]
# MR1 a1 helix references
hMR1 = {'sequence': 'RTHSLRYFRLGVSDPIHGVPEFISVGYVDSHPITTYDSVTRQKEPRAPWMAENLAPDHWERYTQLLRGWQQMFKVELKRLQRHYNHS',
        'BG': 'PWMAENLAPDHWERYTQLLRGWQQMFKVELKRLQRHYN',  # IMGT helix defn PWMAENLAPDHWERYTQLLRGWQQMFKVELKRLQRHYNHS
        'score': 0, 'structSeq': ''}
mMR1 = {'sequence': 'RTHSLRYFRLAVSDPGPVVPEFISVGYVDSHPITTYDSVTRQKEPKAPWMAENLAPDHWERYTQLLRGWQQTFKAELRHLQRHYNHS',
        'BG': 'PWMAENLAPDHWERYTQLLRGWQQTFKAELRHLQRHYN',  # IMGT helix defn PWMAENLAPDHWERYTQLLRGWQQTFKAELRHLQRHYNHS
        'score': 0, 'structSeq': ''}
Alpha = [hMR1, mMR1]

# set alignment parameters
matrix = substitution_matrices.load(name="BLOSUM62")
gap_open = -10
gap_extend = -0.5

## get MHC sequences in PDB Structure ##
MHCSeq = cmd.get_fastastr(selection=PDB + '_' + MHCa, state=-1, quiet=1)
# print(MHCSeq)
MHCSeq = MHCSeq.split("\n")
# print(MHCSeq)
MHCSeq = ''.join([MHCSeq[i] for i in range(len(MHCSeq)) if i != 0])
# print("\n")
# print(MHCSeq)
if "?" in MHCSeq:
    MHCSeq = MHCSeq.replace("?", "X")

## Align A1 helix
scores = []
for i in Alpha:
    algn = pairwise2.align.globalds(i['sequence'], MHCSeq, matrix, gap_open, gap_extend, one_alignment_only=True)
    # print(algn)
    i['sequence'], i['structSeq'], i['score'], start, end = algn[0]
    # print(i['score'])
    scores.append(i['score'])
print(scores)
indA = scores.index(max(scores))
# print(Alpha[indA]['structSeq'])

# determine BG1 sequence
startA = 0
hold = Alpha[indA]['sequence']
BGbool = 1
while (BGbool == 1):
    if (Alpha[indA]['BG'] in hold.replace("-", "")):
        hold = hold[1:len(hold)]
        # print(startA)
        # print(hold)
        startA = startA + 1
    else:
        BGbool = 0
        # print(startA)
        # print(Alpha[indA]['sequence'][(startA-1):len(Alpha[indA]['sequence'])])

endA = len(Alpha[indA]['sequence'])
hold = Alpha[indA]['sequence']
BGbool = 1
while (BGbool == 1):
    if (Alpha[indA]['BG'] in hold.replace("-", "")):
        hold = hold[0:len(hold) - 1]
        # print(endA)
        # print(hold)
        endA = endA - 1
    else:
        BGbool = 0
        # print(endA)
        # print(Alpha[indA]['sequence'][0:(endA + 1)])

# print(Alpha[indA]['sequence'][(startA-1):(endA+1)])
# print(Alpha[indA]['structSeq'][(startA-1):(endA+1)])
BG1 = Alpha[indA]['structSeq'][(startA - 1):(endA + 1)]
# print(BG1)
BG1 = BG1.replace("-", "")
print(BG1)

## Align A2 helix
scores = []
for i in Beta:
    algn = pairwise2.align.globalds(i['sequence'], MHCSeq, matrix, gap_open, gap_extend, one_alignment_only=True)
    # print(algn)
    i['sequence'], i['structSeq'], i['score'], start, end = algn[0]
    # print(i['score'])
    scores.append(i['score'])
print(scores)
indB = scores.index(max(scores))
# print(Alpha[indA]['structSeq'])

startB = 0
hold = Beta[indB]['sequence']
BGbool = 1
while (BGbool == 1):
    if (Beta[indB]['BG'] in hold.replace("-", "")):
        hold = hold[1:len(hold)]
        # print(startB)
        # print(hold)
        startB = startB + 1
    else:
        BGbool = 0
        # print(startB)
        # print(Beta[indB]['sequence'][(startB-1):len(Beta[indB]['sequence'])])

endB = len(Beta[indB]['sequence'])
hold = Beta[indB]['sequence']
BGbool = 1
while (BGbool == 1):
    if (Beta[indB]['BG'] in hold.replace("-", "")):
        hold = hold[0:len(hold) - 1]
        # print(endB)
        # print(hold)
        endB = endB - 1
    else:
        BGbool = 0
        # print(endB)
        # print(Beta[indB]['sequence'][0:(endB + 1)])

# print(Beta[indB]['sequence'][(startB-1):(endB+1)])
# print(Beta[indB]['structSeq'][(startB-1):(endB+1)])
BG2 = Beta[indB]['structSeq'][(startB - 1):(endB + 1)]
# print(BG2)
BG2 = BG2.replace("-", "")
print(BG2)

## select binding groove residues
cmd.select('BGheli', '(pepseq %s in %s) or (pepseq %s in %s)' % (BG1, PDB + "_" + MHCa, BG2, PDB + "_" + MHCa))
cmd.color('salmon', 'BGheli')
cmd.select('BGheli_Ca', 'bca. BGheli')

MHC_helix_Ca_coords = []
cmd.iterate_state(1, 'BGheli_Ca', 'MHC_helix_Ca_coords.append([x,y,z])')

# fit line in 3D space for binding groove
points = Points(MHC_helix_Ca_coords)
line_fit = Line.best_fit(points)

# plot BGvector
t = [-20,-10,10,20]
linePoints = [0,0,0,0]
holder = []
index = 0
for i in t:
    holder = line_fit.point + i * line_fit.direction
    holder = [holder[0],holder[1],holder[2]]
    linePoints[index] = holder
    index = index + 1
BGpoint = [line_fit.point[0], line_fit.point[1], line_fit.point[2]]
cmd.pseudoatom('t0', pos = BGpoint)
cmd.show_as('spheres', 't0')
cmd.color('green', 't0')

# calculate angles between BGvector and TCR planes
cmd.select('firstMHC_Ca', 'first BGheli_Ca')

firstMHC = cmd.centerofmass('firstMHC_Ca')
d_BGneg = numpy.linalg.norm(numpy.array(firstMHC) - numpy.array(linePoints[0]))
d_BGpos = numpy.linalg.norm(numpy.array(firstMHC) - numpy.array(linePoints[3]))
if(d_BGneg > d_BGpos):
    BGvec = [-line_fit.direction[0], -line_fit.direction[1], -line_fit.direction[2]]
    index = 0
    for i in linePoints:
        if t[index] < 0:
            name = 't' + 'neg' + str(abs(t[index]))
        else:
            name = 't' + str(t[index])
        cmd.pseudoatom(name, pos=i)
        if t[index] < 0:
            cmd.show_as('spheres', 't' + 'neg' + str(abs(t[index])))
            cmd.color('red', 't' + 'neg' + str(abs(t[index])))
        else:
            cmd.show_as('spheres', 't' + str(abs(t[index])))
            cmd.color('salmon', 't' + str(abs(t[index])))
        index = index + 1
else:
    BGvec = [line_fit.direction[0],line_fit.direction[1],line_fit.direction[2]]
    index = 0
    for i in linePoints:
        if t[index] < 0:
            name = 't' + 'neg' + str(abs(t[index]))
        else:
            name = 't' + str(t[index])
        cmd.pseudoatom(name, pos=i)
        if t[index] < 0:
            cmd.show_as('spheres', 't' + 'neg' + str(abs(t[index])))
            cmd.color('salmon', 't' + 'neg' + str(abs(t[index])))
        else:
            cmd.show_as('spheres', 't' + str(abs(t[index])))
            cmd.color('red', 't' + str(abs(t[index])))
        index = index + 1

###################################################################################################
############# fit line in 3D space for CDRs of Beta chain and Alpha chain  ########################
############# Using All Atoms of CDR Residues                              ########################
############# to model Alpha and Beta Planes                               ########################

CDRa_coords = []
cmd.iterate_state(1, 'CDR1a', 'CDRa_coords.append([x,y,z])')
cmd.iterate_state(1, 'CDR2a', 'CDRa_coords.append([x,y,z])')
cmd.iterate_state(1, 'CDR3a', 'CDRa_coords.append([x,y,z])')
CDRa_points = Points(CDRa_coords)
CDRa_line_fit = Line.best_fit(CDRa_points)

cdr = [-20,-10,10,20]
CDRa_linePoints = [0,0,0,0]
holder = []
index = 0
for i in cdr:
    holder = CDRa_line_fit.point + i * CDRa_line_fit.direction
    holder = [holder[0],holder[1],holder[2]]
    CDRa_linePoints[index] = holder
    index = index + 1

CDRa_fit = getDirPlane(CDRa_linePoints[0], CDRa_linePoints[3], TRA_COM, CDR_COM_list[1], pointsBool = False)
####make_plane_points('CDRa_fit', TRA_COM, CDRa_linePoints[0], CDRa_linePoints[3])
BG_CDRa_fitPlane = getAngle(BGvec, CDRa_fit, "PlaneVector")


CDRb_coords = []
cmd.iterate_state(1, 'CDR1b', 'CDRb_coords.append([x,y,z])')
cmd.iterate_state(1, 'CDR2b', 'CDRb_coords.append([x,y,z])')
cmd.iterate_state(1, 'CDR3b', 'CDRb_coords.append([x,y,z])')
CDRb_points = Points(CDRb_coords)
CDRb_line_fit = Line.best_fit(CDRb_points)

cdr = [-20,-10,10,20]
CDRb_linePoints = [0,0,0,0]
holder = []
index = 0
for i in cdr:
    holder = CDRb_line_fit.point + i * CDRb_line_fit.direction
    holder = [holder[0],holder[1],holder[2]]
    CDRb_linePoints[index] = holder
    index = index + 1

CDRb_fit = getDirPlane(CDRb_linePoints[0], CDRb_linePoints[3], TRB_COM, CDR_COM_list[4], pointsBool = False)
###make_plane_points('CDRb_fit', TRB_COM, CDRb_linePoints[0], CDRb_linePoints[3])
BG_CDRb_fitPlane = getAngle(BGvec, CDRb_fit, "PlaneVector")

#######################################################################################
########## GERMLINE (just CDR1 and CDR2)                                      #########
Agerm_coords = []
cmd.iterate_state(1, 'CDR1a', 'Agerm_coords.append([x,y,z])')
cmd.iterate_state(1, 'CDR2a', 'Agerm_coords.append([x,y,z])')
Agerm_points = Points(Agerm_coords)
Agerm_line_fit = Line.best_fit(Agerm_points)

Agerm = [-20,-10,10,20]
Agerm_linePoints = [0,0,0,0]
holder = []
index = 0
for i in Agerm:
    holder = Agerm_line_fit.point + i * Agerm_line_fit.direction
    holder = [holder[0],holder[1],holder[2]]
    Agerm_linePoints[index] = holder
    index = index + 1

Agerm_fit = getDirPlane(Agerm_linePoints[0], Agerm_linePoints[3], TRA_COM, CDR_COM_list[1], pointsBool = False)
###make_plane_points('Agerm_fit', TRA_COM, Agerm_linePoints[0], Agerm_linePoints[3])
BG_Agerm_fitPlane = getAngle(BGvec, Agerm_fit, "PlaneVector")


Bgerm_coords = []
cmd.iterate_state(1, 'CDR1b', 'Bgerm_coords.append([x,y,z])')
cmd.iterate_state(1, 'CDR2b', 'Bgerm_coords.append([x,y,z])')
Bgerm_points = Points(Bgerm_coords)
Bgerm_line_fit = Line.best_fit(Bgerm_points)

Bgerm = [-20,-10,10,20]
Bgerm_linePoints = [0,0,0,0]
holder = []
index = 0
for i in Bgerm:
    holder = Bgerm_line_fit.point + i * Bgerm_line_fit.direction
    holder = [holder[0],holder[1],holder[2]]
    Bgerm_linePoints[index] = holder
    index = index + 1

Bgerm_fit = getDirPlane(Bgerm_linePoints[0], Bgerm_linePoints[3], TRB_COM, CDR_COM_list[4], pointsBool = False)
###make_plane_points('Bgerm_fit', TRB_COM, Bgerm_linePoints[0], Bgerm_linePoints[3])
BG_Bgerm_fitPlane = getAngle(BGvec, Bgerm_fit, "PlaneVector")

TCRgerm_coords = []
cmd.iterate_state(1, 'CDR1a', 'TCRgerm_coords.append([x,y,z])')
cmd.iterate_state(1, 'CDR2a', 'TCRgerm_coords.append([x,y,z])')
cmd.iterate_state(1, 'CDR1b', 'TCRgerm_coords.append([x,y,z])')
cmd.iterate_state(1, 'CDR2b', 'TCRgerm_coords.append([x,y,z])')
TCRgerm_points = Points(TCRgerm_coords)
TCRgerm_line_fit = Line.best_fit(TCRgerm_points)

TCRgerm = [-20,-10,10,20]
TCRgerm_linePoints = [0,0,0,0]
holder = []
index = 0
for i in TCRgerm:
    holder = TCRgerm_line_fit.point + i * TCRgerm_line_fit.direction
    holder = [holder[0],holder[1],holder[2]]
    TCRgerm_linePoints[index] = holder
    index = index + 1

TCRgerm_fit = getDirPlane(TCRgerm_linePoints[0], TCRgerm_linePoints[3], allTRX_COM, TRA_COM, pointsBool = False)
###make_plane_points('TCRgerm_fit', allTRX_COM, TCRgerm_linePoints[0], TCRgerm_linePoints[3])
BG_TCRgerm_fitPlane = getAngle(BGvec, TCRgerm_fit, "PlaneVector")

######################################################################################

    # print results
print("PDB ID: %s" % PDB)
print('%d points were used in the fitting of the binding groove'% len(MHC_helix_Ca_coords))
###print("The angle between the Binding Groove vector model and the CDRa TRA plane model is:")
print("%f" % BG_CDRa_fitPlane)
###print("The angle between the Binding Groove vector model and the TRA Germline plane model is:")
print("%f" % BG_Agerm_fitPlane)
###print("The angle between the Binding Groove vector model and the CDRb TRB plane model is:")
print("%f" % BG_CDRb_fitPlane)
###print("The angle between the Binding Groove vector model and the TRB Germline plane model is:")
print("%f" % BG_Bgerm_fitPlane)


###################################################################################################
######################### fit line in 3D space for all CDRs Using Ca  #############################
### print("CDR Ca")

CDR_Ca_coords = []
cmd.select('allCDR_Ca', 'bca. %s' % 'allCDR')
cmd.iterate_state(1, 'allCDR_Ca', 'CDR_Ca_coords.append([x,y,z])')
CDR_Capoints = Points(CDR_Ca_coords)
CDR_Caline_fit = Line.best_fit(CDR_Capoints)

# plot CDRvector
cdr = [-20,-10,10,20]
CDR_CalinePoints = [0,0,0,0]
holder = []
index = 0
for i in cdr:
    holder = CDR_Caline_fit.point + i * CDR_Caline_fit.direction
    holder = [holder[0],holder[1],holder[2]]
    CDR_CalinePoints[index] = holder
    index = index + 1
#index = 0
#for i in CDR_CalinePoints:
#    if cdr[index] < 0:
#        name = 'cdr' + 'neg' + str(abs(cdr[index]))
#    else:
#        name = 'cdr' + str(cdr[index])
#    cmd.pseudoatom(name, pos = i)
#    if t[index] < 0:
#        cmd.show_as('spheres', 'cdr' + 'neg' + str(abs(cdr[index])))
#        cmd.color('magenta', 'cdr' + 'neg' + str(abs(cdr[index])))
#    else:
#        cmd.show_as('spheres', 'cdr' + str(abs(cdr[index])))
#        cmd.color('magenta', 'cdr' + str(abs(cdr[index])))
#    index = index + 1
#cmd.pseudoatom('cdr0', pos = [CDR_Caline_fit.point[0],CDR_Caline_fit.point[1],CDR_Caline_fit.point[2]])
#cmd.show_as('spheres', 'cdr0')
#cmd.color('magenta', 'cdr0')

CDR_Cavec = [CDR_Caline_fit.direction[0],CDR_Caline_fit.direction[1],CDR_Caline_fit.direction[2]]

CDR_Cafit = getDirPlane(CDR_CalinePoints[0], CDR_CalinePoints[3], allTRX_COM, TRA_COM, pointsBool = False)
####make_plane_points('CDR_Cafit', allTRX_COM, CDR_CalinePoints[0], CDR_CalinePoints[3])
BG_CDR_CafitPlane = getAngle(BGvec, CDR_Cafit, "PlaneVector")

###print("The angle between the Binding Groove vector model and the CDRs Ca plane model is:")
print("%f" % BG_CDR_CafitPlane)



###################################################################################################
############# fit line in 3D space for all CDRs Using All Atoms of CDR Residues ###################
### print("CDR ALL ATOMS")

CDR_atoms_coords = []
cmd.iterate_state(1, 'allCDR', 'CDR_atoms_coords.append([x,y,z])')
CDR_atomspoints = Points(CDR_atoms_coords)
CDR_atomsline_fit = Line.best_fit(CDR_atomspoints)

# plot CDRvector
cdr = [-20,-10,10,20]
CDR_atomslinePoints = [0,0,0,0]
holder = []
index = 0
for i in cdr:
    holder = CDR_atomsline_fit.point + i * CDR_atomsline_fit.direction
    holder = [holder[0],holder[1],holder[2]]
    CDR_atomslinePoints[index] = holder
    index = index + 1
#index = 0
#for i in CDR_atomslinePoints:
#    if cdr[index] < 0:
#        name = 'cdr' + 'neg' + str(abs(cdr[index]))
#    else:
#        name = 'cdr' + str(cdr[index])
#    cmd.pseudoatom(name, pos = i)
#    if t[index] < 0:
#        cmd.show_as('spheres', 'cdr' + 'neg' + str(abs(cdr[index])))
#        cmd.color('magenta', 'cdr' + 'neg' + str(abs(cdr[index])))
#    else:
#        cmd.show_as('spheres', 'cdr' + str(abs(cdr[index])))
#        cmd.color('magenta', 'cdr' + str(abs(cdr[index])))
#    index = index + 1
#cmd.pseudoatom('cdr0', pos = [CDR_atomsline_fit.point[0],CDR_atomsline_fit.point[1],CDR_atomsline_fit.point[2]])
#cmd.show_as('spheres', 'cdr0')
#cmd.color('magenta', 'cdr0')

#CDR_atomsvec = [CDR_atomsline_fit.direction[0],CDR_atomsline_fit.direction[1],CDR_atomsline_fit.direction[2]]

CDR_atomsfit = getDirPlane(CDR_atomslinePoints[0], CDR_atomslinePoints[3], allTRX_COM, TRA_COM, pointsBool = False)
####make_plane_points('CDR_atomsfit', allTRX_COM, CDR_atomslinePoints[0], CDR_atomslinePoints[3])
BG_CDR_atomsfitPlane = getAngle(BGvec, CDR_atomsfit, "PlaneVector")

###print("The angle between the Binding Groove vector model and the TCR plane (CDRs all atoms) model is:")
print("%f" % BG_CDR_atomsfitPlane)
###print("The angle between the Binding Groove vector model and the TCR germline plane is:")
print("%f" % BG_TCRgerm_fitPlane)


########################################################################################################################
########################## fit line in 3D space for all contact residues using Ca ######################################
if peptide != "Error":
    ContactCa_35_coords = []
    cmd.iterate_state(1, 'AlphaContactCa_3.5', 'ContactCa_35_coords.append([x,y,z])')
    cmd.iterate_state(1, 'BetaContactCa_3.5', 'ContactCa_35_coords.append([x,y,z])')
    ###print('\n%d points were used in the fitting of the ContactCa_35 Vector'% len(ContactCa_35_coords))
    print(len(ContactCa_35_coords))

    ContactCa_35_points = Points(ContactCa_35_coords)
    ContactCa_35_line_fit = Line.best_fit(ContactCa_35_points)

    # plot contact vector
    c = [-20,-10,10,20]
    ContactCa_35_linePoints = [0,0,0,0]
    holder = []
    index = 0
    for i in c:
        holder = ContactCa_35_line_fit.point + i * ContactCa_35_line_fit.direction
        holder = [holder[0],holder[1],holder[2]]
        ContactCa_35_linePoints[index] = holder
        index = index + 1
    #index = 0
    #for i in ContactCa_35_linePoints:
    #    if c[index] < 0:
    #        name = 'c' + 'neg' + str(abs(c[index]))
    #    else:
    #        name = 'c' + str(c[index])
    #    cmd.pseudoatom(name, pos = i)
    #    if t[index] < 0:
    #        cmd.show_as('spheres', 'c' + 'neg' + str(abs(c[index])))
    #        cmd.color('slate', 'c' + 'neg' + str(abs(c[index])))
    #    else:
    #        cmd.show_as('spheres', 'c' + str(abs(c[index])))
    #        cmd.color('slate', 'c' + str(abs(c[index])))
    #    index = index + 1
    #cmd.pseudoatom('c0', pos = [ContactCa_35_line_fit.point[0],ContactCa_35_line_fit.point[1],ContactCa_35_line_fit.point[2]])
    #cmd.show_as('spheres', 'c0')
    #cmd.color('slate', 'c0')

    ContactCa_35 = getDirPlane(ContactCa_35_linePoints[0], ContactCa_35_linePoints[3], allTRX_COM, TRA_COM, pointsBool= False)
    ####make_plane_points('ContactCa_3.5', allTRX_COM, ContactCalinePoints[0], ContactCalinePoints[3])
    BG_ContactCa_35_Plane = getAngle(BGvec, ContactCa_35, "PlaneVector")

    ###print("The angle between the Binding Groove vector model and the Contact Residue within 3.5 Angstroms Ca plane model is:")
    print("%f" % BG_ContactCa_35_Plane)
else:
    ###print("The angle between the Binding Groove vector model and the Contact Residue within 3.5 Angstroms Ca plane model is:")
    print("N/A")
    print("N/A")


if peptide != "Error":
    ContactCa_45_coords = []
    cmd.iterate_state(1, 'AlphaContactCa_4.5', 'ContactCa_45_coords.append([x,y,z])')
    cmd.iterate_state(1, 'BetaContactCa_4.5', 'ContactCa_45_coords.append([x,y,z])')
    ###print('\n%d points were used in the fitting of the ContactCa_45 Vector'% len(ContactCa_45_coords))
    print(len(ContactCa_45_coords))

    ContactCa_45_points = Points(ContactCa_45_coords)
    ContactCa_45_line_fit = Line.best_fit(ContactCa_45_points)

    # plot contact vector
    c = [-20,-10,10,20]
    ContactCa_45_linePoints = [0,0,0,0]
    holder = []
    index = 0
    for i in c:
        holder = ContactCa_45_line_fit.point + i * ContactCa_45_line_fit.direction
        holder = [holder[0],holder[1],holder[2]]
        ContactCa_45_linePoints[index] = holder
        index = index + 1
    #index = 0
    #for i in ContactCa_45_linePoints:
    #    if c[index] < 0:
    #        name = 'c' + 'neg' + str(abs(c[index]))
    #    else:
    #        name = 'c' + str(c[index])
    #    cmd.pseudoatom(name, pos = i)
    #    if t[index] < 0:
    #        cmd.show_as('spheres', 'c' + 'neg' + str(abs(c[index])))
    #        cmd.color('slate', 'c' + 'neg' + str(abs(c[index])))
    #    else:
    #        cmd.show_as('spheres', 'c' + str(abs(c[index])))
    #        cmd.color('slate', 'c' + str(abs(c[index])))
    #    index = index + 1
    #cmd.pseudoatom('c0', pos = [ContactCa_45_line_fit.point[0],ContactCa_45_line_fit.point[1],ContactCa_45_line_fit.point[2]])
    #cmd.show_as('spheres', 'c0')
    #cmd.color('slate', 'c0')

    ContactCa_45 = getDirPlane(ContactCa_45_linePoints[0], ContactCa_45_linePoints[3], allTRX_COM, TRA_COM, pointsBool= False)
    ####make_plane_points('ContactCa_4.5', allTRX_COM, ContactCalinePoints[0], ContactCalinePoints[3])
    BG_ContactCa_45_Plane = getAngle(BGvec, ContactCa_45, "PlaneVector")

    ###print("The angle between the Binding Groove vector model and the Contact Residue within 4.5 Angstroms Ca plane model is:")
    print("%f" % BG_ContactCa_45_Plane)
else:
    ###print("The angle between the Binding Groove vector model and the Contact Residue within 4.5 Angstroms Ca plane model is:")
    print("N/A")
    print("N/A")

########################################################################################################################
################### fit line in 3D space for all contact residues using All Atoms of Contact Residues ##################
### print('CONTACT ALL ATOMS')

if peptide != "Error":
    Contact_35_coords = []
    cmd.iterate_state(1, 'AlphaContact_3.5', 'Contact_35_coords.append([x,y,z])')
    cmd.iterate_state(1, 'BetaContact_3.5', 'Contact_35_coords.append([x,y,z])')
    ###print('\n%d points were used in the fitting of the Contact_35 Vector'% len(Contact_35_coords))
    print(len(Contact_35_coords))
    
    Contact_35_points = Points(Contact_35_coords)
    Contact_35_line_fit = Line.best_fit(Contact_35_points)

    # plot contact vector
    c = [-20,-10,10,20]
    Contact_35_linePoints = [0,0,0,0]
    holder = []
    index = 0
    for i in c:
        holder = Contact_35_line_fit.point + i * Contact_35_line_fit.direction
        holder = [holder[0],holder[1],holder[2]]
        Contact_35_linePoints[index] = holder
        index = index + 1
    #index = 0
    #for i in Contact_35_linePoints:
    #    if c[index] < 0:
    #        name = 'c' + 'neg' + str(abs(c[index]))
    #    else:
    #        name = 'c' + str(c[index])
    #    cmd.pseudoatom(name, pos = i)
    #    if t[index] < 0:
    #        cmd.show_as('spheres', 'c' + 'neg' + str(abs(c[index])))
    #        cmd.color('slate', 'c' + 'neg' + str(abs(c[index])))
    #    else:
    #        cmd.show_as('spheres', 'c' + str(abs(c[index])))
    #        cmd.color('slate', 'c' + str(abs(c[index])))
    #    index = index + 1
    #cmd.pseudoatom('c0', pos = [Contact_35_line_fit.point[0],Contact_35_line_fit.point[1],Contact_35_line_fit.point[2]])
    #cmd.show_as('spheres', 'c0')
    #cmd.color('slate', 'c0')

    Contact_35 = getDirPlane(Contact_35_linePoints[0], Contact_35_linePoints[3], allTRX_COM, TRA_COM, pointsBool=False)
    ####make_plane_points('Contact_3.5', allTRX_COM, ContactlinePoints[0], ContactlinePoints[3])
    BG_Contact_35_Plane = getAngle(BGvec, Contact_35, "PlaneVector")

    ###print("The angle between the Binding Groove vector model and the Contact Residue within 3.5 Ansgtroms all atoms plane model is:")
    print("%f" % BG_Contact_35_Plane)
else:
    ###print("The angle between the Binding Groove vector model and the Contact Residue within 3.5 Ansgtroms plane model is:")
    print("N/A")
    print("N/A")


if peptide != "Error":
    Contact_45_coords = []
    cmd.iterate_state(1, 'AlphaContact_4.5', 'Contact_45_coords.append([x,y,z])')
    cmd.iterate_state(1, 'BetaContact_4.5', 'Contact_45_coords.append([x,y,z])')
    ###print('\n%d points were used in the fitting of the Contact_45 Vector'% len(Contact_45_coords))
    print(len(Contact_45_coords))

    Contact_45_points = Points(Contact_45_coords)
    Contact_45_line_fit = Line.best_fit(Contact_45_points)

    # plot contact vector
    c = [-20,-10,10,20]
    Contact_45_linePoints = [0,0,0,0]
    holder = []
    index = 0
    for i in c:
        holder = Contact_45_line_fit.point + i * Contact_45_line_fit.direction
        holder = [holder[0],holder[1],holder[2]]
        Contact_45_linePoints[index] = holder
        index = index + 1
    #index = 0
    #for i in Contact_45_linePoints:
    #    if c[index] < 0:
    #        name = 'c' + 'neg' + str(abs(c[index]))
    #    else:
    #        name = 'c' + str(c[index])
    #    cmd.pseudoatom(name, pos = i)
    #    if t[index] < 0:
    #        cmd.show_as('spheres', 'c' + 'neg' + str(abs(c[index])))
    #        cmd.color('slate', 'c' + 'neg' + str(abs(c[index])))
    #    else:
    #        cmd.show_as('spheres', 'c' + str(abs(c[index])))
    #        cmd.color('slate', 'c' + str(abs(c[index])))
    #    index = index + 1
    #cmd.pseudoatom('c0', pos = [Contact_45_line_fit.point[0],Contact_45_line_fit.point[1],Contact_45_line_fit.point[2]])
    #cmd.show_as('spheres', 'c0')
    #cmd.color('slate', 'c0')

    Contact_45 = getDirPlane(Contact_45_linePoints[0], Contact_45_linePoints[3], allTRX_COM, TRA_COM, pointsBool=False)
    ####make_plane_points('Contact_4.5', allTRX_COM, ContactlinePoints[0], ContactlinePoints[3])
    BG_Contact_45_Plane = getAngle(BGvec, Contact_45, "PlaneVector")

    ###print("The angle between the Binding Groove vector model and the Contact Residue within 4.5 Ansgtroms all atoms plane model is:")
    print("%f" % BG_Contact_45_Plane)
else:
    ###print("The angle between the Binding Groove vector model and the Contact Residue within 4.5 Ansgtroms plane model is:")
    print("N/A")
    print("N/A")

#################################################################################################################
# torsion between Va and Vb
test4 = getAngle(CDRa_fit, CDRb_fit, "PlanePlane")
###print("The angle between the TRA and TRB planes is:")
print(test4)

######################################################################################################################
####
#### Distance Calculations for the BG center and its intersection with a given Plane
####
######################################################################################################################

Aplane = getPlane(CDRa_fit, TRA_COM)
Bplane = getPlane(CDRb_fit, TRB_COM)
CDR_Caplane = getPlane(CDR_Cafit, CDR_CalinePoints[0])
CDR_atomsplane = getPlane(CDR_atomsfit, CDR_atomslinePoints[0])
if peptide != "Error":
    #### 35 Angstroms ####
    CONTCaplane35 = getPlane(ContactCa_35, ContactCa_35_linePoints[0])
    CONT_atomsplane35 = getPlane(Contact_35, Contact_35_linePoints[0])
    BG_CONTCainters35 = solveIntersect(CONTCaplane35, BGvec, BGpoint)
    BG_CONTatomsinters35 = solveIntersect(CONT_atomsplane35, BGvec, BGpoint)
    BG_CONTCadist35 = getDist(BG_CONTCainters35, BGpoint)
    BG_CONT_atomsdist35 = getDist(BG_CONTatomsinters35, BGpoint)
    #### 45 Angstroms ####
    CONTCaplane45 = getPlane(ContactCa_45, ContactCa_45_linePoints[0])
    CONT_atomsplane45 = getPlane(Contact_45, Contact_45_linePoints[0])
    BG_CONTCainters45 = solveIntersect(CONTCaplane45, BGvec, BGpoint)
    BG_CONTatomsinters45 = solveIntersect(CONT_atomsplane45, BGvec, BGpoint)
    BG_CONTCadist45 = getDist(BG_CONTCainters45, BGpoint)
    BG_CONT_atomsdist45 = getDist(BG_CONTatomsinters45, BGpoint)

BG_Ainters = solveIntersect(Aplane, BGvec, BGpoint)
#cmd.pseudoatom('A_Inters', pos=BG_Ainters)
#cmd.show_as('spheres', 'A_Inters')
#cmd.color('limon', 'A_Inters')
BG_Binters = solveIntersect(Bplane, BGvec, BGpoint)
#cmd.pseudoatom('B_Inters', pos=BG_Binters)
#cmd.show_as('spheres', 'B_Inters')
#cmd.color('marine', 'B_Inters')
BG_CDR_Cainters = solveIntersect(CDR_Caplane, BGvec, BGpoint)
#cmd.pseudoatom('CDR_Ca_Inters', pos=BG_CDR_Cainters)
#cmd.show_as('spheres', 'CDR_Ca_Inters')
#cmd.color('olive', 'CDR_Ca_Inters')
BG_CDR_atomsinters = solveIntersect(CDR_atomsplane, BGvec, BGpoint)
#cmd.pseudoatom('CDR_atoms_Inters', pos=BG_CDR_atomsinters)
#cmd.show_as('spheres', 'CDR_atoms_Inters')
#cmd.color('olive', 'CDR_atoms_Inters')


BG_Adist = getDist(BG_Ainters, BGpoint)
BG_Bdist = getDist(BG_Binters, BGpoint)
BG_CDR_Cadist = getDist(BG_CDR_Cainters, BGpoint)
BG_CDR_atomsdist = getDist(BG_CDR_atomsinters, BGpoint)


###print("\n")
###print("The distance between the center of the Binding Groove and its intersection with the CDRa TRA plane is:")
print(BG_Adist)
###print("The distance between the center of the Binding Groove and its intersection with the CDRb TRB plane is:")
print(BG_Bdist)
###print("The distance between the center of the Binding Groove and its intersection with the CDRs Ca plane is:")
print(BG_CDR_Cadist)
###print("The distance between the center of the Binding Groove and its intersection with the TCR plane (CDRs all atoms) is:")
print(BG_CDR_atomsdist)
if peptide != "Error":
    ###print("The distance between the center of the Binding Groove and its intersection with the Contact Ca plane (within 3.5 Angstroms) is:")
    print(BG_CONTCadist35)
    ###print("The distance between the center of the Binding Groove and its intersection with the Contact plane (within 3.5 Angstroms) is:")
    print(BG_CONT_atomsdist35)
    ###print("The distance between the center of the Binding Groove and its intersection with the Contact Ca plane (within 4.5 Angstroms) is:")
    print(BG_CONTCadist45)
    ###print("The distance between the center of the Binding Groove and its intersection with the Contact plane (within 4.5 Angstroms) is:")
    print(BG_CONT_atomsdist45)
else:
    print("N/A")
    print("N/A")
    print("N/A")
    print("N/A")

cmd.show('cartoon', PDB + "_" + MHCa)
#cmd.show('cartoon', PDB + "_" + MHCb)
cmd.show('cartoon', PDB + "_" + TRAchain)
cmd.show('cartoon', PDB + "_" + TRBchain)
if peptide != "None" and peptide != "Error":
    cmd.show('cartoon', PDB + "_" + peptide)
