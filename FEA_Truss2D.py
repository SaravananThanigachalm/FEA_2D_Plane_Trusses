import numpy as np
import math

np.set_printoptions(3, suppress= True) #Round-off decimal to 3 significant digits after decimal
print("\nWelcome to 2D Truss FEA Solver V1")
print("By Saravanan T\n")

#<--------------- Elements and Nodes --------------->#
Elements = int(input("Please enter the number of Elements: ")) # Specify the number of Elements in the Truss
Nodes = int(input("Please enter the number of Nodes: ")) # Specify the number of Nodes in the Truss
Xco = []    # X Co-ordinates of Nodes
Yco = []    # Y Co-ordinates of Nodes

for i in range(Nodes):
    X = float(input("Enter the X co-ordinate of the node "+str(i+1)+" in mm : "))
    Y = float(input("Enter the Y co-ordinate of the node "+str(i+1)+" in mm : "))
    Xco.append(X)
    Yco.append(Y)

A = []  # Area of Individual Elements
E = []  # Young's modulus of Individual Elements in N/mm^2
for i in range(Elements):
    a = float(input("Enter the Cross-sectional Area of the element "+str(i+1)+" in mm^2 : "))
    e = float(input("Enter the Young's modulus of the element "+str(i+1)+" in N/mm^2 : "))
    A.append(a)
    E.append(e)

StNode = [] # Starting Nodes of Individual Elements
EdNode = [] # Ending Nodes of Individual Elements
L = [] # Length of Individual Elements
Cons = [] # Constants of Individual Elements
CosNode = [] # Cosine of Individual Elements
SinNode = [] # Sine of Individual Elements
for i in range(Elements):
    s = int(input("Enter the Start node of the Element "+str(i+1)+" : "))
    en = int(input("Enter the End node of the Element "+str(i+1)+" : "))
    x1 = float(Xco[s - 1])
    y1 = float(Yco[s - 1])
    x2 = float(Xco[en - 1])
    y2 = float(Yco[en - 1])
    a = float(A[i])
    e = float(E[i])
    l = math.sqrt((x2-x1)**2+(y2-y1)**2)
    con = (a * e) / l
    cos = (x2-x1) / l
    sin = (y2-y1) / l

    StNode.append(s)
    EdNode.append(en)
    L.append(l)
    Cons.append(con)
    CosNode.append(cos)
    SinNode.append(sin)

#<--------------- Global Element Stiffness Matrix --------------->#
GoElStiff = []
for i in range(Elements):
    c2 = float(CosNode[i]**2)
    s2 = float(SinNode[i]**2)
    cs = float(CosNode[i] * SinNode[i])
    emat = Cons[i] * np.array([[c2, cs, -c2, -cs],
                                [cs, s2, -cs, -s2],
                                [-c2, -cs, c2, cs],
                                [-cs, -s2, cs, s2]]) # This appends the elemental constants into Global Stiffness Matrix form

    GoElStiff.append(emat)

#<--------------- Assembling the Elemental Global Stiffness matrix into Assemblage Stiffness Matrix --------------->#
gelstill = [] # Assemblage stiffness matrix mapping
for i in range(Elements):
    m = StNode[i] * 2
    n = EdNode[i] * 2
    add = [m-1, m, n-1, n]             # Address of columns and rows of gelstill for elemet(i)
                                        # if startnode is 1 and end node is 2 then add=[1,2,3,4]
                                        # if startnode is 1 and end node is 3 then add=[1,2,5,6]
    gmat = np.zeros((Nodes*2, Nodes*2)) 
    elmat = GoElStiff[i]                
    for j in range(4):
        for k in range(4):
            a = add[j]-1
            b = add[k]-1
            gmat[a,b] = elmat[j,k]  # Updating the values in GST matrix with EST matrix of element(i)
    gelstill.append(gmat)           # Storing the resultant matrix in gelstill list

GASM = np.zeros((Nodes*2, Nodes*2))
for mat in gelstill:
    GASM = GASM + mat       # Adding all the matrix in the gelstill list

""" This Generates the Assemblage Stiffness Matrix of the respective Truss Structure """
# print("\n Assemblage Stiffness Matrix :-\n")
# print(GASM)

#<--------------- Displacement and Force Matrices --------------->#
Disp = []
Force = []
for i in range(Nodes):
    a = "u"+str(i+1)
    Disp.append(a)
    b = "v"+str(i+1)
    Disp.append(b)
    c = "fx"+str(i+1)
    Force.append(c)
    Disp.append(b)
    d = "fy"+str(i+1)
    Force.append(d)

#<--------------- Support Constraints --------------->#
print("\nSupport Specification\n")
ResMat = np.ones((Nodes*2,1))
TSupNo = int(input("Specify the number of Nodes have supports: "))
tcondition = ["P for Pinned ",
            "H for horizontally constrained but movable in the Vertical direction",
            "V for Vertically constrained but movable in the Horizontal direction"]

for i in range(TSupNo):
    SupNo = int(input("Enter the Node number that is supported: "))
    for j in tcondition:
        print(j)
    condition = input("Specify the condition at the node: ")
    if condition in ['P','p']:
        ResMat[SupNo*2-1, 0] = 0
        ResMat[SupNo*2-2, 0] = 0
    elif condition in ['H','h']:
        ResMat[SupNo*2-2, 0] = 0
    elif condition in ['V','v']:
        ResMat[SupNo*2-1, 0] = 0
    else:
        print("Please enter a valid option!!!")
        
#<--------------- Load Constraints --------------->#
print("\nLoading Specification\n")
FMat = np.zeros((Nodes*2,1))
tload = int(input("Enter the number of Loaded Nodes: "))
for i in range(tload):
    load = int(input("Enter the Node number that is Loaded: "))
    fx = float(input("Enter the Horizontal Load at that Node: "))
    fy = float(input("Enter the Vertical Load at that Node: "))
    FMat[load*2 - 2, 0] = fx
    FMat[load*2 - 1, 0] = fy

#<--------------- Matrix reduction --------------->#
rcdlist = []
for i in range(Nodes*2):
    if ResMat[i,0] == 0:
        rcdlist.append(i)

rrgsm = np.delete(GASM, rcdlist, 0) #row reduction
crgsm = np.delete(rrgsm, rcdlist, 1) #column reduction
rgsm = crgsm #reduced Assembleage stiffness matrix
rforcemat = np.delete(FMat, rcdlist, 0) #reduced force mat
rdispmat = np.delete(ResMat, rcdlist, 0) #reduced disp mat

#<--------------- Solving --------------->#
dispresult = np.linalg.solve(rgsm, rforcemat) # Using Matrix Inversion Method
rin = 0
for i in range(Nodes*2):
    if ResMat[i,0] == 1:
        ResMat[i,0] = dispresult[rin,0]
        rin = rin+1

forceresult = np.matmul(GASM, ResMat)

print("\n\nAssemblage Stiffness Matrix of the Truss :-\n")
print(GASM)
print("\n\nDisplacement matrix of nodes :-\n")
print(ResMat)
print("\n\nForce matrix of nodes :-\n")
print(forceresult)

#<--------------- Updating the Nodal Co-ordinates w.r.t Deformation --------------->#
newxco = []     # Updating the respective Nodes w.r.t to the Displacements (i.e) Deformations
newyco = []     # Updating the respective Nodes w.r.t to the Displacements (i.e) Deformations
count = 0
for i in range(Nodes):
    k = Xco[i]+ResMat[count,0]
    newxco.append(k)
    count = count+1
    l = Yco[i]+ResMat[count,0]
    newyco.append(l)
    count = count+1

#<--------------- Updating the Length of Elements w.r.t Nodal Co-ordinates --------------->#
NewL = [] # Generating New lengths of Elements for new Co-ordinates
for i in range(Elements):
    a, b = StNode[i], EdNode[i]
    x1 = float(newxco[a-1])
    y1 = float(newyco[a-1])
    x2 = float(newxco[b-1])
    y2 = float(newyco[b-1])
    l = math.sqrt((x2-x1)**2+(y2-y1)**2)
    NewL.append(l)

#<--------------- Solving for Elemental Strains --------------->#
Elstrain = np.zeros((Elements,1))
for i in range(Elements):
    Elstrain[i,0] = (NewL[i]-L[i])/(L[i])

print("\n***Positive is Tensile\nNegetive is Compressive***")
print("\n\nStrain in the elements :-\n")
print(Elstrain)

#<--------------- Solving for Elemental Stresses --------------->#
Elstress = np.zeros((Elements,1))
for i in range(Elements):
    Elstress[i,0] = E[i] * Elstrain[i,0] 

print("\n\nStress in the elements :-\n")
print(Elstress)

#<--------------- Solving for Elemental Forces --------------->#
Eforce = np.zeros((Elements,1))
for i in range(Elements):
    Eforce[i,0] = A[i] * Elstress[i,0]

print("\n\nForce in the element :-\n")
print(Eforce)
