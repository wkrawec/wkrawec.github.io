##
# Copyright (c) 2025

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##

# The following is the code we used to test our quantum PKE implementation for the paper "Quantum Public Key Encryption for NISQ Devices".  It does not use a quantum secure PRF so is not an actual secure implementation.  Nor does it use any actual privacy amplification.  Thus, this program does not actually encrypt or decrypt anything.  Instead, it tests the circuit implementations for the measurement operations and evaluates the noise.

# "runNoiseTest" and "runExperiment" are the main functions.

from azure.quantum import Workspace
workspace = Workspace (
   resource_id = "",
   location = "eastus"
)

from azure.quantum import Workspace
from azure.quantum.qiskit import AzureQuantumProvider
provider = AzureQuantumProvider(workspace)

#quantinuum_api_val_backend = provider.get_backend("quantinuum.sim.h1-1e")
quantinuum_api_val_backend = provider.get_backend("quantinuum.qpu.h1-1")
from qiskit import QuantumCircuit
from qiskit.visualization import plot_histogram

import random

def getNthBit(i, n):
    for k in range(n):
        i = i >> 1
    return i & 1

def flipNthBit(i, n):
    j = 1
    for k in range(n):
        j = j << 1
    return i^j

def chooseRandomInt(maxval):
    return random.randint(0,maxval-1)
    
def chooseDistinctRandomInts(maxval):
    i = random.randint(0,maxval-1)
    j = random.randint(0,maxval-1)
    while j == i:
        j = random.randint(0,maxval-1)

    return (i,j)
    

# Some helper functions
def applyX(base, wire):
    global circuit
    circuit.x(base+wire)
    #print("Apply X to ", base+wire)

def applyH(base, wire):
    global circuit
    circuit.h(base+wire)

def applyT(base, control1, control2, target):
    global circuit
    circuit.ccx(control1, control2, target)

def applyCNOT(base, control, target):
    global circuit
    circuit.cx(base+control, base+target)
    #print("Apply CNOT with c = ", base+control, "; t = ", base+target)
    
def applySwap(base, i, j):
    global circuit
    if i == j:
        return
    circuit.cx(base+i, base+j)
    circuit.cx(base+j, base+i)
    circuit.cx(base+i, base+j)
    #circuit.swap(base+i, base+j)
    #print("Apply SWAP between ", base+i, " and ", base+j)

def interpretResultNew(blockNum, numWires, dataArray, delta):
    # First convert result to an integer:
    x = 0
    p = 1
    for w in range(blockNum*numWires, (blockNum+1)*numWires):
        x = x + p*int(dataArray[w])
        p = p << 1

    if (x%2) == 0:
        i = permInv[delta-1][x]
        j = permInv[delta-1][x+1]
        w = 0
    else:
        i = permInv[delta-1][x-1]
        j = permInv[delta-1][x]
        w = 1

    return(i, j, w)

def interpretResult(blockNum, numWires, dataArray):
    allZero = True
    for w in range(blockNum*numWires+1, (blockNum+1)*numWires):
        if int(dataArray[w]):
            allZero = False
            break

    print("allzero = ", allZero)
    if allZero:
        x = int(dataArray[blockNum*numWires]) ## this should be F(i) + F(j)
        print("x = ", x)
        return x

    return -1

def ConvertStrToBits(message):
    res = ''.join(format(ord(i), '08b') for i in message)

    return res

def ConvertBitsToStr(bits):
    res = ''

    for i in range(0, len(bits), 8):
        temp = int(bits[i:i+8], 2)
        res = res + chr(temp)

    return res


# numWires is the total number of wires in a block (not counting the extra aux wire)
# Thus it is realy the security parameter; input to PRF
# block num is the block index (with one block being numWires+1)
def buildIJSwap(i, j, numWires, blockNum):
    for k in range(numWires):
        if getNthBit(i, k) == 1:
            applyX(blockNum*(numWires+1), k)
            j = flipNthBit(j, k)

    minK = numWires+1
    for k in range(numWires):
        if getNthBit(j, k) == 1:
            minK = k
            break
    if minK > numWires:
        print("Internal Error in buildIJSwap: i = j")
        return 0

    for l in range(k+1, numWires):
        if getNthBit(j, l) == 1:
            applyCNOT(blockNum*(numWires+1),minK, l)

    applySwap(blockNum*(numWires+1), minK, 0)

def PRF(key, input):
    x = getNthBit(input, 2)
    if key[1]:
        x = x ^ getNthBit(input, 1)
    if key[0]:
        x = x ^ 1

    return x
def buildUF(base, numWires):
    global circuit
    circuit.cx(base+2, base+(numWires))
    if(enckey[1]):
        circuit.cx(base+1, base+(numWires))
    if (enckey[0]):
        circuit.x(base+(numWires))


# New code for each delta; we only have up to n=3 qubits:
def buildDelta1(numWires, blockNum):
    base = blockNum*(numWires+1)

def buildDelta2(numWires, blockNum):
    base = blockNum*(numWires+1)
    applySwap(base, 1,0)

def buildDelta4(numWires, blockNum):
    base = blockNum*(numWires+1)
    applySwap(base, 2, 0)

def buildDelta3(numWires, blockNum):
    base = blockNum*(numWires+1)
    applySwap(base, 0, 1)
    applySwap(base, 0, 2)
    applyT(base, 1, 2, 0)
    applyCNOT(base, 0, 1)
    applyCNOT(base, 0, 2)
    applyT(base, 1, 2, 0)

def buildDelta7(numWires, blockNum):
    base = blockNum*(numWires+1)
    applySwap(base, 0, 1)
    applyCNOT(base, 0, 1)
    applySwap(base, 0, 2)
    applyCNOT(base, 0, 2)
    applyT(base, 1, 2, 0)
    applyT(base, 0, 1, 2)

def buildDelta5(numWires, blockNum):
    base = blockNum*(numWires+1)
    applySwap(base, 0, 1)
    applySwap(base, 0, 2)
    applyCNOT(base, 0, 1)

def buildDelta6(numWires, blockNum):
    base = blockNum*(numWires+1)
    applySwap(base, 0, 1)
    applySwap(base, 0, 2)
    applyCNOT(base, 0, 2)

def buildDelta(d, numWires, blockNum):
    match d:
        case 1: buildDelta1(numWires, blockNum)
        case 2: buildDelta2(numWires, blockNum)
        case 3: buildDelta3(numWires, blockNum)
        case 4: buildDelta4(numWires, blockNum)
        case 5: buildDelta5(numWires, blockNum)
        case 6: buildDelta6(numWires, blockNum)
        case 7: buildDelta7(numWires, blockNum)

permInv = []
def buildPermutations():
    global permInv
    permInv = []
    # delta = 1:
    perm1 = [0,1,2,3,4,5,6,7]
    perm2 = [0,2,1,3,4,6,5,7]
    perm3 = [0,3,1,6,2,5,4,7]
    perm4 = [0,4,2,6,1,5,3,7]
    perm5 = [0,5,1,4,2,7,3,6]
    perm6 = [0,6,1,7,2,4,3,5]
    perm7 = [0,7,1,2,3,4,5,6]

    permInv.append(perm1)
    permInv.append(perm2)
    permInv.append(perm3)
    permInv.append(perm4)
    permInv.append(perm5)
    permInv.append(perm6)
    permInv.append(perm7)

# numWires is the security parameter
# numBlocks is how many groups of key to create
# total qubits required = numBlocks*(numWires+1)
circuit = 0
totalQubits = 0
# ijChoices = []
deltaChoices = []

def setupCircuit(numWires, numBlocks):
    global totalQubits
    global circuit
    totalQubits = numBlocks*(numWires+1)
    circuit = QuantumCircuit(numBlocks*(numWires+1), numBlocks*numWires)

def createPubKey(numWires, numBlocks):
    global circuit
    # First setup Aux wires to have |1>:

    for i in range(numBlocks):
        applyX(i*(numWires+1), numWires)

    ## apply H everywhere:
    for i in range(numBlocks*(numWires+1)):
        circuit.h(i)

    circuit.barrier(range(totalQubits))
    ## apply UF to each block:
    for i in range(numBlocks):
        buildUF(i*(numWires+1), numWires)

    circuit.barrier(range(totalQubits))

def EncryptCircuit(numWires, numBlocks):
    global circuit
    global deltaChoices
    deltaChoices = []
    for b in range(numBlocks):
        #(i, j) = chooseDistinctRandomInts(2**numWires)
        delta = chooseRandomInt((2**numWires) - 1)+1
        #while (delta == 3) or (delta == 7):
        ##while (delta == 7):
        delta = chooseRandomInt((2**numWires) - 1) + 1
        #delta = 3
        buildDelta(delta, numWires, b)
        #buildIJSwap(i, j, numWires, b)
        deltaChoices.append( delta )  ## save for later - not all will be used due to measurement failure

    circuit.barrier(range(totalQubits))

    for b in range(numBlocks):
        circuit.h(b*(numWires+1))

    circuit.barrier(range(totalQubits))
    ## now measure:
    for b in range(numBlocks):
        for w in range(numWires):
            circuit.measure(b*(numWires+1) + w, b*numWires + w)

def Decrypt(ct_ij, ct_pad, key):
    pad = []
    for i in range(len(ct_ij)):
        pad.append(PRF(key, ct_ij[i][0]) ^ PRF(key, ct_ij[i][1]))

    f = hashFct(pad, len(ct_pad))
    pt = ''
    for i in range(len(ct_pad)):
        b = ct_pad[i] ^ f[i]
        if b == 0:
            pt += '0'
        else:
            pt += '1'

    return ConvertBitsToStr(pt)

## This will run the circuit and interpret measurement results

ctIJ = []
deltaChoices = []
ctPad = []
prePad = []
messageArray = [1,0,1]
prePadSize = 0

enckey = [0,0,0]
deckey = [0,0,0]

def SetupCipherText(message):
    global ctIJ
    global ctPad
    global prePad
    global messageArray
    global prePadSize
    global circuit
    global deltaChoices

    ctIJ = []
    ctPad = []
    prePad = []
    deltaChoices = []
    messageArray = ConvertStrToBits(message)
    print(messageArray)
    prePadSize = 1*len(messageArray)
    
def CreateNextCipherText(numWires, numBlocks):
    global ctIJ
    global ctPad
    global prePad
    global messageArray
    global prePadSize
    global circuit
    global deltaChoices

    buildPermutations()

    ## run the job:
    quantinuum_api_val_backend = provider.get_backend("quantinuum.sim.h1-1e")

    job = quantinuum_api_val_backend.run(circuit, shots=1)
    print("Job id:", job.id())

    result = job.result().results[0].data.counts

    plot_histogram(job.result().get_counts(circuit), title="Result")

    result_arr = max(result, key= lambda key: result[key])[::-1]

    print("MEASUREMENT RESULT = ", result_arr)

    ## now interpret the result block-by-block:
    for b in range(numBlocks):
        i,j,pad = interpretResultNew(b, numWires, result_arr, deltaChoices[b])
        if pad >= 0:
            ctIJ.append([i,j])
            ##m = int(messageArray[messageIndex])
            prePad.append(pad)

            print("i = ", i)
            print("j = ", j)
            print("pad = ", pad)

            if len(prePad) >= prePadSize:
                return True
            
    return False ## more to enc

def hashFct(x, desiredSize):
    t = []
    for i in range(desiredSize):
        t.append(x[i])

    return t
        
def FinalizeCiphertext():
    global ctIJ
    global ctPad
    global prePad
    global messageArray
    global prePadSize
    global circuit
    f = hashFct(prePad, len(messageArray))

    for i in range(len(messageArray)):
        ctPad.append(f[i] ^ int(messageArray[i]))

def runExperiment(securityParameter, numPubKeysPerTrial, message, key):
    global ctIJ
    global ctPad
    global prePad
    global messageArray
    global prePadSize
    global circuit
    #global ijChoices
    global deltaChoices
    SetupCipherText(message)

    print(messageArray)
    while len(prePad) < prePadSize:
        setupCircuit(securityParameter, numPubKeysPerTrial)
        #circuit.x(0)
        #for b in range(numPubKeysPerTrial):
        #    for w in range(securityParameter):
        #        circuit.measure(b*(securityParameter+1) + w, b*securityParameter + w)
        createPubKey(securityParameter, numPubKeysPerTrial)
        EncryptCircuit(securityParameter, numPubKeysPerTrial)

        print(deltaChoices)

        print(circuit.draw(fold=300))

        CreateNextCipherText(securityParameter, numPubKeysPerTrial)
        #break

    FinalizeCiphertext()

    ##print("prepad = ", prePad)
    print("ctIJ = ", ctIJ)
    print("ctPad = ", ctPad)

    return (ctIJ, ctPad)

def runNoiseTest_Circuit(numWires, numBlocks):
    global deltaNoise
    global deltaCounts

    global noiseTotalCount
    global ctIJ
    global ctPad
    global prePad
    global messageArray
    global prePadSize
    global circuit
    #global ijChoices
    global deltaChoices

    buildPermutations()

    ## run the job:
    quantinuum_api_val_backend = provider.get_backend("quantinuum.sim.h1-1e")

    job = quantinuum_api_val_backend.run(circuit, shots=1)
    print("Job id:", job.id())

    result = job.result().results[0].data.counts

    result_arr = max(result, key= lambda key: result[key])[::-1]

    ## now interpret the result block-by-block:
    for b in range(numBlocks):
        i,j,pad = interpretResultNew(b, numWires, result_arr, deltaChoices[b])
        if pad >= 0:
            ctIJ.append([i,j])
            ##m = int(messageArray[messageIndex])
            prePad.append(pad)
            epad = PRF(enckey,i) ^ PRF(enckey, j)

            #epad = chooseRandomInt(2)
            print("i = ", i)
            print("j = ", j)
            print("pad = ", pad)
            print("Expected pad = ", epad)

            if pad != epad:
                deltaNoise[deltaChoices[b]] += 1
            deltaCounts[deltaChoices[b]] += 1

def runNoiseTest(numRounds):
    global deltaNoise
    global deltaCounts

    global ctIJ
    global ctPad
    global prePad
    global messageArray
    global prePadSize
    global circuit
    #global ijChoices
    global deltaChoices
    

    noiseTotalCount = 0
    deltaNoise=[0,0,0,0,0,0,0,0,0]
    deltaCounts=[0,0,0,0,0,0,0,0,0]
    
    p = 3 # num parallel instances

    for rnd in range(numRounds):
        setupCircuit(3, p)
        createPubKey(3, p)
        EncryptCircuit(3, p)

        runNoiseTest_Circuit(3,p)

        print("Done with ",rnd+1, " rounds.")
        print("Noise = ", deltaNoise)
        print("Delta Counts = ", deltaCounts)

        totalError = 0
        totalPads = 0

        for i in range(len(deltaNoise)):
            totalError += deltaNoise[i]
            totalPads += deltaCounts[i]
            if deltaCounts[i] > 0:
                print("Error for Delta = ", i, ":\t", float(deltaNoise[i])/float(deltaCounts[i]))

        print("Total Error = ", float(totalError)/float(totalPads))

CT = runExperiment(3, 1, "hi", [1,0,0])
