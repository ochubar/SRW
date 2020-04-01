
from mpi4py import MPI
from array import *
import subprocess
import sys
import time

#-------------------------------------------------------------

comMPI = MPI.COMM_WORLD #MPI stuff
rank = comMPI.Get_rank()
nProc = comMPI.Get_size()

nmProtein = sys.argv[1] #./smiles_dock.sh script arguments
nmFile2Proc = sys.argv[2]
sCrd = sys.argv[3]
sGrid = sys.argv[4]

arInds = array('i', [0])
arRes = array('i', [0])

if(rank == 0): #Master
    if(nProc == 1): #There are no Workers, do all the work foe them without splitting the main file
        s2run = './smiles_dock.sh ' + nmProtein + ' ' + nmFile2Proc + ' \'' + sCrd + '\' \'' + sGrid + '\''
        print('rank/size', rank, '/', nProc, ' command to run: ', s2run)
        subprocess.run(['./smiles_dock.sh', nmProtein, nmFile2Proc, sCrd, sGrid]) #Do Processing
        exit()

    #split the main input file, save sub-files, send job orders to Workers, and receive reports from them
    f = open(nmFile2Proc, 'r')
    lines = f.readlines()
    f.close()
    nLines = len(lines)

    nWorkersBusy = 0
    for i in range(nLines): #Send messages to Workers to process sub-files
        nmSubFile = nmFile2Proc + '.' + repr(i)
        fs = open(nmSubFile, 'w')
        fs.write(lines[i])
        fs.close()

        arInds[0] = i
        if(i < nProc - 1): #First, just send sub-file index to a Worker
            comMPI.Send([arInds, MPI.INT], dest=i+1)
            nWorkersBusy += 1
        else: #If there are more sub-files than Workers, wait to hear from anyone that his previous job is done, and send a new one
            comMPI.Recv([arRes, MPI.INT], source=MPI.ANY_SOURCE) #Receive OK from any Worker that sub-file is processed
            rankWorker = arRes[0]
            comMPI.Send([arInds, MPI.INT], dest=rankWorker)

    arInds[0] = -1
    for i in range(nWorkersBusy): #Wait for Workers to finish and Send them messages to exit 
        comMPI.Recv([arRes, MPI.INT], source=MPI.ANY_SOURCE) #Receive OK from any Worker that sub-file is processed
        rankWorker = arRes[0]
        comMPI.Send([arInds, MPI.INT], dest=rankWorker)

    for i in range(nLines, nProc): #Send messages to Workers to exit 
        comMPI.Send([arInds, MPI.INT], dest=i)

else: #Workers: receive job orders from Master and report back after each job is done
    arRes[0] = rank
    maxNumSubFiles = 1000000 #Very Large Number
    for i in range(maxNumSubFiles):
        
        comMPI.Recv([arInds, MPI.INT], source=0) #Receive next job order from Master
        if(arInds[0] < 0): #Order to exit
            print('rank/size', rank, '/', nProc, ' exiting')
            break
        else:
            nmSubFile = nmFile2Proc + '.' + repr(arInds[0])

            s2run = './smiles_dock.sh ' + nmProtein + ' ' + nmSubFile + ' \'' + sCrd + '\' \'' + sGrid + '\''
            print('rank/size', rank, '/', nProc, ' command to run: ', s2run)
            subprocess.run(['./smiles_dock.sh', nmProtein, nmSubFile, sCrd, sGrid]) #Do Processing
            #Analyze results? 
            #time.sleep(2.) #Mimic execution (for debugging)

            comMPI.Send([arRes, MPI.INT], dest=0) #Report to Master that job was done
