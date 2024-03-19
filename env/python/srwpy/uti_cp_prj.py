# -*- coding: utf-8 -*-
#############################################################################
# SRWLib Utility for copying SRW files
# v 0.01
# Authors: O.C.
#############################################################################

import argparse
import os
import shutil
#import copy

try: #HG24052023
    from . import uti_io
except:
    import uti_io
import uti_io

#****************************************************************************
#Usage example.
#To create "srw_python" directory in the directory "C:\SoftwareDevelopments\SRW_Dev\local\dev_contrib\from_Rebecca",
#run this from a work "srw_python" directory:
#python uti_cp_prj.py -f proj_files.txt -d C:\SoftwareDevelopments\SRW_Dev\local\dev_contrib\from_Rebecca --nsbd --nex --ntx
#****************************************************************************

if __name__=='__main__':

    p = argparse.ArgumentParser()
    
    p.add_argument('-f', '--flist', dest='file_list', type=str, default='proj_files.txt', help='file containing list of files to be copied')
    p.add_argument('-d', '--dest', dest='dest_dir', type=str, default='', help='destination directory project files to be copied to')
    p.add_argument('--nsbd', dest='no_sub_dirs', default=False, action='store_true', help='copy or not sub-directories')
    p.add_argument('--nex', dest='no_examples', default=False, action='store_true', help='copy or not examples')
    p.add_argument('--ntx', dest='no_txt', default=False, action='store_true', help='copy or not text files')

    v = p.parse_args()

    if not os.path.isfile(v.file_list):
        p.error('File name containing list of files is not specified or is incorrect. Use -f option to specify the file name with path.')

    workDir = os.getcwd()
    newDirPath = os.path.join(v.dest_dir, os.path.basename(workDir))

    if(not os.path.exists(newDirPath)): os.mkdir(newDirPath)

    fileList = uti_io.read_ascii_data_cols(v.file_list, '\n', _float=False)[0]

    def checkIfInList(_list, _item):
        for i in range(len(_list)):
            if(_list[i] == _item): return True
        return False
    
    subDirName = None
    listSubDirsCreated = []
    for i in range(len(fileList)):
        curFileName = fileList[i]
        #print(curFileName)

        if(v.no_examples):
            indExam = curFileName.lower().find('example')
            if(indExam >= 0): continue

        indSlash = curFileName.find('/')
        if(indSlash >= 0):
            if(not v.no_sub_dirs):
                subDirName = curFileName[0:indSlash]
                #curFileName = curFileName[(indSlash+1):]
                curDestDir = os.path.join(newDirPath, subDirName)
        
                if(not checkIfInList(listSubDirsCreated, subDirName)):
                    newSubDirPath = os.path.join(newDirPath, subDirName)
                    if(not os.path.exists(newSubDirPath)): os.mkdir(newSubDirPath)
                    listSubDirsCreated.append(subDirName)

                #print(curFileName)
                if(os.path.exists(curFileName)): shutil.copy(curFileName, curDestDir)
            continue

        if(v.no_txt):
            lenFileName = len(curFileName)
            if(lenFileName >= 4):
                if(curFileName.lower()[lenFileName-4:lenFileName] == '.txt'): continue

        #print(curFileName)
        if(os.path.exists(curFileName)): shutil.copy(curFileName, newDirPath)
