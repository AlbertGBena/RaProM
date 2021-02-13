
import numpy as np
import calendar
import datetime
import time

import glob
import os
import sys
import shutil

##import matplotlib.dates as mdates

def date2unix(date):
    return calendar.timegm(date.timetuple())
def unix2date(unix):
    return datetime.datetime.utcfromtimestamp(unix)

def CorrectorFile(fid):
    NameFile=fid 

    FileCorre=NameFile[:-4]+'-corrected.raw' #create a new file

    FileCorre2=NameFile[:-4]+'-corrected-2.raw' #create a new file

    

    folderName='Moved'
    if not os.path.exists(folderName):
        os.mkdir(folderName)
    
    f1=open(FileCorre,'w+')

    f=open(NameFile,'r')
    file_length = len(f.read().split('\n'))
    totallines=(file_length-7)/67
    
    m=0
    timeList=list()
    f.close()
########CHECKING FOR ERRORS IN THE RAW FILE
    f=open(NameFile,'r')
    cond=1
    lineCount=0
    LineCount=0
    countF=0
    
    TotalLinesFile=0
    while cond:
        
        line=f.readline()
        TotalLinesFile+=1#lines counter
        if line=='':
            
            break

    
    
        line2=line.strip()
        columns=line2.split()

        valor=columns[0]

                
        
        if valor[0]=='M' or valor[0]=='H' or valor[0]=='T' or valor[0]=='F':
            LineCount+=1
            if valor[0]=='M':
                if countF!=64 or countF==0:
                    
                    if countF!=0:
                        lineCount=lineCount+LineCount#sumo el nombre de lines
                    
                    
                else:
                    f1.write(Line)
                
                Line=line
                LineCount=0#reset numero lines
                countF=0#rest el contador de les F
                

            else:
                Line=Line+line
                

            
            if valor[0]=='F':
                countF+=1
            if valor[0]=='T':
                LineTF=line
        else:
            
            
            if valor[0]=='C' and len(columns)>4:

                if columns[4]=='"TF':
                    Line=Line+LineTF
##                    print(line[25:-3])
##                    print(TotalLinesFile)
                else:
                    lineCount+=1
            else:
                lineCount+=1
            
            
    
    f1.close()
    f.close()

    if lineCount==0:
        os.remove(FileCorre)
        OutName=NameFile
        f1=open(FileCorre,'w+')
        f=open(NameFile,'r')
        
        
    else:
        shutil.copy(os.path.join('folder', NameFile), folderName)
        os.remove(NameFile)
        OutName=FileCorre

        f1=open(FileCorre2,'w+')
        f=open(FileCorre,'r')

########CHECKING FOR JUMPS TIME FOR THE SYNCHRONIZATION

    

    while cond:
        
        line=f.readline()
        if line=='':
            break

    
        line2=line.strip()
        columns=line2.split()
    
        if len(columns)==1:
            longStr=2
        else:
            Date=columns[1]
            longStr=len(str(Date))

       
        
            
        
        if len(timeList)==0:
            dat = datetime.datetime(year = 2000+int(Date[0:2]), month = int(Date[2:4]), day = int(Date[4:6]), hour = int(Date[6:8]), minute = int(Date[8:10]), second = int(Date[10:12]))
            dat=int(date2unix(dat))
            for j in range(67):
                
                f1.write(line)
                
                if j<66:
                    line=f.readline()
            timeList.append(dat)
        else:

            if longStr==12:
                dat = datetime.datetime(year = 2000+int(Date[0:2]), month = int(Date[2:4]), day = int(Date[4:6]), hour = int(Date[6:8]), minute = int(Date[8:10]), second = int(Date[10:12]))
                dat=int(date2unix(dat))
                if timeList[-1]<dat:
                    
                    for j in range(67):
                        f1.write(line)
                        if j<66:
                            line=f.readline()
                    timeList.append(dat)

            
            else:
                m=m+1

          
    f1.close()
    f.close()
    if m==0 and lineCount==0:
        os.remove(FileCorre)

        OutName=NameFile
    if m==0 and lineCount!=0:

        os.remove(FileCorre2)
        OutName=FileCorre

    if m!=0 and lineCount!=0:

        shutil.copy(os.path.join('folder', FileCorre), folderName)
        os.remove(FileCorre)
        
        OutName=FileCorre2[:-6]+'.raw'
        
##        print('From file ',FileCorre2[:-6]+'.raw',' ',m, 'rows were deleted for time jumps\n',)
##        print('A new file with the same name but finished as -corrected is created ')
        os.rename(FileCorre2,FileCorre2[:-6]+'.raw')
    print('In ',NameFile,' the ratio of correction is ',round(100.*(m+lineCount)/TotalLinesFile,2),'% where ',m+lineCount,' rows were deleted from ',TotalLinesFile,' rows')
##    print('total lines',TotalLinesFile,' and deleted ',m+lineCount)
    return OutName

np.warnings.filterwarnings('ignore')#to avoid the error messages



print('Insert the path where the raw are --for instance d:\Mrrdata/')
Root=input()  #input from the user 
os.chdir(Root)

folder=Root
dircf=glob.glob(Root+'*.raw')
dircf=np.sort(dircf)
print('In this folder there are '+str(len(dircf))+' raw files')
##print('The script generate netcdf file with the same name of raw files or include correted at the end\n')
print('A new file with the same name but finished as -corrected is created \n')

for name in dircf:
    
    NameFile=name 

    count=0
    NameFile=CorrectorFile(NameFile)
    
