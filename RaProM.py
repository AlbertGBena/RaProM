###script for process the raw data from Mrr radar of Mettek
###Autor: Albert Garcia Benadi
###ORCID: 0000-0002-5560-4392


import numpy as np
import calendar
import datetime
import time
import miepython as mp
from math import e
from netCDF4 import Dataset
import glob
import os
import sys
import shutil


def CorrectorFile(fid):
    NameFile=fid 

    FileCorre=NameFile[:-4]+'-corrected.raw' #create a new file

    

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
    f=open(NameFile,'r')

    for i in range(totallines):
        line=f.readline()
    
    
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
    if m==0:
        os.remove(FileCorre)
        OutName=NameFile
    else:
        shutil.copy(os.path.join('folder', NameFile), folderName)
        os.remove(NameFile)
        OutName=FileCorre
        print('A new file with the same name but finished as -corrected is created ')
        print('From file ',NameFile,' ',m, 'rows were deleted',)

    return OutName

    
            
    

    
    
def Promig(vector):
    rep=len(vector)
    Out=[] 
    for i in range(31):
        H=[]
        T=[]
        for j in range(rep):
            h1=vector[j]

            h2=h1[i]

            H.append(h2) 
    
        for j in range(64):
            S=[]
            for k in range(rep):
                h3=H[k]
            
                h4=h3[j]
            
                S.append(h4)
        
            NoNul=np.count_nonzero(S)
            if NoNul>len(S)/(100./Ocurrence):
                Sum=float(np.sum(S))/float(NoNul)
            else:
                Sum=np.nan
        
            T.append(Sum)
        Out.append(T)
    return np.asarray(Out)


def group(a,indexcentral,Nnan,d):
    
    d=np.asarray(d)
    a=np.asarray(a)
    b=np.where(a>=0)
    acut=np.asarray(a[64:128])
    bcut=np.where(acut>=0)
    

    c=b[0]
    ccut=bcut[0]

    if indexcentral<=80 or indexcentral>=112:

        for i in range(np.size(b)-1):
            if c[i]-indexcentral<=0 and c[i+1]-indexcentral>=0:
                index=c[i+1]
                break
            else:
                index=indexcentral
        
        cond=True
        cont=0
        incr1=0#starts at 0

        while cond:
            if cont>=Nnan or index+incr1>=len(a)-1:
                cond=False
            
            if np.isnan(a[index+incr1]):
                cont+=1
            else:
                cont=0


            incr1+=1

        cont=0;
        incr2=1#starts at 1
        
        cond=True
        while cond:
            if cont>=Nnan or index-incr2<=0:
                cond=False
            
            if np.isnan(a[index-incr2]):
                cont+=1
            
            else:
                cont=0


            incr2+=1
        vf2=np.copy(a);xf2=np.copy(d)
        vf2[0:index-incr2+1]=np.nan
        vf2[index+incr1:]=np.nan
        xf2[0:index-incr2+1]=np.nan
        xf2[index+incr1:]=np.nan
            

    else:
        
        for i in range(np.size(bcut)-1):
            if ccut[i]-indexcentral<=0 and ccut[i+1]-indexcentral>=0:
                index=ccut[i+1]
                break
            else:
                index=indexcentral-64
        
        cond=True
        cont=0
        incr1=0#starts at 0
    
        while cond:
            if cont>=Nnan or index+incr1>=len(acut)-1:
                cond=False
            
            if np.isnan(acut[index+incr1]):
                cont+=1
            else:
                cont=0


            incr1+=1

        cont=0;
        incr2=1#starts at 1
        
        cond=True
        while cond:
            if cont>=Nnan or index-incr2<=0:
                cond=False
            
            if np.isnan(acut[index-incr2]):
                cont+=1
            
            else:
                cont=0


            incr2+=1
        blanckv=np.nan*np.ones(len(acut))
        vf1=np.copy(acut);xf1=np.copy(d[64:128])
        vf1[0:index-incr2+1]=np.nan
        vf1[index+incr1:]=np.nan
        vff1=np.concatenate((blanckv,vf1))
        vf2=np.concatenate((vff1,blanckv))
                               
        xf1[0:index-incr2+1]=np.nan
        xf1[index+incr1:]=np.nan
        xff1=np.concatenate((blanckv,xf1))
        xf2=np.concatenate((xff1,blanckv))
        


    return vf2,xf2
    
    



def Process(matrix,he,temps,D):#This is the core from the preocessing
    neta,etan,etaNdb=HildrenS(matrix)
    
    Etan=FindRealPeaks(etan)
    etaN=np.multiply(Etan,Cte)
    etaV=etaN/Deltav#convert eta(n) in eta(v)
    state=[]#get the values 10 to water,-10 to snow and 0 if it is impossible to
    
    zewater=[];Ni=[];VT=[];Z=[];Z_da=[];Vhail=[]
    Noise=np.multiply(neta,Cte)
    
    for m in range(len(etaV)): 
        
        proba=np.where(~np.isnan(etaN[m]))
        leN=len(etaN[m])
                ##                DELETE THE SPORADIC VALUES IN THE LAST HEIGTHS
        if len(proba[0])!=0 and m==30:
            tes1=np.where(~np.isnan(etaN[m-1]))
            tes2=np.where(~np.isnan(etaN[m-2]))
                    
            if len(tes1[0])==0 or len(tes2[0])==0:
                etaV[m]=np.ones(leN)*np.nan
                etaN[m]=np.ones(leN)*np.nan
        

        zewater.append(10**18*lamb**4*np.nansum(etaV[m])*Deltav/K2w)#Calculate by its definition of equivalent reflectivy
        nde=[];vt=[];velHail=[]
        for n in range(len(etaV[0])):#Calculate the Ze from every gate without PIA                    
            value=6.18*etaV[m][n]*dv[m]*e**(-1*0.6*D[m][n])
            
            value3=(9.65-10.3*e**(-1*0.6*D[m][n]))*dv[m]

            velHail.append(13.96*np.sqrt(10*D[m][n]))#vel from Hail Ulbrich and atlas 1982
            
            sbk=SigmaScatt[m][n]
            
            value2=(10**6)*(value/sbk)#N in m-3 mm-1
            
            nde.append(value2)#units mm-1 m-3
            
            vt.append(value3)#terminal speed in function heigh and diameter
        Vhail.append(velHail)    
        VT.append(vt)    
        Ni.append(nde)
        
        
    z,lwc,rr,ze=Parameters(Ni,D,VT,0)
    
    #stratiform case (M-P)
    vwaterR=2.65*np.power(zewater,.114)#values a,b from Atlas et al. 1973
    vsnowR=.817*np.power(zewater,.063)#values a,b from Atlas et al. 1973
    vwaterMiestr=2.65*np.power(ze,.114)
    #Thunderstorm Rain (S-S)
    
    vwaterMieconv=4.13*np.power(ze,.062)

    vwaterMie=np.nanmean([vwaterMieconv,vwaterMiestr],axis=0)
    
    speeddeal=np.arange(-64*fNy,2*64*fNy,fNy)
    NewM=[];state=[];mov=[];VerTur=[]
    W=[];Sig=[];Sk=[];lwc=[];rr=[];Z_da=[];SnowRate=[];N_da=[];Snr=[];Kurt=[];dm=[];nw=[]
    if np.sum(np.nanmean([vwaterR,vsnowR],axis=0))!=0:
        
        limitValueDeal=4#Limit value to activate the dealiasing
        DealMatrix=[]
        
        for o in range(len(etan)): 
            Snr.append(10*np.log10(1+((np.nanmax(etan[o])/neta[o]))))
                
            if o==0 or o==len(etan)-1:
                if o==0:
                    N_deal=np.ones(len(etaN[o]))*np.nan
                    n1=np.concatenate((N_deal,etaN[o]),axis=None)
                    etaN_da=np.concatenate((n1,etaN[o+1]),axis=None)
                if o==len(etan)-1:
                    N_deal=np.ones(len(etaN[o]))*np.nan
                    n1=np.concatenate((etaN[o],N_deal),axis=None)
                    etaN_da=np.concatenate((etaN[o-1],n1),axis=None)
                        
            else:
                N_deal=etaN[o-1]
                n1=np.concatenate((N_deal,etaN[o]),axis=None)
                etaN_da=np.concatenate((n1,etaN[o+1]),axis=None)

            DealMatrix.append(etaN_da) 

    ##            LOOP TO THE PRECIPITATION TYPE ESTIMATION
            Indvel=speeddeal*(etaN_da/etaN_da)
                
            
            av=np.where(etaN_da>=0.)
            av2=np.where(etaN[o]>=0.)

            
            if np.size(av)==0 or np.size(av2)==0 or np.size(av2)==64:
                
                ReVect=etaN_da*np.nan
                INewV=Indvel*np.nan
            else:
                I=np.nanargmax(np.asarray(etaN[o]))+64
                

                ReVect,INewV=group(etaN_da,I,5,Indvel)#function to find the group of values




            if np.isnan(ReVect).all():
                S=np.nan
                L=np.nan
                sigma3=np.nan
            else:
                if INewV[np.nanargmax(ReVect)]<0:
                    S=.5
                    L=5.

                else:
                                        
                    PT3=np.nansum(ReVect)
                    w3=np.nansum(np.prod([ReVect,speeddeal],axis=0))/PT3#estimated velocity
                    sigma3=np.sqrt(np.nansum(np.prod([ReVect,np.power(speeddeal-w3,2)],axis=0))/PT3)# spectral witdh
                    
                    S=(INewV[np.nanargmax(ReVect)]-vsnowR[o])
                    L=(INewV[np.nanargmax(ReVect)]-vwaterMie[o])
                    comment='nothing'
                    

            

                

            if abs(S)<=(2*abs(sigma3)) and abs(L)>(2*abs(sigma3)):#case not liquid, possible snow

                state.append(-10)
                if np.nanmin(INewV)<0 or np.nanmax(INewV)>12:
                    if np.nanmin(INewV)<0:
                        mov.append(-1)
                    else:
                        mov.append(1)
                else:
                    mov.append(np.nan)
                    
            if abs(S)>(2*abs(sigma3)) and abs(L)<=(2*abs(sigma3)):#case liquid

                state.append(10)
                if np.nanmin(INewV)<0 or np.nanmax(INewV)>12:
                    if np.nanmin(INewV)<0:
                        mov.append(-1)
                    else:
                        mov.append(1)
                else:
                    mov.append(np.nan)
            if np.isnan(L) and np.isnan(S):
                state.append(np.nan)
                if np.nanmin(INewV)<0 or np.nanmax(INewV)>12:
                    if np.nanmin(INewV)<0:
                        mov.append(-1)
                    else:
                        mov.append(1)
                else:
                    mov.append(np.nan)
            if abs(S)==abs(L) or (abs(L)<=(2*abs(sigma3)) and abs(S)<=(2*abs(sigma3))):#case mixed
                if INewV[np.nanargmax(ReVect)]<vwaterMie[o] and INewV[np.nanargmax(ReVect)]> vsnowR[o]:
                    state.append(0)#cas mixed
                else:
                    if INewV[np.nanargmax(ReVect)]< vsnowR[o]:
                        state.append(-10)#cas snow
                    if INewV[np.nanargmax(ReVect)]> vwaterMie[o]:
                        state.append(10)#cas rain
                        
                    
                if np.nanmin(INewV)<0 or np.nanmax(INewV)>12:
                    if np.nanmin(INewV)<0:
                        mov.append(-1)
                    else:
                        mov.append(1)
                else:
                    mov.append(np.nan)

            if np.isnan(S) and ~np.isnan(L):#case liquid, but possible wrong election
                state.append(10)
                if abs(L)>=limitValueDeal:
                    
                    LonCut=len(etaN[o])
                    n1=np.ones(LonCut)*np.nan
                    
                    n2=etaN_da[LonCut:LonCut*2]
                    nc1=np.concatenate((n1,n2))
                    ReVect=np.concatenate((nc1,n1))#vector corrected
                    mov.append(1)
                else:
                    mov.append(np.nan)
                    

                
            if ~np.isnan(S) and np.isnan(L):#case not liquid, but possible wrong election
                state.append(-10)

                if S>=limitValueDeal:
                    
                    LonCut=len(etaN[o])/2
                    n1=np.ones(LonCut)*np.nan
                    n3=np.ones(3*LonCut)*np.nan
                    n2=etaN_da[LonCut:LonCut+len(etaN[o])]
                    nc1=np.concatenate((n1,n2))
                    ReVect=np.concatenate((nc1,n3))#vector corrected
                    mov.append(-1)
                else:
                    mov.append(np.nan)

                
                
                    

            if abs(L)>(2*abs(sigma3)) and abs(S)>(2*abs(sigma3)):#difficult case

                if abs(L)>=limitValueDeal and abs(S)>=limitValueDeal:#wrong election
                    if L-limitValueDeal<=0 or S-limitValueDeal<=0:#v expected is bigger than the found. Shift the vector from 128-32 to 128+32

                        if np.isnan(etaN_da[96:160]).all():
                            ReVect=np.ones(len(etaN_da))*np.nan
                        else:
                            I=np.nanargmax(np.asarray(etaN_da[96:160]))+96
                            ReVect,INewV=group(etaN_da,I,5,Indvel)

                        
                        if np.isnan(ReVect).all():
                            S1=np.nan
                            L1=np.nan
                            mov.append(np.nan)
                        else:
                            mov.append(1)
                            S1=(INewV[np.nanargmax(ReVect)]-vsnowR[o])
                            L1=(INewV[np.nanargmax(ReVect)]-vwaterMie[o])
                            PT4=np.nansum(ReVect)
                            w4=np.nansum(np.prod([ReVect,speeddeal],axis=0))/PT4#estimated velocity
                            sigma4=np.sqrt(np.nansum(np.prod([ReVect,np.power(speeddeal-w4,2)],axis=0))/PT4)# spectral witdh
                        
                            
                            if (abs(L1)<=(2*abs(sigma4)) and abs(S1)>(2*abs(sigma4))):
                                state.append(10)
                            if ~np.isnan(S1) and abs(L1)>(2*abs(sigma4)):
                                state.append(-10)
                            if abs(S1)<=(2*abs(sigma4)) and abs(L1)<=(2*abs(sigma4)):
                                state.append(0)
                            

                            if np.isnan(S1) and ~np.isnan(L1):
                                state.append(10)#liquid
                            if np.isnan(L1) and ~np.isnan(S1):
                                state.append(-10)#not liquid

                        if np.isnan(S1) and np.isnan(L1):
                            state.append(np.nan)
                            
                    if L-limitValueDeal>0 or S-limitValueDeal>0:#v expected is lower than the found. Shift the vector from 64-32 to 64+32




                        if np.isnan(etaN_da[32:96]).all():
                            ReVect=np.ones(len(etaN_da))*np.nan
                        else:
                            I=np.nanargmax(np.asarray(etaN_da[32:96]))+32
                            ReVect,INewV=group(etaN_da,I,5,Indvel)

                            
                        if np.isnan(ReVect).all():
                            S1=np.nan
                            L1=np.nan
                            mov.append(np.nan)
                        else:
                            mov.append(-1)
                            S1=(INewV[np.nanargmax(ReVect)]-vsnowR[o])
                            L1=(INewV[np.nanargmax(ReVect)]-vwaterMie[o])
                            PT5=np.nansum(ReVect)
                            w5=np.nansum(np.prod([ReVect,speeddeal],axis=0))/PT5#estimated velocity
                            sigma5=np.sqrt(np.nansum(np.prod([ReVect,np.power(speeddeal-w5,2)],axis=0))/PT5)# spectral witdh
                        

                            if (abs(L1)<=(2*abs(sigma5)) and abs(S1)>(2*abs(sigma5))):
                                state.append(10)
                            if ~np.isnan(S1) and abs(L1)>(2*abs(sigma5)):
                                state.append(-10)
                            if abs(S1)<=(2*abs(sigma5)) and abs(L1)<=(2*abs(sigma5)):
                                state.append(0)
                            

                            if np.isnan(S1) and ~np.isnan(L1):
                                state.append(10)#liquid
                            if np.isnan(L1) and ~np.isnan(S1):
                                state.append(-10)#not liquid

                        if np.isnan(S1) and np.isnan(L1):
                            state.append(np.nan)
                if (abs(L)-limitValueDeal<0) and (abs(S)-limitValueDeal<0):
                    if abs(S)>abs(L):
                        state.append(10)
                        mov.append(np.nan)
                    if abs(S)<abs(L):
                        state.append(-10)
                        mov.append(np.nan)
                    if abs(S)==abs(L):
                        state.append(0)
                        mov.append(np.nan)
                if (abs(L)-limitValueDeal<=0) and (abs(S)>limitValueDeal):
                    
                    state.append(10)
                    mov.append(np.nan)
                if (abs(S)-limitValueDeal<=0) and (abs(L)>limitValueDeal):
                    
                    state.append(-10)
                    mov.append(np.nan)
                    
                        
                    
            NewM.append(ReVect)
        

        valuemax=[]
        for m in range(len(NewM)): 
               

            iden=np.where(NewM[m]>=0.)
    
            if np.size(iden)==0:
                valuemax.append(np.nan)
            else:
                valuemax.append(speeddeal[np.nanargmax(NewM[m])])
     

        Dife=np.diff(valuemax)
        
        indole=np.argwhere(abs(np.diff(valuemax))>8.)
        Inindole=np.copy(indole)

            
        CountIndole=[]
        if len(indole)!=0:
            CountIndole.append(int(indole[0]))
        while len(indole)!=0:
                
            Vector=DealMatrix[int(indole[0]+1)]

            indx=np.nanargmax(NewM[int(indole[0]+1)])

                
            Indvel=speeddeal*(Vector/Vector)
                
            if Dife[int(indole[0])]>0:
                if indx>(64+32):
                    IDX=indx-64
                else:
                    IDX=indx
                CoVector,daig=group(Vector,IDX,5,Indvel)#function to find the group of values
                NewM[int(indole[0]+1)]=CoVector
                    
            else:
                if indx<(64+32):
                    IDX=indx+32
                else:
                    IDX=indx+64
                CoVector,daig=group(Vector,IDX,5,Indvel)#function to find the group of values
                NewM[int(indole[0]+1)]=CoVector
            newZe=10**18*lamb**4*np.nansum(CoVector)*Deltav/K2w
            NvwaterR=2.65*np.power(newZe,.114)#values a,b from Atlas et al. 1973
            NvsnowR=.817*np.power(newZe,.063)#values a,b from Atlas et al. 1973
            if np.isnan(CoVector).all():
                S=np.nan
                L=np.nan
            else:
                S=(speeddeal[np.nanargmax(CoVector)]-NvsnowR)
                L=(speeddeal[np.nanargmax(CoVector)]-NvwaterR)
                PT2=np.nansum(CoVector)
                w2=np.nansum(np.prod([CoVector,speeddeal],axis=0))/PT2#estimated velocity
                sigma2=np.sqrt(np.nansum(np.prod([CoVector,np.power(speeddeal-w2,2)],axis=0))/PT2)# spectral witdh

                
                if abs(L)<=(abs(sigma2)) and abs(S)>(abs(sigma2)):
                    state[int(indole[0]+1)]=10.
                if abs(S)<=(abs(sigma2)) and abs(L)>(abs(sigma2)):
                    state[int(indole[0]+1)]=-10.
                    comment='snow'
                if abs(L)>(abs(sigma2))and abs(S)>(abs(sigma2)):
                    state[int(indole[0]+1)]=20.
                    comment='unkown'
                if abs(L)<=(abs(sigma2)) and abs(S)<=(abs(sigma2)):
                    if w2>NvsnowR and w2<NvwaterR:
                        state[int(indole[0]+1)]=0.
                        comment='mixed'
                    else:
                        if w2<NvsnowR:
                            state[int(indole[0]+1)]=-10.
                            comment='snow'
                        if w2 > NvwaterR:
                            state[int(indole[0]+1)]=+10.
                            comment='rain'
                        

            valuemax=[]
            for m in range(len(NewM)): 
                   
   
                iden=np.where(NewM[m]>=0.)
   
                if np.size(iden)==0:
                    valuemax.append(np.nan)
                else:
                    valuemax.append(speeddeal[np.nanargmax(NewM[m])])
                    

                
            indole=np.argwhere(abs(np.diff(valuemax))>8.)
            Dife=np.diff(valuemax)
            if len(indole)!=0:

                if int(indole[0])==CountIndole[-1]:
                    indole=[]
                else:
                    CountIndole.append(int(indole[0]))
        

        for m in range(len(state)):
            if m!=0 and m!=len(state)-1:
                s1=state[m-1]
                s2=state[m]
                s3=state[m+1]
                if s2==0 and s1==-10 and s3==-10:
                    state[m]=-10
                if s2==0 and s1==10 and s3==10:
                    state[m]=10
                if s2==20 and s1==10 and s3==10:
                    state[m]=10
                if s2==20 and s1==-10 and s3==-10:
                    state[m]=-10
        Mwater=[]
        Msnow=[]
        Mmixed=[]
        Mhail=[]
        MDriz=[]
        Munk=[]
        ZE=[]

        roW=10**6 #water density g/m3

        vel=np.copy(speeddeal)
        Nde=[]
        PIA=[]

        for m in range(len(NewM)):
            vector=NewM[m]
            

            PT=np.nansum(vector)
            w=np.nansum(np.prod([vector,vel],axis=0))/PT#estimated velocity
            sigma=np.sqrt(np.nansum(np.prod([vector,np.power(vel-w,2)],axis=0))/PT)# spectral witdh
            sk=np.nansum(np.prod([vector,np.power(vel-w,3)],axis=0))/(PT*pow(sigma,3))# spectral witdh
            Kur=np.nansum(np.prod([vector,np.power(vel-w,4)],axis=0))/(PT*pow(sigma,4))# Kurtossis
            ValueZe=(10**18*lamb**4*np.nansum(NewM[m]))/(np.pi**5*K2w)
            if ValueZe<=0 or np.isnan(ValueZe):
                ZE.append(np.nan)
            else:
                ZE.append(10*np.log10(ValueZe))
            if w==0.:
                w=np.nan

            
            W.append(w)
            
            Sig.append(sigma)
            Sk.append(sk)
            Kurt.append(Kur)

            if state[m]==0:#Mixed case
                Mmixed.append(NewM[m])
                Value=10**18*lamb**4*np.nansum(NewM[m])/(np.pi**5*K2w)
                
                Z_da.append(np.nan)
                lwc.append(np.nan)
                dm.append(np.nan)
                nw.append(np.nan)
                rr.append(np.nan)
                N_da.append(np.nan)
                PIA.append(np.nan)

                SnowRate.append(np.nan)
                Nde.append(np.ones(len(etaN[1]))*np.nan)
                if mov[m]==-1:#case mixed and upward
                      VerTur.append((.817*np.power(ValueZe,.063))-w)
                if mov[m]==1:#case mixed and downward
                      VerTur.append(w-(.817*np.power(ValueZe,.063)))
                if np.isnan(mov[m]):
                      VerTur.append(np.nan)
                            
            else:
                Mmixed.append(NewM[m]*np.nan)
            
            
            if state[m]==10:#rain case
                
                Mwater.append(NewM[m])
                SnowRate.append(np.nan)

                dif=[]#diference between diameters for N
                dif2=[]#diference between diameters for Z
                nde=[]
                indexFinded=[]
                for n in range(len(D[m])):
                    if n==0 or n==len(D[m])-1:
                        if n==0:
                            dif2.append(D[m][n+1]-D[m][n])
                            dif.append(D[m][n+1]-D[m][n])
                        if n==len(D[m])-1:
                            dif.append(abs(D[m][n-1]-D[m][n]))
                            dif2.append(abs(D[m][n]-D[m][n-1]))
                    else:
                        dif2.append(abs((D[m][n+1]-D[m][n])))
                        dif.append(abs((D[m][n+1]-D[m][n-1]))/2.)
                        condFH=speed[n]-w
                        
                        if condFH>=0:
                            indexFinded.append(n)
                        
                
                    EtaV=NewM[m][64:64*2]#interval for water choosed
                    value=6.18*EtaV[n]*dv[m]*e**(-1*0.6*D[m][n])
                    s=SigmaScatt[m][n]
                    value2=(10**6)*(value/s)#N in m-3 mm-1
                    nde.append(value2)#units mm-1 m-3
                
                #Calculate the diamater from the mean vel found
                diaWork=D[m]
                if w<=0 or w>=11:
                    if w<=0:#cas snow
                        diamHail=1
                    if w>=11:#case hail
                        diamHail=6
                else:
                    if len(indexFinded)==0:
                        diamHail=3#case liquid
                    else:
                        
                        diamHail=diaWork[indexFinded[0]]
                
                    
                
                ##                APPLY THE ATTENUATTION
                if m==0:
                    PIA.append(1)
                if m!=0 and len(PIA)==0:
                    PIA.append(np.ones(m+1))
                    
                Np=np.multiply(nde,PIA[-1])
                Pro=[]
                for k in range(len(Np)):
                    pro=SigmaExt[m][k]*Np[k]*dif[k]
            
                    Pro.append(pro)
           
                kp=np.nansum(Pro)*10**-6
        
                num=2*kp*he[m]
                N=-np.multiply(Np,np.log(1-num)/num)
                Pro2=[]
                for k in range(len(N)):
                    pro2=SigmaExt[m][k]*N[k]*dif[k]
                    Pro2.append(pro2)
            
                Kr=np.nansum(Pro2)*10**-6
                pia=PIA[-1]*e**(2*Kr*he[m])
                if pia>=10. or num==0.:
                    pia=10.
                    LastN=nde
                else:
                    LastN=N
                PIA.append(pia)
                

                Nde.append(LastN)
                value=np.nansum(np.prod([np.power(D[m],6),LastN,dif2],axis=0))
                value2=np.nansum(np.prod([np.power(D[m],3),LastN,dif2],axis=0))
                value3=np.nansum(np.prod([np.power(D[m],3),LastN,dif2,w],axis=0))
                value4=np.nansum(np.prod([np.power(D[m],4),LastN,dif2],axis=0))
                if np.nansum(nde)<=0.:
                    N_da.append(np.nan)
                else:
                    N_da.append(np.log10(np.nansum(LastN)))

                if diamHail>=5:#Hail case
                    
                    Mhail.append(NewM[m])
                    state[m]=-20.
                    
                else:
                    Mhail.append(NewM[m]*np.nan)

                if diamHail<0.5:#drizzle case
                    
                    MDriz.append(NewM[m])
                    state[m]=5.
                    
                else:
                    MDriz.append(NewM[m]*np.nan)
                        
                        
                        
                if value<=0. or np.isnan(value):
                    Z_da.append(np.nan)
                else:
                    Z_da.append(10*np.log10(value))
                if value2==0.:
                    lwc.append(np.nan)
                else:
                    lwc.append(roW*value2*(np.pi/6.)*(10**-9))
                if value3==0.:
                    rr.append(np.nan)
                else:
                    rr.append(value3*(np.pi/6.)*(10**-9)*1000.*3600.)
                if value4==0:
                    dm.append(np.nan)
                    nw.append(np.nan)
                else:
                    dm.append(value4/value2)
                    nw.append(np.log10(256.*(roW*value2*(np.pi/6.)*(10**-9))/ (np.pi*roW*(value4/value2)**4)))
                if mov[m]==-1:#case rain and upward
                      VerTur.append((2.6*np.power(ValueZe,.107))-w)
                if mov[m]==1:#case rain and downward
                      VerTur.append(w-(2.6*np.power(ValueZe,.107)))
                if np.isnan(mov[m]):#case rain and upward
                      VerTur.append(np.nan)
               
            else:
                Mwater.append(NewM[m]*np.nan)
                
                
            if state[m]==-10:#Snow case
                
                Msnow.append(NewM[m])


                if ValueZe<=0 or np.isnan(ValueZe):
                    
                    SnowRate.append(np.nan)
                else:
                    
                    SnowRate.append(np.power(ValueZe/56.,1/1.2))#following Matrosov (2007) constants - https://link.springer.com/article/10.1007/s00703-011-0142-z#CR15
                Z_da.append(np.nan)
                lwc.append(np.nan)
                nw.append(np.nan)
                dm.append(np.nan)
                rr.append(np.nan)
                N_da.append(np.nan)
                PIA.append(np.nan)
                Nde.append(np.ones(len(etaN[1]))*np.nan)
                if mov[m]==-1:#case snow and upward
                      VerTur.append((.817*np.power(ValueZe,.063))-w)
                if mov[m]==1:#case snow and downward
                      VerTur.append(w-(.817*np.power(ValueZe,.063)))
                if np.isnan(mov[m]):
                      VerTur.append(np.nan)
                    
                
                
            else:
                Msnow.append(NewM[m]*np.nan)
                
            

            if np.isnan(state[m]):
                Mwater.append(NewM[m]*np.nan)
                Msnow.append(NewM[m]*np.nan)
                Mmixed.append(NewM[m]*np.nan)
                Mhail.append(NewM[m]*np.nan)
                MDriz.append(NewM[m]*np.nan)
                Munk.append(NewM[m]*np.nan)
                Z_da.append(np.nan)
                lwc.append(np.nan)
                nw.append(np.nan)
                dm.append(np.nan)
                rr.append(np.nan)
                N_da.append(np.nan)
                SnowRate.append(np.nan)
                Nde.append(np.ones(len(etaN[1]))*np.nan)
                VerTur.append(np.nan)
                PIA.append(np.nan)

                
            if state[m]==20:#cas unknown
                Munk.append(NewM[m])
                
                
                Value=10**18*lamb**4*np.nansum(NewM[m])/(np.pi**5*K2w)#use the rayleight estimation
                
                Z_da.append(np.nan)
                lwc.append(np.nan)
                nw.append(np.nan)
                dm.append(np.nan)
                rr.append(np.nan)
                N_da.append(np.nan)
                PIA.append(np.nan)

                SnowRate.append(np.nan)
                Nde.append(np.ones(len(etaN[1]))*np.nan)
                if mov[m]==-1:#case hail and upward, using the same for snow
                      VerTur.append((.817*np.power(ValueZe,.063))-w)
                if mov[m]==1:#case hail and downward, using the same for snow
                      VerTur.append(w-(.817*np.power(ValueZe,.063)))
                if np.isnan(mov[m]):
                      VerTur.append(np.nan)
                            
            else:
                Munk.append(NewM[m]*np.nan)
                            
        for m in range(len(state)):
            if m!=0 and m!=len(state)-1:
                s1=state[m-1]
                s2=state[m]
                s3=state[m+1]
                if s2==-20 and s1==10 and s3==10:
                    state[m]=10        
                
                

        

    else:#There is not Signal
        
        blanck=np.nan*np.ones(len(etaN))
        state=blanck
        NewM=(etaN*np.nan)
        Z_da=blanck
        lwc=blanck
        nw=blanck
        dm=blanck
        rr=blanck
        N_da=blanck
        W=blanck
        Sig=blanck
        Sk=blanck
        Kurt=blanck
        SnowRate=blanck
        Nde=(etaN*np.nan)
        ZE=blanck
        mov=blanck
        VerTur=blanck
        Snr=blanck
        PIA=np.nan*np.ones(len(etaN)+1)
    if len(PIA)==31:
        PIA.append(np.nan)
    return state,NewM,Z_da,lwc,rr,SnowRate,W,Sig,Sk,Noise,N_da,Nde,ZE,mov,VerTur,Snr,Kurt,PIA,nw,dm

       
    
def group_consecutives(vals, step=1):
    """Return list of consecutive lists of numbers from vals (number list)."""
    run = []
    result = [run]
    expect = None
    for v in vals:
        if (v == expect) or (expect is None):
            run.append(v)
        else:
            run = [v]
            result.append(run)
        expect = v + step
    return result

def Parameters(n,d,v,da):#the diffrences between diameter aren't constant
    Z=[];lwc=[];rr=[];ze=[]
    
    roW=10**6 #water density g/m3
    for i in range(len(n)):#entrem en una altura
        D=d[i]
        N=n[i]
        w=v[i]
        dif=[]
        
        
        if da==1:#estem amb l'eix dealiased
            for m in range(len(D)):#creo el vector de diferencies
                if m==0 or m==len(D)-1:
                    if m==0:
                        dif.append(d[i][m+1]-d[i][m])
                    if m==len(D)-1:
                        dif.append(abs(d[i][m-1]-d[i][m]))
                else:
                    if m<len(w)/2.:
                        dif.append((d[i][m+1]-d[i][m]))
                    else:
                        dif.append(abs((d[i][m]-d[i][m+1])))

        else:
            
            for m in range(len(D)):#creo el vector de diferencies

                if m==0 or m==len(D)-1:
                    if m==0:
                        dif.append(d[i][m+1]-d[i][m])
                    if m==len(D)-1:
                        dif.append(d[i][m]-d[i][m-1])
                else:
                    dif.append((d[i][m+1]-d[i][m]))
        value=np.nansum(np.prod([np.power(D,6),N,dif],axis=0))
        value2=np.nansum(np.prod([np.power(D,3),N,dif],axis=0))
        value3=np.nansum(np.prod([np.power(D,3),N,dif,w],axis=0))
        if value==0.:
            Z.append(np.nan)
            ze.append(np.nan)
        else:
            Z.append(10*np.log10(value))
            ze.append(value)
        if value2==0.:
            lwc.append(np.nan)
        else:
            lwc.append(roW*value2*(np.pi/6.)*(10**-9))
        if value3==0.:
            rr.append(np.nan)
        else:
            rr.append(value3*(np.pi/6.)*(10**-9)*1000.*3600.)
        
    return Z,lwc,rr,ze

def FindRealPeaks(matrix):#function to detect real peaks, wherein a peak has a minuim 3 consecutives values
    Matrix=[]
    for i in range(len(matrix)):
        vector=matrix[i]
        Vector=np.ones(len(vector))*np.nan
        Ind=np.argwhere(~np.isnan(vector))

        out=group_consecutives(Ind)

        Out=[]
        for j in range(len(out)):
            if len(out[j])>=3:
                Out.append(out[j])
        
        if len(Out)>0:
            for j in range(len(Out)):
                value=Out[j]
                for k in range(len(value)):
                    np.put(Vector,value[k],vector[value[k]])
        
            
        Matrix.append(Vector)
                    
        
                
        
    return Matrix

def ScatExt(diameter,longW):#for 1 height gate
    ag_lam=longW*1000.#mm Convert lamb from m to mm

    ag_mre=6.417
    ag_mim=2.758
    m = ag_mre + 1.0j * ag_mim
    
    scatt=[];extinct=[]
    for i in range(len(diameter)):
        r=diameter[i]/2.

        if np.isnan(r):
            scatt.append(np.nan)
            extinct.append(np.nan)

        else:
            
            x = 2*np.pi*r/ag_lam;#is non dimension, so the ag_lam and r have the same units
            qext, qsca, qback, g = mp.mie(m,x)
            absorb  = (qext - qsca) * np.pi * r**2
            scatt.append(qsca * np.pi * r**2)
            extinct.append(qext* np.pi * r**2)
    return scatt,extinct
            

        




def HildrenS(matrix):
    Pot=matrix
    
    PotHS=[]
    Noise=[]
    PotWithOutHS=[]
    for j in range(len(Pot)):
        v=Pot[j]#is the spectrum bins for one height
        v2=np.ones(len(v))*np.nan#create the vector result
        meanv=np.nanmean(v)
        varv=np.nanvar(v)
        if (meanv**2/(varv))>58 or varv==0:#bad signal
            v2=np.ones(len(v))*np.nan
            soroll=np.nan
            PotWithOutHS.append(v2)
        else:
                
            indMax=Peak(v)
            VMAX=[]
            for m in range(len(indMax)):
                Vmax=v[indMax[m]]
                VMAX.append(Vmax)

            maxim=np.argmax(VMAX)
        
            
            v[v>v[indMax[maxim]]]=np.nan
            PotWithOutHS.append(v)
            
            condition=True
            while condition:
                if (np.power(np.nanmean(v),2)/np.nanvar(v))>58 or np.nanvar(v)==0.:
                    soroll=np.nanmean(v)
                    condition=False

                np.put(v2,np.nanargmax(v),np.nanmax(v))
                np.put(v,np.nanargmax(v),np.nan)
        v2[v2<=1.2*soroll]=np.nan#Record the signal up 1.2 the noise level found

        Noise.append(soroll)
        PotHS.append(v2-soroll)
    return Noise,PotHS,PotWithOutHS


             





def date2unix(date):
    return calendar.timegm(date.timetuple())
def unix2date(unix):
    return datetime.datetime.utcfromtimestamp(unix)

def Peak(vector):
    InMax=[]
    
    L=len(vector)
    for i in range(L):
        j=i+1
        if j<(L-1):
            if vector[j]>=vector[j+1] and vector[j]>=vector[j-1]:
                InMax.append(j)
                
    
    return InMax 






velc=299792458.#light speed 
lamb=velc/(24.23*1e9)  #La frequency of radar is 24.23 GHz
ag_lam=lamb
fsampling=125000#Hz valor de la frecuecnia de muestreo
fNy=fsampling*lamb/(2*2*32*64) 
K2w=0.92
K2i=0.18
K2s=np.mean((K2w,K2i))
Deltaf=fsampling/(2*32*64) #ronda els 30 Hz
CetNtoetaV=2./(Deltaf*lamb)
Deltav=Deltaf*lamb/2.
Ocurrence=50.#value in % of ocurrence in the averaging matrix



np.warnings.filterwarnings('ignore')#to avoid the error messages



print('Insert the path where the raw are --for instance d:\Mrrdata/')
Root=raw_input()  #input from the user 
os.chdir(Root)

print('Insert the number of seconds for integration (usually 60 seconds)')
IntTime=raw_input()
IntTime=int(IntTime)

folder=Root
dircf=glob.glob(Root+'*.raw')
dircf=np.sort(dircf)
print('In this folder there are '+str(len(dircf))+' raw files')
print('The script generate netcdf file with the same name of raw files\n')

for name in dircf:
    
    NameFile=name 

    
    


    f=open(NameFile,'r')
    NumberLines=0
    while True:
        line=f.readline()

        line=line.strip()
        if line=='':
            break
        NumberLines+=1
    NumberSpectre=NumberLines/67.#divide by 67 because is the number of lines from spectrum


    
    f.close()

    



    count=0
    NameFile=CorrectorFile(NameFile)
    print('File in process  '+str(NameFile[:-4]))
    filenameplot=NameFile[:-4]+'-processed'

    f=open(NameFile,'r')
    a=f.readline()
    Hini=f.readline()
    f.close()
    Hini=Hini.strip()
    HIcolum=Hini.split()
    HIcolum=map(int,HIcolum[1:len(HIcolum)])#Get the height values and change to integer

    ##Found the parameters dv in function of the height (mrr physics equation)
    dv=[]
    for i in range(len(HIcolum)):
        if i>=1:
            dv.append(1+3.68*10**-5*HIcolum[i]+1.71*10**-9*HIcolum[i]**2)

    speed=np.arange(0,64*fNy,fNy)
    speed21=np.arange(0,32*fNy,fNy)
    speed22=np.arange(-32*fNy,0,fNy)
    speed2=np.concatenate((speed21,speed22),axis=0)#gis vector speed is to evaluate the upward
    


    ##Found the diameters in function of height and speed
    D=[]
    for i in range(len(dv)):
        d=[]
        for j in range(len(speed)):
            
            b=speed[j]/dv[i]
    ##        print(b,speed[j],dv[i])
            if b>=0.002 and b<=9.37:#Condition of diameter is good for 0.109 mm< D< 6 mm
                d.append(np.log((9.65-b)/10.3)*(-1/0.6))
            else:
                d.append(np.nan)
        D.append(d)#dimension 31 x 64 in mm


    dataset=Dataset(filenameplot+'.nc','w',format='NETCDF4')
    dataset.description='Data processed by MRR radar'
    dataset.author='Albert Garcia Benad'+u'\xed'
    dataset.orcid='0000-0002-5560-4392 '

    dataset.createDimension('DropSize',len(D[0]))

    dataset.createDimension('Height',len(HIcolum[1:]))

    dataset.createDimension('PIA_Height',len(HIcolum))

    dataset.createDimension('time',None)

    nc_times=dataset.createVariable('Time','float64',('time',))
    nc_ranges_H=dataset.createVariable('Height','f',('Height',))
    nc_ranges_H_PIA=dataset.createVariable('PIA_Height','f',('PIA_Height',))
    nc_ranges_DropSize=dataset.createVariable('DropSize','f',('DropSize',))

    nc_times.units = 'UNIX Time Stamp, SECOND SINCE 1970-01-01'
    nc_times.description='Time in unix format'

    nc_ranges_H.units = 'm'
    nc_ranges_H.description = 'Heights in meters a.g.l.'

    nc_ranges_H_PIA.units = 'm'
    nc_ranges_H_PIA.description = 'Heights in meters a.g.l.'

    nc_ranges_DropSize.units = 'mm'
    nc_ranges_DropSize.description = 'Size of the water drops'

    
    ##Found the scatter and extint sections in function of height and speed
    SigmaScatt=[]
    SigmaExt=[]
    for i in range(len(D)):#entry in height
        sig1,sig2=ScatExt(D[i],lamb)
        SigmaScatt.append(sig1)
        SigmaExt.append(sig2)
        




    #Start the script

    f=open(NameFile,'r')
    DataT=[]
    Mdades=[[]]
    o=0


    ########REFRACTION INDEX FROM WATERm=6.417+i*2.758 cited by Segelstein 1981
    ag_mre=6.417
    ag_mim=2.758
    Waterm = ag_mre + 1.0j * ag_mim


    #parametres que s'empran en el processat, que s'han de inicialitzar
    timeList=list()

    co=0


    PotCorrSum=[]


    PotSumWN=np.empty(shape=[32,64])
    PotSum=np.empty(shape=[31,64])#there is 32 height but the firts height is deleted
    

    Timecount=0
    TimeCounter=0
    Cont=0
    #START THE DATA READING
    print(datetime.datetime.now())
    ContPlot=0
    countwork=0
    

    while True:
        
        
        if countwork==7:
            sys.stdout.write('\b\b\b\b\b\b\b-------\r')
            countwork=0
        if countwork==0:
            sys.stdout.write('working')
            
        if countwork==6:
            sys.stdout.write('\b\b\b\b\b\b\bw-rking\r')
        
        if countwork==5:
            sys.stdout.write('\b\b\b\b\b\b\bwo-king\r')
        
        if countwork==4:
            sys.stdout.write('\b\b\b\b\b\b\bwor-ing\r')
        
        if countwork==3:
            sys.stdout.write('\b\b\b\b\b\b\bwork-ng\r')
        
        if countwork==2:
            sys.stdout.write('\b\b\b\b\b\b\bworki-g\r')
        
        if countwork==1:
            sys.stdout.write('\b\b\b\b\b\b\bworkin-\r')
        
        countwork+=1
            
        
        
        line=f.readline()

        line=line.strip()
        if line=='':
            if len(PotCorrSum)!=0:

                TimeCounter+=1

                timeVec=timeList[0]+(IntTime*TimeCounter)
                nc_times[Timecount:Timecount+1]=timeVec#add 1 second, beacuse the timestapmps from Mrr is the last second
    
    

                proeta=Promig(PotCorrSum)
                
                estat,NewMatrix,z_da,lwc,rr,SnowRate,w,sig,sk,Noi,DSD,NdE,Ze,Mov,velTur,snr,kur,PiA,NW,DM=Process(proeta,Hcolum[1:],timeVec,D)

                nc_state[Timecount,:]=np.array(np.ma.masked_invalid(estat),dtype='f')
                nc_w[Timecount,:]=np.array(np.ma.masked_invalid(w),dtype='f')
                nc_sig[Timecount,:]=np.array(np.ma.masked_invalid(sig),dtype='f')
                nc_sk[Timecount,:]=np.array(np.ma.masked_invalid(sk),dtype='f')
                nc_kur[Timecount,:]=np.array(np.ma.masked_invalid(kur),dtype='f')

                nc_PIA[Timecount,:]=np.array(np.ma.masked_invalid(PiA),dtype='f')
                
                
                nc_LWC[Timecount,:]=np.array(np.ma.masked_invalid(lwc),dtype='f')
                nc_RR[Timecount,:]=np.array(np.ma.masked_invalid(rr),dtype='f')

                nc_nw[Timecount,:]=np.array(np.ma.masked_invalid(NW),dtype='f')
                nc_dm[Timecount,:]=np.array(np.ma.masked_invalid(DM),dtype='f')
                
                nc_Z_da[Timecount,:]=np.array(np.ma.masked_invalid(z_da),dtype='f')
                nc_Z_e[Timecount,:]=np.array(np.ma.masked_invalid(Ze),dtype='f')
                nc_VerMov[Timecount,:]=np.array(np.ma.masked_invalid(Mov),dtype='f')

                nc_N_daTH[Timecount,:,:]=np.array(np.ma.masked_invalid(np.log10(NdE)),dtype='f')
                
                nc_SnowR[Timecount,:]=np.array(np.ma.masked_invalid(SnowRate),dtype='f')
                nc_Noi[Timecount,:]=np.array(np.ma.masked_invalid(Noi),dtype='f')
                nc_SNR[Timecount,:]=np.array(np.ma.masked_invalid(snr),dtype='f')
                nc_N_da[Timecount,:]=np.array(np.ma.masked_invalid(DSD),dtype='f')

                nc_VelTur[Timecount,:]=np.array(np.ma.masked_invalid(velTur),dtype='f')
               
            break


        
            
        columns=line.split()

        Date=columns[1] 

        
        dat = datetime.datetime(year = 2000+int(Date[0:2]), month = int(Date[2:4]), day = int(Date[4:6]), hour = int(Date[6:8]), minute = int(Date[8:10]), second = int(Date[10:12]))
        
        dat=int(date2unix(dat))

        timeList.append(dat-int(unix2date(dat).second))#fixs the time at 0 seconds
        TypeDate=columns[2] 
        DVS=columns[4] 
        DSN=columns[6] #serial number
        BW=columns[8] #Witdh band
        CC=int(columns[10]) #Calibration constant
        MDQ=columns[12:15] 
                        

        TypeFile=columns[16] 

        

        #Read the heigh parameters (second line from raw file)
        H=f.readline()
        H=H.strip()
        Hcolum=H.split()
        Hcolum=map(int,Hcolum[1:len(Hcolum)])#Get the height values and change to integer
        DeltaH=Hcolum[5]-Hcolum[4]#Height difference
        #
        #Read the tranference function (third line from raw file)
        FT=f.readline()
        FT=FT.strip()
        FTcolum=FT.split()
        FTcolum=map(float,FTcolum[1:len(FTcolum)])
        vectorV=np.arange(0,64*fNy,fNy)
        

        #constant value to conevrt F to f, include all constants
        Cte=DeltaH*float(CC)/(10**20)
        
        
        #Read the values F from file
        
        FQ=[]
        DataT=[]
        
        
        Number_bins=64#is the number of heights from file
        for j in range(Number_bins):
            Data=f.readline()
            Data=Data.strip()
            Data=Data.split()
            Dades=map(int,Data[1:len(Data)]) #extract the title for eac line F00, F01, etc
            
            FQ.append(Data[0])
            DataT.append(Dades)
            
        Pot=[]
        
        

        for k in range(len(Hcolum)):#the result is  amtrix crrected by transfer fuction and height, and dimension (LenHcolum)-1)
            COL=[row[k] for row in DataT]

            if k>=1:
                quo=FTcolum[k]/(k**2)
                Pot.append(np.divide(COL,quo))#pot is the spectra potence (v,i) except by muliply a constant where eac array is a column, or gate height


        hcor=np.asarray(Hcolum[1:])
        


        timeVec=timeList[0]+(IntTime*TimeCounter)
        if (dat-timeVec)>=IntTime:
            
            TimeCounter+=1
            timeVec=timeList[0]+(IntTime*TimeCounter)
            if len(PotCorrSum)==0:
                proeta=Pot*np.nan
            else:

                proeta=Promig(PotCorrSum)

            nc_times[Timecount:Timecount+1]=timeVec#add 1 second, beacuse the timestapmps from Mrr is the last second
            nc_ranges_H[:]=np.array(Hcolum[1:],dtype='f4')
            nc_ranges_H_PIA[:]=np.array(Hcolum,dtype='f4')
            nc_ranges_DropSize[:]=np.array(np.ma.masked_invalid(D[0]),dtype='f')
            ncShape2D = ('time','Height',)
            ncShape2D_PIA = ('time','PIA_Height',)#PIA starts in h=0, for this has 32 heights
            ncShape3D = ('time','Height','DropSize',)
            if Timecount==0:
                ##################create netcdf############

                nc_w=dataset.createVariable('W','f',ncShape2D)
                nc_w.description='vertical speed with aliasing correction'
                nc_w.units='m/s'
                

                nc_sig=dataset.createVariable('spectral width','f',ncShape2D)
                nc_sig.description='spectral width with dealiasing'
                nc_sig.units='m/s'
                

                nc_sk=dataset.createVariable('Skewness','f',ncShape2D)
                nc_sk.description='skewness with dealiasing'
                nc_sk.units='none'

                nc_kur=dataset.createVariable('Kurtosis','f',ncShape2D)
                nc_kur.description='Kurtosis with dealiasing'
                nc_kur.units='none'

                nc_PIA=dataset.createVariable('PIA','f',ncShape2D_PIA)
                nc_PIA.description='Path Integrated Attenuation'
                nc_PIA.units='none'

                nc_state=dataset.createVariable('Type','f',ncShape2D)
                nc_state.description='Indicate the type from hydrometeor as unknown(20), rain water (10), drizzle (5), mixed (0), snow(-10) and hail (-20) '
                nc_state.units=''
                

                

                nc_LWC=dataset.createVariable('LWC','f',ncShape2D)
                nc_LWC.description='Liquid water content'
                nc_LWC.units='g m-3'
                

                nc_RR=dataset.createVariable('RR','f',ncShape2D)
                nc_RR.description='Rain Rate'
                nc_RR.units='mm hr-1'

                nc_nw=dataset.createVariable('Nw','f',ncShape2D)
                nc_nw.description='Intercept parameter of the gamma distribution normalized to the liquid water content'
                nc_nw.units='log10(mm-1 m-3)'

                nc_dm=dataset.createVariable('Dm','f',ncShape2D)
                nc_dm.description='mean mas-weighted raindrop diameter'
                nc_dm.units='mm'

                nc_SnowR=dataset.createVariable('SR','f',ncShape2D)
                nc_SnowR.description='Snow Rate'
                nc_SnowR.units='mm hr-1'
                
                
                

                nc_Z_da=dataset.createVariable('Z','f',ncShape2D)
                nc_Z_da.description='Reflectivity considering only liquid drops'
                nc_Z_da.units='dBZ'

                nc_Z_e=dataset.createVariable('Ze','f',ncShape2D)
                nc_Z_e.description='Equivalent Reflectivity'
                nc_Z_e.units='dBZ'

                nc_VerMov=dataset.createVariable('Vmov','f',ncShape2D)
                nc_VerMov.description='Verical movement +1 downward -1 upward'
                nc_VerMov.units='None'
                
                nc_N_da=dataset.createVariable('N(D)','f',ncShape2D)
                nc_N_da.description='Drop Size Distribution'
                nc_N_da.units='log10(m-3 mm-1)'

                nc_N_daTH=dataset.createVariable('N(D) in function of time and height','f',ncShape3D)
                nc_N_daTH.description='Drop Size Distributionin function of time and height'
                nc_N_daTH.units='log10(m-3 mm-1)'
                
                nc_SNR=dataset.createVariable('SNR','f',ncShape2D)
                nc_SNR.description='Signal noise relation from signal without dealiasing'
                nc_SNR.units='dB'

                nc_Noi=dataset.createVariable('Noise','f',ncShape2D)
                nc_Noi.description='Noise'
                nc_Noi.units='m-1'

                nc_VelTur=dataset.createVariable('Fall speed variability','f',ncShape2D)
                nc_VelTur.description='Estimate the fall speed variability'
                nc_VelTur.units='m/s'


            PotCorrSum=[]#empty the array    
            estat,NewMatrix,z_da,Lwc,Rr,SnowRate,w,sig,sk,Noi,DSD,NdE,Ze,Mov,velTur,snr,kur,PiA,NW,DM=Process(proeta,Hcolum[1:],timeVec,D)

            nc_state[Timecount,:]=np.array(np.ma.masked_invalid(estat),dtype='f')
            nc_w[Timecount,:]=np.array(np.ma.masked_invalid(w),dtype='f')
            nc_sig[Timecount,:]=np.array(np.ma.masked_invalid(sig),dtype='f')
            nc_sk[Timecount,:]=np.array(np.ma.masked_invalid(sk),dtype='f')
            nc_kur[Timecount,:]=np.array(np.ma.masked_invalid(kur),dtype='f')
            nc_PIA[Timecount,:]=np.array(np.ma.masked_invalid(PiA),dtype='f')

            nc_N_daTH[Timecount,:,:]=np.array(np.ma.masked_invalid(np.log10(NdE)),dtype='f')
            
            
            nc_LWC[Timecount,:]=np.array(np.ma.masked_invalid(Lwc),dtype='f')
            nc_RR[Timecount,:]=np.array(np.ma.masked_invalid(Rr),dtype='f')

            nc_nw[Timecount,:]=np.array(np.ma.masked_invalid(NW),dtype='f')
            nc_dm[Timecount,:]=np.array(np.ma.masked_invalid(DM),dtype='f')
            
            nc_Z_da[Timecount,:]=np.array(np.ma.masked_invalid(z_da),dtype='f')
            nc_Z_e[Timecount,:]=np.array(np.ma.masked_invalid(Ze),dtype='f')
            nc_VerMov[Timecount,:]=np.array(np.ma.masked_invalid(Mov),dtype='f')
            
            nc_SnowR[Timecount,:]=np.array(np.ma.masked_invalid(SnowRate),dtype='f')
            nc_Noi[Timecount,:]=np.array(np.ma.masked_invalid(Noi),dtype='f')
            nc_SNR[Timecount,:]=np.array(np.ma.masked_invalid(snr),dtype='f')
            nc_N_da[Timecount,:]=np.array(np.ma.masked_invalid(DSD),dtype='f')

            nc_VelTur[Timecount,:]=np.array(np.ma.masked_invalid(velTur),dtype='f')

            Timecount=Timecount+1
            
            
        PotCorrSum.append(Pot)#add matrix

    f.close()
    print('\n')
    print(datetime.datetime.now())                    




    dataset.close()
        


            

            
                                
                                
            
