import numpy as np
from netCDF4 import Dataset
import netCDF4 as nc4
import matplotlib.pyplot as plt
import datetime
import calendar
import matplotlib.dates as mdates
import matplotlib as mpl
import os
import glob


def date2unix(date):
    return calendar.timegm(date.timetuple())
def unix2date(unix):
    return datetime.datetime.utcfromtimestamp(unix)



print('Insert the path where the raw are --for instance d:\Mrrdata/')
Root=raw_input()  #input from the user 
os.chdir(Root)
print('select the file without extension')
dircf=glob.glob(Root+'*.nc')
for i in dircf:
    print(i)
    print(' ')
    

NameF=raw_input()

NameFile=Root+NameF+'.nc'

f=nc4.Dataset(NameFile,'r')


ZE=f.variables['Ze'][:,:]
W=f.variables['W'][:,:]
S=f.variables['Type'][:,:]
ti=f.variables['Time'][:] 
He=f.variables['Height'][:]

f.close()

Zem=np.ma.masked_invalid(ZE)
Wm=np.ma.masked_invalid(W)
Sm=np.ma.masked_invalid(S)

t=ti-60.
row,col=np.shape(Zem)



Data=str(unix2date(ti[0]-60.))


#lines to convert the time array
xti2=np.array([unix2date(t[i]) for i in range(row)])
xfmt2 = mdates.DateFormatter('%H:%M') #format de l'eix x
Txfmt2=str(unix2date(t[0]))#format de la finestra



row,col=np.shape(W)
#lines to convert the time array
xti=np.array([unix2date(ti[i]) for i in range(row)])
xfmt = mdates.DateFormatter('%H:%M') #format de l'eix x
Txfmt=str(unix2date(ti[0]))#format de la finestra


####GRAFIC DE STATE DEL NSTRE MeTODE
cmap2=mpl.colors.ListedColormap(['purple','silver','khaki','aqua','black'])
ax1=plt.subplot(1,1,1)
plt.gcf().autofmt_xdate()
plt.pcolor(xti,He,Sm.T,cmap=cmap2,vmin=-20.,vmax=20.)#xti,y,
ax1=plt.gca()
ax1.xaxis.set_major_formatter(xfmt)
plt.title('Hydrometeor type at '+Data[0:10])#INTENSE RAIN DAY

plt.ylabel('Height (m)')
plt.ylim([0,He[-1]])
cbar = plt.colorbar()
cbar.set_ticks([-15., -10., 0.,5.,15])
cbar.set_ticklabels(['Hail', 'Snow', 'Mixed','Rain','Unkown'])
plt.show()

######GRAFIC DE REFLECTIVITAT DEL NSTRE MeTODE 

ax1=plt.subplot(1,1,1)
plt.gcf().autofmt_xdate()
plt.pcolor(xti,He,Zem.T,cmap='jet',vmin=np.min(Zem),vmax=np.max(Zem))#xti,y,
ax1=plt.gca()
ax1.xaxis.set_major_formatter(xfmt)
plt.title('Equivalent Reflectivity (dBZ) at '+Data[0:10])#intense rain day

plt.ylabel('Height (m)')
plt.ylim([0,He[-1]])
plt.colorbar()

plt.show()


######GRAFIC DE VELOCITAT VERTICAL DEL NOSTRE MeTODE 
ax1=plt.subplot(1,1,1)
plt.gcf().autofmt_xdate()
plt.pcolor(xti,He,Wm.T,cmap='jet',vmin=np.nanmin(W),vmax=np.nanmax(W))#xti,y,
ax1=plt.gca()
ax1.xaxis.set_major_formatter(xfmt)
plt.title('Fall speed (m/s) at '+Data[0:10])#intense rain day

plt.ylabel('Height (m)')
plt.ylim([0,He[-1]])
plt.colorbar()

plt.show()
##
