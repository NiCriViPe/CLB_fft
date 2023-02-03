# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 15:34:42 2023
@author: vidal
"""
import os
## Prorgama de integración de mediciones
########### Inicialización  ########
file_list = []                          #Lista de archivos
s = input('Nombre de archivo: \n');     #Ruta de primer archivo
s = s.replace("'",'')

file_list.append(s)                     #Agregar a la lista de archivos

print('Se considerarán los siguientes archivos \n') #GUI
print(s, end='\n')                                  #GUI

flag = True
for i in range(1,10):
    if flag:            
        if os.path.isfile(s.replace('.csv','_00'+str(i)+'.csv')):
            aux1 =  s.replace('.csv','_00'+str(i)+'.csv')
            print(aux1)        
            file_list.append(aux1)
        else:
            flag = False
            print('No se encontraron más arhivos en el directorio')


#####       Procesado       ########

Integracion = open(s.replace('.csv','_MERGED.csv'),'w')  #Archivo final, abrir

for i in range(len(file_list)):
    file = open(file_list[i],'r')
    counter = 0
    #print(i)
    for line in file:
        if i == 0:
            Integracion.write(line)
        if i>0 and counter >0:
            Integracion.write(line)
        #if i>0 and counter ==0:
            # print('Se ha omitido la siguiente línea:')
            # print(line)
        counter+=1
    file.close()
Integracion.close()
################## FIN      ##########
print('La integración de las mediciones ha finalizado.')
    
# Nombre_archivo_salida =input('Ingrese el nombre del archivo de salida: ')
print('Comenzará la obtención de la FFT, por favor espere.  ')
print('El proceso tardará menos de dos minutos.\n')
Nombre_archivo_salida   = s.replace('.csv','.pdf')
print('El archivo de salida es:',Nombre_archivo_salida)

from pandas import read_csv
from numpy import linspace,size,array
#import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, subplots,xlim,ylim,grid,minorticks_on
from matplotlib.pyplot import savefig,show,plot
from scipy.fftpack import fft

data = read_csv(s.replace('.csv','_MERGED.csv'),header=0)
########## Script ##########
# ### Cálculo de FFT ###
Fs      =   20e3        ## Frecuencia de muestreo.
Ts      =   1/Fs        ## Tiempo de muestreo.
a       = (data['I1'])  ## Corriente de estator.
t       = linspace(0,size(data,0)*Ts,size(data,0))
L       = size(t,0)     ## Número de muestras.
Ts      =(t[2]-t[1])    ## ALternativa de tiempo de muestreo.
Fs      =1/Ts           ## Alt. Frecuencia de muestreo.
Fn      =Fs/2           ## Frecuencia de Nyquist
LastFreq =Fn            ## última frecuecnia FFT

N = L                   ## Number of samplepoints
T = Ts                  ## sample spacing
x = t
y = array(a)
yf = fft(y)
xf = linspace(0.0, 1.0/(2.0*T), N//2)
### Plot ###
#fig, ax = plt.figure(dpi=1200)
fig, ax = subplots()
# ax.plot(xf, 2.0/N * np.abs(yf[:N//2]))
ax.plot(xf, 2.0/N * abs(yf[:N//2]))
#ax.stem(xf, 2.0/N * abs(yf[:N//2]))
# plt.xlim(0, 100)
xlim(20, 80)
ylim(0, 0.005)
grid(True)
# plt.grid(b=True, which='major', color='b', linestyle='-')
grid(visible=True, which='minor', color='k', linestyle='--',alpha = 0.5)
minorticks_on()
# plt.savefig('filename.png',dpi=1200)
#savefig(Nombre_archivo_salida+'.pdf')
savefig(Nombre_archivo_salida)
show()

print('El proceso ha termindo exitosamente, ya puede cerrar esta ventana')