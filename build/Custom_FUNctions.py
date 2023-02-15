# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 13:17:48 2023
@author: vidal
"""
from pandas import read_csv
from pandas import DataFrame
from numpy import linspace,size,array,diff,log
from scipy.fftpack import fft
from scipy.signal import find_peaks
from matplotlib.pyplot import subplots,xlim,ylim,xlabel,ylabel,grid,minorticks_on,tight_layout
from matplotlib.pyplot import show,subplots_adjust,vlines,text,semilogy,subplot2grid,legend,title
from tkinter import filedialog,simpledialog
import time

####         Futuras ideas
###     Generar un .csv con los puntos de la FFT para no ejecutar el cálculo más de una vez
###     Definir un número fijo de puntos para la FFT "independiente" del tamaño de la muestra
###     Obtención de un estimador del deslizamiento en base a snom_0
###     Generaciónd de secuencias de armónicas (para posterior etiquetado)
###     Etiquetado automático de armónicas en base a criterios de ubicación en función del deslizamiento
###     Definir una función para hacer seguimiento del comportamiento de las armónicas
###     Implementar "frec_criticas" para ranuración, excentricidad, etc...
###     Implementar "Persistencia" en las figuras para evaluar en la gráfica los cambios



## Funciones para obtención de armónicas
def slip2frec(x):
    return (1-2*x)*50
def slip2frecOsc(x):
    return (1+2*x)*50
def frec2slip(x):
    return (1 -x/50)/2

def CalcularFFT(v):
    print(v)
    start = time.time()
    #### Cálculo de FFT ###
    #### Estos datos son fijos
    #### La frecuencia de muestreo es la del equipo analizador
    #### Los valores se almacenan en la columna "I1"
    #### Recordar ajustar la función si cambia el equipo
    #### Equipo: HIOKI PW3198
    #data = read_csv(v,header=0)
    dat = read_csv(v,header=0)
    data = dat.truncate(after=size(dat,0))
    #data = data.truncate(after=size(data,0)//2)
    a       = (data['I1'])  ## Corriente de estator.
    Fs      =   20e3        ## Frecuencia de muestreo original
    Ts      =   1/Fs        ## Tiempo de muestreo\Sampling.
    t       = linspace(0,size(data,0)*Ts,size(data,0))
    L       = size(t,0)     ## Número de muestras.
    N = L                   ## Number of samplepoints
    print("Points: ",N)
    # nskip   =   20
    # auxlist = range(0,N,nskip)

    # #### Resample
    # a1       =   a[auxlist]
    # y       =   array(a1)
    # Fs      =   Fs/nskip
    # Ts      =   1/Fs
    # t       =   linspace(0,size(a1)*Ts,size(a1))
    # L       = size(t,0)     ## Número de muestras.
    # N = L                   ## Number of samplepoints
    # #print(size(a),size(),size(a.iloc[auxlist]))
    # print("Points: ",N)
    # auxlist1 = range(0,N//2)

    # #### Half-Chop
    # a2       =   a1[auxlist1]
    # y1       =   array(a2)
    # #Fs      =   Fs/nskip
    # #Ts      =   1/Fs
    # t       =   linspace(0,size(a2)*Ts,size(a2))
    # L       =   size(t,0)     ## Número de muestras.
    # N = L                   ## Number of samplepoints
    # #print(size(a),size(),size(a.iloc[auxlist]))
    # print(N)    

    ##Variables finales
    Y   = array(a)
    yf = fft(Y)
    xf = linspace(0.0, 1.0/(2.0*Ts), N//2)
    print(size(xf),size(yf))
    yf2 =2.0/N * abs(yf[:N//2])
    aux1 = max(yf2)
    #yf2 = 100*yf2/aux1
    yf2 = yf2/aux1
    a   = a/max(a)
    end1 = time.time()
    print("FFT elapsed time:",end1-start)
    print("Tamaño x: ",size(yf2))
    print("Tamaño y: ",size(xf))
    return  [xf,100*yf2,t,a/max(a)]

def plotear(xf,yf,t,a):
    start = time.time()
    ##############  Create Figure 
    fig, (ax1, ax2) = subplots(1, 2)
    ## Plot I vs t
    #ax1.set_xlim(0, 5)
    ax1.set_xlim(0, 280)
    ax1.set_ylim(-1.2, 1.2)
    ax1.set_xlabel('Tiempo [s]')
    ax1.set_ylabel('Magnitud relativa [%]')
    ax1.grid(visible=True, which='minor', color='k', linestyle='--',alpha = 0.5)
    ax1.minorticks_on()
    ax1.plot(t,a)

    ##### PLOT FFT
    #ax2.set_xlim(0, 500)
    ax2.set_xlim(0, max(xf))
    ax2.set_ylim(1e-6, 100)
    ax2.set_xlabel('Frecuencia [Hz]')
    ax2.set_ylabel('Magnitud relativa [%]')
    ax2.grid(visible=True, which='minor', color='k', linestyle='--',alpha = 0.5)
    ax2.minorticks_on()
    #ax2.semilogy(xf, yf)
    ax2.plot(xf,yf)
    
    subplots_adjust(left=None, bottom=0.1, right=None, top=0.8, wspace=None, hspace=None)
    tight_layout()
    show(block=False)
    end1 = time.time()
    print("plot elapsed time:",end1-start)
    return [ax2,fig]

def peakplot(ax,x,y,p,d,h):
    ### To do
    ## Definir criterio de parámetros de "find_peaks()"
    peaks, _ = find_peaks(y, distance=d,prominence=p,height=h)
    diff(peaks)
    ax.plot(x[peaks], y[peaks], "x")
    return peaks

def brb_harm(snom,n):
    frec_crit=list()
    for i in range(1,n+1):
        frec_crit.append(slip2frec(snom*i))
        frec_crit.append(slip2frecOsc(snom*i))
    frec_crit=array(frec_crit)

    vlines(x=abs(frec_crit[0]), ymin = 0, ymax =1 , color = 'red',alpha = 0.5)
    for i, x in enumerate(frec_crit):
        if i == 0:
            text(abs(x), 0.5 , "%.2f Hz" % frec_crit[i], rotation=90, verticalalignment='center')
        else:
            text(abs(x), 0.5 , "%.2f Hz" % frec_crit[i], rotation=90, verticalalignment='center')
            vlines(x=abs(frec_crit[i]), ymin = 0, ymax =1 , color = 'purple',alpha = 0.5)
    return

## Armónicas de rodamientos
def rod_harm(arg):
    return


##  Escritura csv de la FFT
def FFT2CSV(xf,yf,v):
    d1  ={'frec':xf,'M':yf}
    df1 =DataFrame(data=d1)
    p1  =v.replace(".csv","_FFT_DATA.csv")
    df1.to_csv(p1,index=False)
    print("Se guardó FFT: ")
    print(p1)
    return
def YT2CSV(t,y,v):
    d1  ={'time':t,'y':y}
    df1 =DataFrame(data=d1)
    p1=v.replace(".csv","_y-t_DATA.csv")
    df1.to_csv(p1,index=False)
    print("Se guardó forma temporal: ")
    print(p1)
    return


def gen_lista_comparar(n):
    lista_de_archivos = list()
    for i in range(n):
        v = filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("csv files","*.csv"),("all files","*.*")))
        lista_de_archivos.append(v)
    return lista_de_archivos

def compara_plot(lista_de_archivos):
    fig, (ax1) = subplots(1, 1)
    for v in lista_de_archivos:
        df = read_csv(v)
        aux = v.split("/")[-1].split("_")
        ax1.plot(df.values[:,0],df.values[:,1],alpha=0.5,label=' '.join(aux[0:4]))

    ax1.set_xlim([0,max(df.values[:,0])])
    ax1.set_ylim([0,100])
    ax1.grid(visible=True, which='minor', color='k', linestyle='--',alpha = 0.25)
    ax1.minorticks_on()
    ax1.set_xlabel('Frecuencia [Hz]')
    ax1.set_ylabel('Magnitud relativa [%]')
    tight_layout()
    show(block=False)
    legend()
    title(' '.join(aux[0:2]))
    print("ploteado")
    return

### Módulo de integración de archivos .csv exportados dede PQONE
## La salida "y1" corresponde a la ruta del archivo integrado
## La salida "y2" es un string compuesto en caso de ser neceasario monitorear la ejecución
## Opt: Revisar posibilidad de no generar un archivo concatenado "..._MERGED"
import os
def Integracion_de_archivos(s):
    if ("MERGE" in s):
        y2 = "El archivo seleccionado ya está integrado"
        y1 = s
        print(y2)
        return [y1,y2]
    elif  (os.path.isfile(s.replace('.csv','_MERGED.csv'))):
        y2 = "El archivo seleccionado ya está integrado"
        y1 = s.replace('.csv','_MERGED.csv')
        print(y2)
        return [y1,y2]
    else:
        file_list = []                          #Lista de archivos
        file_list.append(s)                     #Agregar a la lista de archivos
        
        print('Se considerarán los siguientes archivos \n') #GUI
        print(s, end='\n')                                  #GUI
        ####        Proceso de búsqueda de archivos #####
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
        #####       Procesado/ Integración       ########
        Integracion = open(s.replace('.csv','_MERGED.csv'),'w')  #Archivo final, abrir
        
        for i in range(len(file_list)):
            file = open(file_list[i],'r')
            counter = 0
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
        y2 = "La integración ha sido exitosa" 
        y1= s.replace('.csv','_MERGED.csv')
        print(y2)
        return [y1,y2]
    return 0

