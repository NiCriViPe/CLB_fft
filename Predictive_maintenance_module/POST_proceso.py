##  Time
import time
from pathlib import Path

##  OS
import os
##  Tkinter
import tkinter as tk
##  Pandas
from pandas import read_csv,DataFrame,concat
##  Numpy
from numpy import linspace,size,array,diff
from numpy import log10
##  Scipy
from scipy.fftpack import fft
from scipy.signal import find_peaks,welch
# from scipy.signal.windows import hamming
##  Matplotlib.pyplot
import matplotlib.pyplot as plt
from matplotlib.pyplot import ylim,xlabel,ylabel
from matplotlib.pyplot import grid,minorticks_on,tight_layout
from matplotlib.pyplot import show,vlines,text
from matplotlib.pyplot import legend,title
from matplotlib.pyplot import rcParams



#### Hallazgos para mejora
#Funcion "ejecutar" solo incluye un canal

################################################
# Post proceso,
# Funciones de ejecución de rutinas
################################################
####         Futuras ideas
###     Obtención de un estimador del deslizamiento en base a snom_0
###     Generaciónd de secuencias de armónicas (para posterior etiquetado)
###     Etiquetado automático de armónicas en base a criterios de ubicación
# en función del deslizamiento
###     Definir una función para hacer seguimiento del comportamiento de las armónicas
###     Implementar "frec_criticas" para ranuración, excentricidad, etc...

# class PostProcessing(object):
#     ## Data
#     all = []
#     ## Default paths
#     # Folder Names
    
#     #composing directories
#     def __init__(self,
#                  name:str
#                  ):
#         assert isinstance(name,str),f"Name {name} is a string"
#         self.name=name
#         self.date =None
#         self.model=None
#         self.y_t =None
#         self.FFT =None
#     @staticmethod
def ejecutar(name):
    # En proceso de convertirse en "PREPROCESADO"
    t1=time.time()     
    # [ruta,texto_1]  =   Integracion_de_archivos(v)
    # dat = read_csv(ruta,header=0)

    data = Solo_lectura(name)
    ## Generar serie de tiempo    
    # Ntrunc  =   1
    # data = dat.truncate(after=size(dat,0)//Ntrunc)
    #################
    ##ADVERTENCIA
    ##MODIFICAR PARA INCLUIR MÁS COLUMNAS DE CORRIENTE
    #################
    a       = (data['I1'])  ## Corriente de estator.
    Fs      =   20e3        ## Frecuencia de muestreo original
    Ts      =   1/Fs        ## Tiempo de muestreo\Sampling
    t       = linspace(0,(size(data,0)-1)*Ts,size(data,0))
    y = a
    t = t
    ## Escribir valores instantáneos
    YT2CSV(t,y,name)
    # show(block=True)
    t2=time.time()
    exe_time = t2-t1
    print('tiempo de ejecución: %3.2f segundos' %exe_time)
    print('ESTADO: EJECUTAR OK')
    # tk.messagebox.showinfo(title="EJECUTAR",message="Finalizado" )
    return [t,y]
def YT2CSV(t,y,name):
    d1  ={'time':t,'y':y}
    df1 =DataFrame(data=d1)
    p1=name.replace(".csv","_y-t_DATA.csv")
    df1.to_csv(p1,index=False)
    print("Se guardó forma temporal: ")
    print(p1)
    return    
def PSD(v):
    if  (os.path.isfile(v.replace(".csv","_y-t_DATA.csv"))) is False:
        ejecutar(v)

    yt  =   v.replace(".csv","_y-t_DATA.csv")

    t1  =   time.time()
    dat =   read_csv(yt,header=0)
    x1  =  dat.values[:,0]
    dt  = x1[1] - x1[0]
    fs  = 1 / dt
    #fs  =   20e3
    y1  = dat.values[:,1]
    Ntrunc =    28
    print('Ntrunc:',Ntrunc)
    sy  =   len(y1)
    
    # f, Pxxf = welch(y1,fs=fs,nperseg=sy/Ntrunc,detrend=False,average='median')
    f, Pxxf = welch(y1,fs=fs,nperseg=sy/Ntrunc,detrend=False,average='median')
    label   = v.split("/")[-1]
    Title   = v.split("/")[-2]
    Pxxf    = Pxxf/max(Pxxf)
    fig, (ax1) = plt.subplots(1, 1)
    ax1.plot(f,10*log10(Pxxf),alpha=1,color='#0000CC',label=label)
    # ax1.plot(f,10*log10(Pxxf),alpha=1,color='#0000CC',label=label)
    ax1.set_title(Title)
    grid(True,'both')
    minorticks_on()
    ylim([-120,0])
    ylabel("PSD [dB/Hz]")
    xlabel("Frequency [Hz]")
    print('del_f:',f[2]-f[1],'[Hz]')
    PSD2CSV(f,10*log10(Pxxf),v)
    t2=time.time()
    exe_time= t2-t1
    # ax = gca()    
    # p   =   -50
    p   =   None
    d   =   10*100
    h   =   -90
    peakplot(ax1,f,10*log10(Pxxf),p,d,h)
    ax1.legend()
    show(block=True)

    print('tiempo de ejecución: %3.2f segundos' %(exe_time))
    return

def low_freq_fft(v):
    yt  =   v.replace(".csv","_y-t_DATA.csv")
    start = time.time()
    print(v)
    dat =   read_csv(yt,header=0)
    #### Cálculo de FFT ###
    #### Recordar ajustar la función si cambia el equipo
    #### Equipo: HIOKI PW3198
    ### Ayuda: df  =  f_sample / N
    #dat = read_csv(v,header=0)
    # #### Half-Chop
    # Ntrunc  =   2
    # dat = dat.truncate(after=size(dat,0)//Ntrunc)
    Nfinal = int(4e6)-1
    dat = dat.truncate(after= Nfinal)
    #data = data.truncate(after=size(data,0)//2)
    a       =   dat.values[:,1]  ## Corriente de estator.
    t       =   dat.values[:,0]
    dt      =   t[1]-t[0]
    Fs      =   1/dt
    # Fs      =   20e3        ## Frecuencia de muestreo original
    # Ts      =   1/Fs        ## Tiempo de muestreo\Sampling.
    # t       = linspace(0,size(data,0)*Ts,size(data,0))
    L       = size(t,0)     ## Número de muestras.
    N       = L                   ## Number of samplepoints
    print("Points: ",N)
    ReSample_flag = True
    #### Resample
    if ReSample_flag:
        nskip   =   20
        # nskip   =   16     #Aceptable
        # nskip   = 10    #Recomendado
        auxlist = range(0,N,nskip)
        #
        # a1       =   a[auxlist]
        # y       =   array(a1)
        a       =   a[auxlist]
        a       =   array(a)
        Fs      =   Fs/nskip
        Ts      =   1/Fs
        t       =   linspace(0,(size(a)-1)*Ts,size(a))
        L       = size(t,0)     ## Número de muestras.
        N = L                   ## Number of samplepoints
        #print(size(a),size(),size(a.iloc[auxlist]))
        print("Points: ",N)
        # auxlist1 = range(0,N//2)
    ##Variables finales
    Y       = array(a)
    yf      = fft(Y)
    xf      = linspace(0.0, 1.0/(2.0*Ts), N//2)
    print(size(xf),size(yf))
    yf2     =2.0/N * abs(yf[:N//2])
    aux1    = len(xf)
    #yf2 = 100*yf2/aux1
    yf2     = yf2/aux1
    a       = a/max(a)
    end1    = time.time()
    fft_time=end1-start
    print("FFT elapsed time: %3.2f[s]" %fft_time)
    print("Tamaño x: ",size(yf2),'pts')
    print("Tamaño y: ",size(xf),'pts')
    print("Máxima frecuencia:",(xf[-1]),'[Hz]')
    print("Resolución frecuencia:",(xf[2]-xf[1]),'[Hz]')

    ##############  Create Figure 
    fig, (ax2) = plt.subplots(1, 1)
    ## Plot I vs t
    #ax1.set_xlim(0, 5)
    #ax1.set_xlim(0, t[size(t)-1])
    # ax1.set_xlim(0, max(t))
    # ax1.set_ylim(-1.2, 1.2)
    # ax1.set_xlabel('Tiempo [s]')
    # ax1.set_ylabel('Magnitud relativa [%]')
    # ax1.grid(visible=True, which='minor', color='k', linestyle='--',alpha = 0.5)
    # ax1.minorticks_on()
    # #ax1.plot(t,a)
    # ax1.plot(t,a,'r',alpha = 0.8)
    # ax1.locator_params(axis='y', nbins=10)
    # ax1.locator_params(axis='x', nbins=10)
    ##### PLOT FFT
    ax2.set_xlim(0, max(xf))
    ax2.set_ylim(-70, 0)
    # ax2.set_ylim(1e-6, 100)
    ax2.set_xlabel('Frecuencia [Hz]')   
    # ax2.set_ylabel('Magnitud relativa [%]')
    ax2.set_ylabel('FFT [dB]')
    ax2.grid(visible=True, which='minor',
            color='k', linestyle='--',alpha = 0.5)
    ax2.minorticks_on()
    #ax2.semilogy(xf, yf)
    #ax2.plot(xf,yf)
    # ax2.plot(xf,yf2,'r',alpha=0.95,)
    ax2.plot(xf,10*log10(yf2),'r',alpha=0.95,)
    # ax2.stem(xf,yf,'r')
    ax2.locator_params(axis='y', nbins=10)
    ax2.locator_params(axis='x', nbins=10)
    
    plt.subplots_adjust(
        left=None, bottom=0.1,
        right=None, top=0.8,
        wspace=None, hspace=None)
    tight_layout()
    show(block=False)
    end1 = time.time()
    exe_time=end1-start
    print("plot elapsed time: %3.2f [s]" %exe_time)

    return 

def high_freq_fft(v):
    yt  =   v.replace(".csv","_y-t_DATA.csv")
    start = time.time()
    print(v)
    dat =   read_csv(yt,header=0)
    #### Cálculo de FFT ###
    #### Recordar ajustar la función si cambia el equipo
    #### Equipo: HIOKI PW3198
    ### Ayuda: df  =  f_sample / N
    #dat = read_csv(v,header=0)
    # #### Half-Chop
    # Ntrunc  =   2
    # dat = dat.truncate(after=size(dat,0)//Ntrunc)
    Nfinal = int(4e6)-1
    dat = dat.truncate(after= Nfinal)
    #data = data.truncate(after=size(data,0)//2)
    a       =   dat.values[:,1]  ## Corriente de estator.
    t       =   dat.values[:,0]
    dt      =   t[1]-t[0]
    Fs      =   1/dt
    # Fs      =   20e3        ## Frecuencia de muestreo original
    Ts      =   1/Fs        ## Tiempo de muestreo\Sampling.
    # t       = linspace(0,size(data,0)*Ts,size(data,0))
    L       = size(t,0)     ## Número de muestras.
    N       = L                   ## Number of samplepoints
    print("Points: ",N)
    ReSample_flag = False
    #### Resample
    if ReSample_flag:
        nskip   =   20
        # nskip   =   16     #Aceptable
        # nskip   = 10    #Recomendado
        auxlist = range(0,N,nskip)
        #
        # a1       =   a[auxlist]
        # y       =   array(a1)
        a       =   a[auxlist]
        a       =   array(a)
        Fs      =   Fs/nskip
        Ts      =   1/Fs
        t       =   linspace(0,(size(a)-1)*Ts,size(a))
        L       = size(t,0)     ## Número de muestras.
        N = L                   ## Number of samplepoints
        #print(size(a),size(),size(a.iloc[auxlist]))
        print("Points: ",N)
        # auxlist1 = range(0,N//2)
    ##Variables finales
    Y       = array(a)
    yf      = fft(Y)
    xf      = linspace(0.0, 1.0/(2.0*Ts), N//2)
    print(size(xf),size(yf))
    yf2     =2.0/N * abs(yf[:N//2])
    aux1    = max(yf2)
    #yf2 = 100*yf2/aux1
    yf2     = yf2/aux1
    a       = a/max(a)
    end1    = time.time()
    fft_time=end1-start
    print("FFT elapsed time: %3.2f[s]" %fft_time)
    print("Tamaño x: ",size(yf2),'pts')
    print("Tamaño y: ",size(xf),'pts')
    print("Máxima frecuencia:",(xf[-1]),'[Hz]')
    print("Resolución frecuencia:",(xf[2]-xf[1]),'[Hz]')

    ##############  Create Figure 
    fig, (ax2) = plt.subplots(1, 1)
    ## Plot I vs t
    #ax1.set_xlim(0, 5)
    #ax1.set_xlim(0, t[size(t)-1])
    # ax1.set_xlim(0, max(t))
    # ax1.set_ylim(-1.2, 1.2)
    # ax1.set_xlabel('Tiempo [s]')
    # ax1.set_ylabel('Magnitud relativa [%]')
    # ax1.grid(visible=True, which='minor', color='k', linestyle='--',alpha = 0.5)
    # ax1.minorticks_on()
    # #ax1.plot(t,a)
    # ax1.plot(t,a,'r',alpha = 0.8)locator_params
    # ax1.locator_params(axis='y', nbins=10)
    # ax1.locator_params(axis='x', nbins=10)
    ##### PLOT FFT
    ax2.set_xlim(0, max(xf))
    ax2.set_ylim(-70, 0)
    # ax2.set_ylim(1e-6, 100)
    ax2.set_xlabel('Frecuencia [Hz]')
    # ax2.set_ylabel('Magnitud relativa [%]')
    ax2.set_ylabel('FFT [dB]')
    ax2.grid(visible=True, which='minor', color='k', linestyle='--',alpha = 0.5)
    ax2.minorticks_on()
    #ax2.semilogy(xf, yf)
    #ax2.plot(xf,yf)
    # ax2.plot(xf,yf2,'r',alpha=0.95,)
    ax2.plot(xf,10*log10(yf2),'r',alpha=0.95,)
    # ax2.stem(xf,yf,'r')
    ax2.locator_params(axis='y', nbins=10)
    ax2.locator_params(axis='x', nbins=10)
    
    plt.subplots_adjust(
        left=None, bottom=0.1,
        right=None, top=0.8,
        wspace=None, hspace=None)
    tight_layout()
    show(block=False)
    end1 = time.time()
    exe_time=end1-start
    print("plot elapsed time: %3.2f [s]" %exe_time)

    return

def harm_growth():
    return

def comparar_FFT(n):
    lista_FFT   =   gen_lista_comparar(n)
    compara_plot(lista_FFT)
    # while len(v) != 0:
    ### Función para superponer dos o más espectros
    ### Se debe seleccionar un espectro
    return

def comparar_PSD(n):
    lista_PSD   =   gen_lista_comparar(n)
    compara_plot_PSD(lista_PSD)
    # while len(v) != 0:
    ### Función para superponer dos o más espectros
    ### Se debe seleccionar un espectro
    return






################################################
################################################
# Funciones requeridas para las rutinas
# (Ex "Custom_FUNctions.py)
################################################
################################################

rcParams['font.family'] = 'monospace'   #Estilo de texto


### Módulo de integración de archivos .csv exportados dede PQONE
## La salida "y1" corresponde a la ruta del archivo integrado
## La salida "y2" es un string compuesto
## en caso de ser neceasario monitorear la ejecución
## Opt: Revisar posibilidad de no generar un archivo concatenado "..._MERGED"

def Solo_lectura(s):
    t1  =   time.time()
    # file_list   =   []                        #Lista de archivos
    df_list     =   []                          #Lista de dataframes
    # file_list.append(s)                     #Agregar a la lista de archivos
    df_list.append(read_csv(s))
    print('Se considerarán los siguientes archivos \n') #GUI
    print(s, end='\n')                                  #GUI
    ####        Proceso de búsqueda de archivos #####
    flag = True
    for i in range(1,10):
        if flag:            
            if os.path.isfile(s.replace('.csv','_00'+str(i)+'.csv')):
                aux1 =  s.replace('.csv','_00'+str(i)+'.csv')
                print(aux1)        
                # file_list.append(aux1)
                df_list.append(read_csv(aux1))
            else:
                flag = False
                print('No se encontraron más archivos en el directorio')
    #####       Procesado/ Integración       ########
    a   =   concat(df_list)
    #Integracion = open(s.replace('.csv','_MERGED.csv'),'w')  #Archivo final, abrir
    # print(len(df_list),len(df_list[0]),len(a))
    t2  =   time.time()
    exe_time = t2-t1
    print("lectura en %3.2f [s]" %exe_time)
    return a


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
        #Archivo final, abrir
        Integracion = open(s.replace('.csv','_MERGED.csv'),'w')
        
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



## Funciones para obtención de armónicas
def slip2frec(x):
    return (1-2*x)*50
def slip2frecOsc(x):
    return (1+2*x)*50
def frec2slip(x):
    return (1 -x/50)/2

# def CalcularFFT(v):
#     print(v)
#     start = time.time()
#     #### Cálculo de FFT ###
#     #### Estos datos son fijos
#     #### La frecuencia de muestreo es la del equipo analizador
#     #### Los valores se almacenan en la columna "I1"
#     #### Recordar ajustar la función si cambia el equipo
#     #### Equipo: HIOKI PW3198

#     ### Nota: truncando 1/8, con Nskip = 20 funciona bien
#     ### Nota: truncando 1/1, con Nskip = 39 límite bueno.
#     dat = read_csv(v,header=0)
#     # #### Half-Chop
#     Ntrunc  =   1
#     data = dat.truncate(after=size(dat,0)//Ntrunc)
#     ## 
#     #data = data.truncate(after=size(data,0)//2)
#     a       = (data['I1'])  ## Corriente de estator.
#     Fs      =   20e3        ## Frecuencia de muestreo original
#     Ts      =   1/Fs        ## Tiempo de muestreo\Sampling.
#     t       = linspace(0,size(data,0)*Ts,size(data,0))
#     L       = size(t,0)     ## Número de muestras.
#     N = L                   ## Number of samplepoints
#     print("Points: ",N)

#     ReSample_flag = False
#     #### Resample
#     if ReSample_flag:
#         nskip   =   20
#         # nskip   =   16     #Aceptable
#         # nskip   = 10    #Recomendado
#         # nskip   =   16
#         # nskip   =   30
#         # nskip   =   40     #Demasiado
#         auxlist = range(0,N,nskip)
#         #
#         # a1       =   a[auxlist]
#         # y       =   array(a1)
#         a       =   a[auxlist]
#         # y       =   array(a)
#         Fs      =   Fs/nskip
#         Ts      =   1/Fs
#         t       =   linspace(0,size(a)*Ts,size(a))
#         L       = size(t,0)     ## Número de muestras.
#         N = L                   ## Number of samplepoints
#         #print(size(a),size(),size(a.iloc[auxlist]))
#         print("Points: ",N)
#         # auxlist1 = range(0,N//2)


#     ##Variables finales
#     Y       = array(a)
#     yf      = fft(Y)
#     xf      = linspace(0.0, 1.0/(2.0*Ts), N//2)
#     print(size(xf),size(yf))
#     yf2     =2.0/N * abs(yf[:N//2])
#     aux1    = max(yf2)
#     #yf2 = 100*yf2/aux1
#     yf2     = yf2/aux1
#     a       = a/max(a)
#     end1    = time.time()
#     fft_time=end1-start
#     print("FFT elapsed time: %3.2f[s]" %fft_time)
#     print("Tamaño x: ",size(yf2),'pts')
#     print("Tamaño y: ",size(xf),'pts')
#     print("Resolución frecuencia:",(xf[2]-xf[1]),'[Hz]')
#     return  [xf,100*yf2]

def plotear(xf,yf,t,a):
    start = time.time()
    ##############  Create Figure 
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ## Plot I vs t
    #ax1.set_xlim(0, 5)
    #ax1.set_xlim(0, t[size(t)-1])
    ax1.set_xlim(0, max(t))
    ax1.set_ylim(-1.2, 1.2)
    ax1.set_xlabel('Tiempo [s]')
    ax1.set_ylabel('Magnitud relativa [%]')
    ax1.grid(visible=True, which='minor', color='k', linestyle='--',alpha = 0.5)
    ax1.minorticks_on()
    #ax1.plot(t,a)
    ax1.plot(t,a,'r',alpha = 0.8)
    ax1.locator_params(axis='y', nbins=10)
    ax1.locator_params(axis='x', nbins=10)
    ##### PLOT FFT
    ax2.set_xlim(0, max(xf))
    ax2.set_ylim(1e-6, 100)
    ax2.set_xlabel('Frecuencia [Hz]')
    ax2.set_ylabel('Magnitud relativa [%]')
    ax2.grid(visible=True, which='minor', color='k', linestyle='--',alpha = 0.5)
    ax2.minorticks_on()
    #ax2.semilogy(xf, yf)
    #ax2.plot(xf,yf)
    ax2.plot(xf,yf,'r',alpha=0.95,)
    # ax2.stem(xf,yf,'r')
    ax2.locator_params(axis='y', nbins=10)
    ax2.locator_params(axis='x', nbins=10)
    
    plt.subplots_adjust(
        left=None, bottom=0.1,
        right=None, top=0.8,
        wspace=None, hspace=None)
    tight_layout()
    show(block=False)
    end1 = time.time()
    exe_time=end1-start
    print("plot elapsed time: %3.2f [s]" %exe_time)
    return [ax2,fig]

def peakplot(ax,x,y,p,d,h):
    ### To do
    ## Definir criterio de parámetros de "find_peaks()"
    peaks, _ = find_peaks(y, distance=d,prominence=p,height=h)
    diff(peaks)
    ax.plot(x[peaks], y[peaks], "x",color='green')
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
            text(abs(x),
                0.5 ,
                "%.2f Hz" % frec_crit[i],
                rotation=90, verticalalignment='center')
        else:
            text(abs(x),
                0.5 ,
                "%.2f Hz" % frec_crit[i],
                rotation=90, verticalalignment='center')
            vlines(x=abs(frec_crit[i]),
                ymin = 0, ymax =1 ,
                color = 'purple',alpha = 0.5)
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

def PSD2CSV(f,Px,v):
    d1  ={'frec':f,'dB/Hz':Px}
    df1 =DataFrame(data=d1)
    p1  =v.replace(".csv","_PSD_DATA.csv")
    df1.to_csv(p1,index=False)
    print("Se guardó PSD: ")
    print(p1)
    return


OUTPUT_PATH = Path(__file__).parent


def gen_lista_comparar(n):
    lista_de_archivos = list()
    for i in range(n):
        v = tk.filedialog.askopenfilename(
            initialdir = OUTPUT_PATH,
            title = "Select file",
            filetypes = (("csv files","*.csv"),("all files","*.*")))
        lista_de_archivos.append(v)
    return lista_de_archivos

def compara_plot(lista_de_archivos):
    fig, (ax1) = plt.subplots(1, 1)
    for v in lista_de_archivos:
        df = read_csv(v)
        # aux = v.split("/")[-1].split("_")
        aux = v.split("/")
        auxlabel = aux[-4]+aux[-3]
        ax1.plot(df.values[:,0],df.values[:,1],
                 alpha=0.5,label=auxlabel)

    ax1.set_xlim(0,max(df.values[:,0]))
    ax1.set_ylim(0,100)   
    ax1.grid(visible=True, which='minor',
              color='k', linestyle='--',alpha = 0.25)
    ax1.minorticks_on()
    ax1.set_xlabel('Frecuencia [Hz]')
    ax1.set_ylabel("PSD [dB/Hz]")
    # ax1.set_ylabel('Magnitud relativa [%]')
    ax1.locator_params(axis='y', nbins=10)
    ax1.locator_params(axis='x', nbins=8)
    tight_layout()
    fig.set_tight_layout(True)
    show(block=False)
    legend()
    title(' '.join(aux[0:2]))
    print("Listo")
    return


def compara_plot_PSD(lista_de_archivos):
    fig, (ax1) = plt.subplots(1, 1)
    fig.set_tight_layout(True)
    for v in lista_de_archivos:
        df = read_csv(v)
        # aux = v.split("/")[-1].split("_")
        aux = v.split("/")
        auxlabel = ' '.join(aux[-4:-1])
        ax1.plot(df.values[:,0],df.values[:,1],
                alpha=0.5,label=auxlabel)
                #  alpha=0.5,label=' '.join(aux[0:4]))

    ax1.set_xlim(0,max(df.values[:,0]))
    #ax1.set_ylim([0,100])   
    ax1.grid(visible=True, which='minor',
            color='k', linestyle='--',alpha = 0.25)
    ax1.minorticks_on()
    ax1.set_xlabel('Frecuencia [Hz]')
    #ax1.set_ylabel('Magnitud relativa [%]')
    ax1.set_ylabel("PSD [dB/Hz]")
    # ax1.set_ylabel('Magnitud relativa [%]')
    ax1.locator_params(axis='y', nbins=10)
    ax1.locator_params(axis='x', nbins=8)
    # tight_layout()
    show(block=False)
    legend()
    title(' '.join(aux[0:2]))
    print("Listo")
    return
################### ChangeLog
## Fixed time print n° digits
## Removed _MERGED file creation
