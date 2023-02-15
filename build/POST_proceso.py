#from Custom_FUNctions import peakplot,brb_harm,Integracion_de_archivos,CalcularFFT,plotear
from Custom_FUNctions import *
from matplotlib.pyplot import show,plot
import tkinter as tk
from tkinter import filedialog,simpledialog

def ejecutar(v):
    [ruta,texto_1]  =   Integracion_de_archivos(v)
    xf,yf,t,a       =   CalcularFFT(ruta)
    FFT2CSV(xf,yf,v)
    YT2CSV(t,a,v)

    [ax2,fig]       =   plotear(xf,yf,t,a)
    h = 0.01            ## Altura mínima

    d = 100             ## Distancia horizontal para peaks
    p = None            ## Parámetro de "Prominencia"
    #peaks = peakplot(ax,xf,yf,p,d,h)

    snom    =  3591/100306
    n       =  10
    #brb_harm(snom,n)
    show(block=True)
    return

def comparar(n):
    lista_FFT   =   gen_lista_comparar(n)
    compara_plot(lista_FFT)
    # while len(v) != 0:

    ### Función para superponer dos o más espectros
    ### Se debe seleccionar un espectro
    return
