import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter.filedialog import askopenfile
import tkinter.scrolledtext as scrolledtext
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyiast
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import math as math
from PIL import ImageTk,Image 

#initialize general variables
Models = ['Interpolator', 'Langmuir', 'Quadratic', 'Henry', 'BET', 'TemkinApprox', 'DSLangmuir']
colours = ['r', 'b', 'y', 'm', 'g','tab:brown', 'c'] #colours for the corresponding model when fitting
selectedModel = ""
selectedColour = ""
Fit1 = []
Fit2 = []
pressures = []
mole_fractions = []
stored_selectivities = [] #selectivities stored for selectivity vs pressure in  plot 3
axisScale = 'linear' #starting axis scale for the plots

#initialize the headers of the loaded data for the fitting graphs
pressure_key1 = ''
loading_key1 = ''
pressure_key2 = ''
loading_key2 = ''

#initialize the custom fit starting parameter guesses for isotherm fit 1
customLangmuir1M = ''
customLangmuir1K = ''
customHenry1KH = ''
customQuadratic1M=''
customQuadratic1Ka=''
customQuadratic1Kb=''
customBET1M=''
customBET1Ka=''
customBET1Kb=''
customDSLangmuir1M1=''
customDSLangmuir1M2=''
customDSLangmuir1K1=''
customDSLangmuir1K2=''
customTemkin1M=''
customTemkin1K=''
customTemkin1Theta=''
initialParameterGuess1 = {} #dictionary that stores the initial fitting parameters for the different models
customParameterarray1 =['','','','','','','','','','','','','','','','']

#initialize the custom fit starting parameter guesses for isotherm fit 2
customLangmuir2M = ''
customLangmuir2K = ''
customHenry2KH = ''
customQuadratic2M=''
customQuadratic2Ka=''
customQuadratic2Kb=''
customBET2M=''
customBET2Ka=''
customBET2Kb=''
customDSLangmuir2M1=''
customDSLangmuir2M2=''
customDSLangmuir2K1=''
customDSLangmuir2K2=''
customTemkin2M=''
customTemkin2K=''
customTemkin2Theta=''
initialParameterGuess2 = {} #dictionary that stores the initial fitting parameters for the different models
customParameterarray2 =['','','','','','','','','','','','','','','','']

#initialize the memory of the last model to be selected.
lastSelectedModel1 = ''
lastSelectedModel2 = ''
lastSelectedColour1 = ''
lastSelectedColour2 = ''

#______________________________________________________________________________
#functions

def import_file1(): #opens a window where the user can choose a csv file with data to load for pure isotherm 1
    global df_gas1 
    filename1 = filedialog.askopenfilename(title = "Select file",filetypes = (("CSV Files","*.csv"),))
    if filename1 != '':
        df_gas1 = pd.read_csv(filename1,  sep='[:;,|_\s+\t]', engine='python')
        return df_gas1
    else:
        return

def import_file2(): #opens a window where the user can choose a csv file with data to load for pure isotherm 2
    global df_gas2 
    filename2 = filedialog.askopenfilename(title = "Select file",filetypes = (("CSV Files","*.csv"),))
    if filename2 != '':
        df_gas2 = pd.read_csv(filename2,  sep='[:;,|_\s+\t]', engine='python')
        return df_gas2
    else:
        return

def entryfield_to_array(entryfield): #takes the input of an entryfield and converts it into an array filled with floats
    string = str(entryfield.get())
    parsedString = string.split(',')
    array = []
    for elem in parsedString:
        try:
            float(elem)
        except:
            error = "Wrongful input in entryfield"
            popup_errormsg(error)
            exit
        number = float(elem)
        array.append(number)
    return array

def setTextInput(textbox, text): #can add a textmessage to the end of a (scrollable) textbox
    textbox.configure(state='normal')
    textbox.insert("end", text)
    textbox.insert("end", '\r\n')
    textbox.insert("end", '\r\n')
    textbox.configure(state='disabled')
    textbox.yview('end')     
     
def fit_isotherm1(): #this function checks which fitting model was selected and tries to fit the loaded data for isotherm 1 with this model. if custom fit parameters have been entered then those will be taken as starting guesses instead of the standard fitting guesses by the model. if successful, the fit will be plotted in the assigned colour for the model in plot 1.
    if len(plot1.axes.lines) > 1: #if a fit is present, it is deleted so it can be replaced with the new fit
        plot1.axes.lines[1].remove()
    
    global initialParameterGuess1
    global axisScale
    global Fit1
    global var
    global Models
    global selectedModel
    global selectedColour
    global lastSelectedModel1
    global pressure_key1
    global loading_key1
    global lastSelectedColour1
    
    position = var.get() #get the selected radiobutton value to determine which model was chosen
    selectedModel = Models[position]
    selectedColour = colours[position]
    try:
        df_gas1
    except:
        error = "No data loaded for Gas 1"
        popup_errormsg(error)
        return
    if selectedModel == "Interpolator":
        Fit1 = pyiast.InterpolatorIsotherm(
                df_gas1,
                loading_key=loading_key1,
                pressure_key=pressure_key1,
                fill_value=df_gas1[loading_key1].max())
        df_gas1array = np.array(df_gas1[loading_key1])
        mse = np.square(np.subtract(df_gas1array, Fit1.loading(df_gas1[pressure_key1]))).mean()
        rmse = math.sqrt(mse) #calculate a root mean squared error for the interpolated data
        firstline = selectedModel + '\r\n'
        secondline = 'RMSE = ' + str(rmse)
        text = firstline + secondline
        setTextInput(textParameters1, text)
    else:
        initialParameterGuess1.clear()
        
        #it will now check for the selected model if custom fit parameters have been entered. if incomplete custom fit parameters have been entered it will generate an error; otherwise it will use the custom fit parameters as a starting gues. if no custom fit parameters have been filled in for the selected model it will use predefined standard fit parameters.
        #the custom fit parameters are loaded from the custom parameter array to the initial parameter guess dictionary so they can be used for fitting
        
        if selectedModel == 'Langmuir' and (customParameterarray1[0] != '' or customParameterarray1[1] != ''):
            if (customParameterarray1[0] != '' and customParameterarray1[1] != ''):
                initialParameterGuess1 = {'M' : float(customParameterarray1[0]), 'K' :float(customParameterarray1[1])}
            else:
                error = "Custom Fit 1: Not all custom guess parameters have been defined for the Langmuir model"
                popup_errormsg(error)
                return
        if selectedModel == 'DSLangmuir' and (customParameterarray1[9] != '' or customParameterarray1[10] != '' or customParameterarray1[11] != '' or customParameterarray1[12] != ''):
           if (customParameterarray1[9] != '' and customParameterarray1[10] != '' and customParameterarray1[11] != '' and customParameterarray1[12] != ''):
                initialParameterGuess1 = {'M1' : float(customParameterarray1[9]), 'M2' :float(customParameterarray1[10]), 'K1' :float(customParameterarray1[11]), 'K2' :float(customParameterarray1[12])}
           else:
                error = "Custom Fit 1: Not all custom guess parameters have been defined for the DSLangmuir model"
                popup_errormsg(error)
                return
        if selectedModel == 'Henry' and (customParameterarray1[2] != ''):
            if (customParameterarray1[2] != ''):
                initialParameterGuess1 = {'KH' : float(customParameterarray1[2])}
            else:
                error = "Custom Fit 1: Not all custom guess parameters have been defined for the Henry model"
                popup_errormsg(error)
                return
        if selectedModel == 'Quadratic' and (customParameterarray1[3] != '' or customParameterarray1[4] != '' or customParameterarray1[5] != ''):
            if (customParameterarray1[3] != '' and customParameterarray1[4] != '' and customParameterarray1[5] != ''):
                initialParameterGuess1 = {'M' : float(customParameterarray1[3]), 'Ka' :float(customParameterarray1[4]), 'Kb' :float(customParameterarray1[5])}
            else:
                error = "Custom Fit 1: Not all custom guess parameters have been defined for the Quadratic model"
                popup_errormsg(error)
                return
        if selectedModel == 'BET' and (customParameterarray1[6] != '' or customParameterarray1[7] != '' or customParameterarray1[8] != ''):
            if (customParameterarray1[6] != '' and customParameterarray1[7] != '' and customParameterarray1[8] != ''):
                initialParameterGuess1 = {'M' : float(customParameterarray1[6]), 'Ka' :float(customParameterarray1[7]), 'Kb' :float(customParameterarray1[8])}
            else:
                error = "Custom Fit 1: Not all custom guess parameters have been defined for the BET model"
                popup_errormsg(error)
                return
        if selectedModel == 'TemkinApprox' and (customParameterarray1[13] != '' or customParameterarray1[14] != '' or customParameterarray1[15] != ''):
            if (customParameterarray1[13] != '' and customParameterarray1[14] != '' and customParameterarray1[15] != ''):
                initialParameterGuess1 = {'M' : float(customParameterarray1[13]), 'K' :float(customParameterarray1[14]), 'theta' :float(customParameterarray1[15])}
            else:
                error = "Custom Fit 1: Not all custom guess parameters have been defined for the TemkinApprox model"
                popup_errormsg(error)
                return
        
        try:
            Fit1 = pyiast.ModelIsotherm(
                df_gas1,
                loading_key=loading_key1,
                pressure_key=pressure_key1,
                model=selectedModel,
                param_guess = initialParameterGuess1)
        except:
            error = "Gas 1: Minimization of Residual Sum of Squares (RSS) isotherm fitting failed." + "\r\n" + "Try different starting parameters in the custom fit menu or select a different model."
            popup_errormsg(error)
            return
        else: #fitting the pure isotherm with the selected model (& initial parameter guesses if entered in the custom fit menu)
            Fit1 = pyiast.ModelIsotherm(
                    df_gas1,
                    loading_key=loading_key1,
                    pressure_key=pressure_key1,
                    model=selectedModel,
                    param_guess = initialParameterGuess1)
                        
            firstline = Fit1.model + '\r\n'
            secondline = ''
            for param, val in Fit1.params.items():
                secondline += "\t%s = %f" % (param, val) + '\r\n'
            thirdline = 'RMSE = ' + str(Fit1.rmse)
            text = firstline + secondline + thirdline
            setTextInput(textParameters1, text)
    plot1.plot(df_gas1[pressure_key1],
               Fit1.loading(df_gas1[pressure_key1],),
               color= selectedColour,
               linewidth=1,
               label='Fit ' + selectedModel)
    plot1.set_xscale(axisScale)
    plot1.legend(loc='best')
    canvasFigure1.draw()
    
    lastSelectedModel1 = selectedModel
    lastSelectedColour1 = selectedColour
    
    return Fit1, lastSelectedModel1
    
def fit_isotherm2(): #this function checks which fitting model was selected and tries to fit the loaded data for isotherm 2 with this model. if custom fit parameters have been entered then those will be taken as starting guesses instead of the standard fitting guesses by the model. if successful, the fit will be plotted in the assigned colour for the model in plot 2.
    if len(plot2.axes.lines) > 1: #if a fit is present, it is deleted so it can be replaced with the new fit
        plot2.axes.lines[1].remove()
    
    global initialParameterGuess2
    global Fit2
    global var
    global Models
    global selectedModel
    global selectedColour
    global lastSelectedModel2
    global lastSelectedColour2
    
    position = var.get() #get the selected radiobutton value to determine which model was chosen
    selectedModel = Models[position]
    selectedColour = colours[position]
    try:
        df_gas2
    except:
        error = "No data loaded for Gas 2"
        popup_errormsg(error)
        return
    if selectedModel == "Interpolator":
        Fit2 = pyiast.InterpolatorIsotherm(
                df_gas2,
                loading_key=loading_key2,
                pressure_key=pressure_key2,
                fill_value=df_gas2[loading_key2].max())
        df_gas2array = np.array(df_gas2[loading_key2])
        mse = np.square(np.subtract(df_gas2array, Fit2.loading(df_gas2[pressure_key2]))).mean()
        rmse = math.sqrt(mse) #calculate a root mean squared error for the interpolated data
        firstline = selectedModel + '\r\n'
        secondline = 'RMSE = ' + str(rmse)
        text = firstline + secondline
        setTextInput(textParameters2, text)
    else:
        initialParameterGuess2.clear()
        
        #it will now check for the selected model if custom fit parameters have been entered. if incomplete custom fit parameters have been entered it will generate an error; otherwise it will use the custom fit parameters as a starting gues. if no custom fit parameters have been filled in for the selected model it will use predefined standard fit parameters.
        #the custom fit parameters are loaded from the custom parameter array to the initial parameter guess dictionary so they can be used for fitting
        
        if selectedModel == 'Langmuir' and (customParameterarray2[0] != '' or customParameterarray2[1] != ''):
            if (customParameterarray2[0] != '' and customParameterarray2[1] != ''):
                initialParameterGuess2 = {'M' : float(customParameterarray2[0]), 'K' :float(customParameterarray2[1])}
            else:
                error = "Custom Fit 2: Not all custom guess parameters have been defined for the Langmuir model"
                popup_errormsg(error)
                return
        if selectedModel == 'DSLangmuir' and (customParameterarray2[9] != '' or customParameterarray2[10] != '' or customParameterarray2[11] != '' or customParameterarray2[12] != ''):
           if (customParameterarray2[9] != '' and customParameterarray2[10] != '' and customParameterarray2[11] != '' and customParameterarray2[12] != ''):
                initialParameterGuess2 = {'M1' : float(customParameterarray2[9]), 'M2' :float(customParameterarray2[10]), 'K1' :float(customParameterarray2[11]), 'K2' :float(customParameterarray2[12])}
           else:
                error = "Custom Fit 2: Not all custom guess parameters have been defined for the DSLangmuir model"
                popup_errormsg(error)
                return
        if selectedModel == 'Henry' and (customParameterarray2[2] != ''):
            if (customParameterarray2[2] != ''):
                initialParameterGuess2 = {'KH' : float(customParameterarray2[2])}
            else:
                error = "Custom Fit 2: Not all custom guess parameters have been defined for the Henry model"
                popup_errormsg(error)
                return
        if selectedModel == 'Quadratic' and (customParameterarray2[3] != '' or customParameterarray2[4] != '' or customParameterarray2[5] != ''):
            if (customParameterarray2[3] != '' and customParameterarray2[4] != '' and customParameterarray2[5] != ''):
                initialParameterGuess2 = {'M' : float(customParameterarray2[3]), 'Ka' :float(customParameterarray2[4]), 'Kb' :float(customParameterarray2[5])}
            else:
                error = "Custom Fit 2: Not all custom guess parameters have been defined for the Quadratic model"
                popup_errormsg(error)
                return
        if selectedModel == 'BET' and (customParameterarray2[6] != '' or customParameterarray2[7] != '' or customParameterarray2[8] != ''):
            if (customParameterarray2[6] != '' and customParameterarray2[7] != '' and customParameterarray2[8] != ''):
                initialParameterGuess2 = {'M' : float(customParameterarray2[6]), 'Ka' :float(customParameterarray2[7]), 'Kb' :float(customParameterarray2[8])}
            else:
                error = "Custom Fit 2: Not all custom guess parameters have been defined for the BET model"
                popup_errormsg(error)
                return
        if selectedModel == 'TemkinApprox' and (customParameterarray2[13] != '' or customParameterarray2[14] != '' or customParameterarray2[15] != ''):
            if (customParameterarray2[13] != '' and customParameterarray2[14] != '' and customParameterarray2[15] != ''):
                initialParameterGuess2 = {'M' : float(customParameterarray2[13]), 'K' :float(customParameterarray2[14]), 'theta' :float(customParameterarray2[15])}
            else:
                error = "Custom Fit 2: Not all custom guess parameters have been defined for the TemkinApprox model"
                popup_errormsg(error)
                return
        
        try:
            Fit2 = pyiast.ModelIsotherm(
                df_gas2,
                loading_key=loading_key2,
                pressure_key=pressure_key2,
                model=selectedModel,
                param_guess = initialParameterGuess2)
        except:
            error = "Gas 2: Minimization of Residual Sum of Squares (RSS) isotherm fitting failed." + "\r\n" + "Try different starting parameters in the custom fit menu or select a different model."
            popup_errormsg(error)
            return        
        else: #fitting the pure isotherm with the selected model (& initial parameter guesses if entered in the custom fit menu)
            Fit2 = pyiast.ModelIsotherm(
                    df_gas2,
                    loading_key=loading_key2,
                    pressure_key=pressure_key2,
                    model=selectedModel,
                    param_guess = initialParameterGuess2)
            firstline = Fit2.model + '\r\n'
            secondline = ''
            for param, val in Fit2.params.items():
                secondline += "\t%s = %f" % (param, val) + '\r\n'
            thirdline = 'RMSE = ' + str(Fit2.rmse)
            text = firstline + secondline + thirdline
            setTextInput(textParameters2, text)
    plot2.plot(df_gas2[pressure_key2],
               Fit2.loading(df_gas2[pressure_key2],),
               color= selectedColour,
               linewidth=1,
               label='Fit ' + selectedModel)
    plot2.legend(loc='best')
    plot2.set_xscale(axisScale)
    canvasFigure2.draw()
    
    lastSelectedModel2 = selectedModel 
    lastSelectedColour2 = selectedColour
    
    return Fit2, lastSelectedModel2

def run_IAST(): #runs the IAST calculation on the fits of isotherm 1 & 2 for the selected pressures & mole fractions. at the end it stores the selectivities calculated from the first mole fraction entered to be plotted in plot 3.
    global Fit1
    global Fit2
    global selectivities
    global entryfieldPressures
    global entryfieldMoleFractions
    global pressures
    global mole_fractions
    global stored_selectivities
    
    if Fit1 == []:
        error = "Gas 1 was not fitted"
        popup_errormsg(error)
        return
    if Fit2 == []:
        error = "Gas 2 was not fitted"
        popup_errormsg(error)
        return
        
    stored_selectivities = []
    total_selectivities = []
    
    try:
        entryfield_to_array(entryfieldPressures)
    except:
        error = "No pressures entered"
        popup_errormsg(error)
        return
    
    try:
        entryfield_to_array(entryfieldMoleFractions)
    except:
        error = "No mole fractions entered"
        popup_errormsg(error)
        return
    
    pressures = entryfield_to_array(entryfieldPressures) #get the entered pressures
    mole_fractions = entryfield_to_array(entryfieldMoleFractions) #get the entered mole fractions
    
    #check if entered mole fractions are allowed
    for i in range(0, len(mole_fractions)): 
        if mole_fractions[i] >= 1.0:
            error = "Mole fractions cannot be higher or equal to 1"
            popup_errormsg(error)
            return
        elif mole_fractions[i] <= 0.0:
            error = "Mole fractions cannot be negative or 0"
            popup_errormsg(error)
            return
        i += 1
    
    #check the maximum pressure of both loaded isotherms to determine how high the entered pressure is allowed to be for the IAST calculations    
    column1 = df_gas1[pressure_key2]
    max_pressure1 = column1.max()
    column2 = df_gas2[pressure_key2]
    max_pressure2 = column2.max()
    if max_pressure1 > max_pressure2:
        max_pressure_final = max_pressure1
    else:
        max_pressure_final = max_pressure2
    
    for i in range(0, len(pressures)):
        if pressures[i] >= max_pressure_final:
            error = "Presssures cannot exceed the highest maximum pressure of the loaded pure adsorption isotherms \r\n Current maximum pressure is: " + str(max_pressure_final) + " bar"
            popup_errormsg(error)
            return
        elif pressures[i] <= 0.0:
            error = "Pressures cannot be negative or 0"
            popup_errormsg(error)
            return
        i += 1
    
    #calculate partial pressures
    n_total_pressures = len(pressures)
    n_mole_fractions = len(mole_fractions)
    n_partial_pressure = n_total_pressures * n_mole_fractions
    partial_pressure = np.eye(n_partial_pressure,2)
    k = 0
    for i in range(n_total_pressures):
        for j in range(n_mole_fractions):
            partial_pressuresA = np.array(pressures[i] * mole_fractions[j]) #partial pressures of gas 1
            partial_pressuresB = np.array(pressures[i] * (1 - mole_fractions[j])) #partial pressures of gas 2
            partial_pressure[k] = ([partial_pressuresA, partial_pressuresB])
            if k <= n_partial_pressure:
                k = k + 1

    # use calculated partial pressures for IAST calculation and selectivity
    selectivities = np.eye(n_partial_pressure, 8)
    m = 0
    i = 0
    j = 0
    for l in range(n_partial_pressure):
        loadings = pyiast.iast(partial_pressure[l], [Fit1, Fit2], verboseflag=False) #calculate loadings for the partial pressures
        selectivity = print_selectivity2(loadings, partial_pressure[l]) #calculate the selectivity (print_selectivity2 was adapted from the pyIAST package, now it gives the selectivities as an output instead of printing to the console)
        total_selectivities.append(selectivity[2]) #store the calculated selectivity
        selectivities[m] = (pressures[i], mole_fractions[j], (1 - mole_fractions[j]), partial_pressure[l,0], partial_pressure[l,1], loadings[0], loadings[1], selectivity[2]) #prepares calculated data for output (pressure, mole fraction gas 1, mole fraction gas 2, partial pressure gas 1, partial pressure gas 2, selectivity)
        if m <= n_partial_pressure and j <= n_mole_fractions:
            m = m + 1
            l = l + 1
            j = j + 1
        if m <= n_partial_pressure and j >=n_mole_fractions:
            j = 0
            i = i + 1
    
    selec = total_selectivities[::len(mole_fractions)]
    for i in range(0, len(selec)):
        stored_selectivities.append(selec[i])
    
    if stored_selectivities == []:
        error = "Something went wrong with the IAST calculation. \r\n Check if the data was correctly imported and fitted (with physically meaningful parameters)"
        popup_errormsg(error)
        return

def print_selectivity2(component_loadings, partial_pressures): #slightly modified function from the PyIAST package that calculates the selectivity for the loadings and pressures input in the entryfields, now it gives the selectivities as an output instead of printing to the console
    """
    Calculate selectivity as a function of component loadings and bulk gas
    pressures

    :param component_loadings: numpy array of component loadings
    :param partial_pressures: partial pressures of components
    """
    n_components = np.size(component_loadings)
    for i in range(n_components):
        for j in range(i + 1, n_components):
            s = (i, j, component_loadings[i] / component_loadings[j] /
                   (partial_pressures[i] / partial_pressures[j]))
    return s

def plot_data1(): #takes the headers of the data for isotherm 1 and plots the isotherm in plot 1
     global pressure_key1
     global loading_key1
     try:
         df_gas1
     except:
         return
     else:
         pressure_key1, loading_key1 = [*df_gas1]
         plot1.axes.clear()
         plot1.plot(df_gas1[pressure_key1],
                 df_gas1[loading_key1],
                 'bo',
                 markersize=3,
                 color='k',
                 label='Gas 1')
         plot1.set_xlabel(pressure_key1, fontsize='x-small')
         plot1.set_ylabel(loading_key1, fontsize='x-small')
         plot1.legend(loc='best')
         plot1.set_xscale(axisScale)
         canvasFigure1.draw()
         plot1.autoscale(enable=False, axis='both', tight=None)
     
def plot_data2(): #takes the headers of the data for isotherm 2 and plots the isotherm in plot 2
     global pressure_key2
     global loading_key2
     try:
         df_gas2
     except:
         return
     else:
         pressure_key2, loading_key2 = [*df_gas2]
         plot2.axes.clear()
         plot2.plot(df_gas2[pressure_key2],
                 df_gas2[loading_key2],
                 'bo',
                 markersize=3,
                 color='k',
                 label='Gas 2')
         plot2.set_xlabel(pressure_key2, fontsize='x-small')
         plot2.set_ylabel(loading_key2, fontsize='x-small')
         plot2.legend(loc='best')
         plot2.set_xscale(axisScale)
         canvasFigure2.draw()
         plot2.autoscale(enable=False, axis='both', tight=None)
     
def plot_data3(): #plots the calculated selectivities vs pressure for the first mole fraction that was put in the mole fraction entryfield
    global pressures
    global stored_selectivities
    
    if stored_selectivities == []:
        return
    
    plot3.axes.clear()
    plot3.plot(pressures,
             stored_selectivities,
             'bo',
             markersize=3,
             color='r',
             label='Selectivity')
    plot3.set_xlabel(pressure_key2, fontsize='x-small')
    plot3.set_ylabel('Selectivity(-)', fontsize='x-small')
    plot3.legend(loc='best')
    plot3.set_xscale(axisScale)
    canvasFigure3.draw()

def Simpletoggle(): #changes the axis of the plots from linear to logarithmic and back
    global axisScale
    global lastSelectedModel2
    
    #from logarithmic to linear, then redo plot 1, 2 & 3 if there is data to be plotted
    if buttonAxis.config('text')[-1] == 'Linear Axis':
        buttonAxis.config(text='Logarithmic Axis')
        axisScale = 'linear'
        try:
            df_gas1
        except:
            pass
        else:
            plot_data1()
        
        try:
            df_gas2
        except:
            pass
        else:
            plot_data2()
        
        if stored_selectivities != []:
            plot_data3()
        if Fit1 != []:
            plot1.plot(df_gas1[pressure_key1],
               Fit1.loading(df_gas1[pressure_key1],),
               color= lastSelectedColour1,
               linewidth=1,
               label='Fit ' + lastSelectedModel1)
            plot1.set_xscale(axisScale)
            plot1.legend(loc='best')
            canvasFigure1.draw()
            
        if Fit2 != []:
            plot2.plot(df_gas2[pressure_key2],
               Fit2.loading(df_gas2[pressure_key2],),
               color= lastSelectedColour2,
               linewidth=1,
               label='Fit ' + lastSelectedModel2)
            plot2.legend(loc='best')
            plot2.set_xscale(axisScale)
            canvasFigure2.draw()
        else:
            return
    
    #from linear to logarithmic, then redo plot 1, 2 & 3 if there is data to be plotted    
    else:
        buttonAxis.config(text='Linear Axis')
        axisScale = 'log'
        try:
            df_gas1
        except:
            pass
        else:
            plot_data1()
        
        try:
            df_gas2
        except:
            pass
        else:
            plot_data2()
        
        if stored_selectivities != []:
            plot_data3()
        if Fit1 != []:
            plot1.plot(df_gas1[pressure_key1],
               Fit1.loading(df_gas1[pressure_key1],),
               color= lastSelectedColour1,
               linewidth=1,
               label='Fit ' + lastSelectedModel1)
            plot1.set_xscale(axisScale)
            plot1.legend(loc='best')
            canvasFigure1.draw()
            
        if Fit2 != []:
            plot2.plot(df_gas2[pressure_key2],
               Fit2.loading(df_gas2[pressure_key2],),
               color= lastSelectedColour2,
               linewidth=1,
               label='Fit ' + lastSelectedModel2)
            plot2.legend(loc='best')
            plot2.set_xscale(axisScale)
            canvasFigure2.draw()
        else:
            return

def popup_customfit1(): # popup menu for entering initial parameter guesses for model fitting for isotherm 1
    customfitwindow1 = tk.Toplevel()
    customfitwindow1.wm_title("Custom Fit 1")
    customfitwindow1.geometry("275x310")
    
    #connect variables to the main window
    global customLangmuir1M
    customLangmuir1M = tk.StringVar(customfitwindow1)
    global customLangmuir1K
    customLangmuir1K = tk.StringVar(customfitwindow1)
    global customHenry1KH
    customHenry1KH = tk.StringVar(customfitwindow1)
    global customQuadratic1M
    customQuadratic1M = tk.StringVar(customfitwindow1)
    global customQuadratic1Ka
    customQuadratic1Ka = tk.StringVar(customfitwindow1)
    global customQuadratic1Kb
    customQuadratic1Kb = tk.StringVar(customfitwindow1)
    global customBET1M
    customBET1M = tk.StringVar(customfitwindow1)
    global customBET1Ka
    customBET1Ka = tk.StringVar(customfitwindow1)
    global customBET1Kb
    customBET1Kb = tk.StringVar(customfitwindow1)
    global customDSLangmuir1M1
    customDSLangmuir1M1 = tk.StringVar(customfitwindow1)
    global customDSLangmuir1M2
    customDSLangmuir1M2 = tk.StringVar(customfitwindow1)
    global customDSLangmuir1K1
    customDSLangmuir1K1 = tk.StringVar(customfitwindow1)
    global customDSLangmuir1K2
    customDSLangmuir1K2 = tk.StringVar(customfitwindow1)
    global customTemkin1M
    customTemkin1M = tk.StringVar(customfitwindow1)
    global customTemkin1K
    customTemkin1K = tk.StringVar(customfitwindow1)
    global customTemkin1Theta
    customTemkin1Theta = tk.StringVar(customfitwindow1)
    
    customlabel = tk.Label(customfitwindow1, text="Custom Fit Parameters", font = 'bold')
    customlabel.grid(row=0, column=0, columnspan=6)
    
    placeholder0 = tk.Label(customfitwindow1, text='\n', height=2, width=5)
    placeholder0.grid(row=14, column=3)

    #button to save and export fitting parameters input by the user
    custombutton = tk.Button(customfitwindow1, text="Save", command= lambda: [customfitwindow1.destroy(), customParameterExport1()], width = 15)
    custombutton.grid(row=13, column=0, columnspan = 6)
    
    #place Langmuir fitting options
    customlabelLangmuir = tk.Label(customfitwindow1, text="Langmuir")
    customlabelLangmuir.configure(font = ('arial', 10, 'bold'))
    customlabelLangmuir.grid(row=1, column=0, columnspan=2, sticky='w')
    customlabelLangmuirM = tk.Label(customfitwindow1, text="M = ")
    customlabelLangmuirM.grid(row=2, column=0)
    entryfieldCustomFitLangmuirM = tk.Entry(customfitwindow1, textvariable = customLangmuir1M, bd=2, width=10, bg='white')
    entryfieldCustomFitLangmuirM.grid(row=2, column=1)
    entryfieldCustomFitLangmuirM.insert(0,str(customParameterarray1[0]))
    customlabelLangmuirK = tk.Label(customfitwindow1, text="K = ")
    customlabelLangmuirK.grid(row=3, column=0)
    entryfieldCustomFitLangmuirK = tk.Entry(customfitwindow1, textvariable = customLangmuir1K, bd=2, width=10, bg='white')
    entryfieldCustomFitLangmuirK.grid(row=3, column=1)
    entryfieldCustomFitLangmuirK.insert(0,str(customParameterarray1[1]))
    
    #place Henry fitting options
    customlabelHenry = tk.Label(customfitwindow1, text="Henry")
    customlabelHenry.configure(font = ('arial', 10, 'bold'))
    customlabelHenry.grid(row=1, column=4, columnspan=2, sticky='w')
    customlabelHenryKH = tk.Label(customfitwindow1, text="KH = ")
    customlabelHenryKH.grid(row=2, column=4)
    entryfieldCustomFitHenry = tk.Entry(customfitwindow1, textvariable = customHenry1KH, bd=2, width=10, bg='white')
    entryfieldCustomFitHenry.grid(row=2, column=5)
    entryfieldCustomFitHenry.insert(0,str(customParameterarray1[2]))
    
    #place Quadratic fitting options
    customlabelQuadratic = tk.Label(customfitwindow1, text="Quadratic")
    customlabelQuadratic.configure(font = ('arial', 10, 'bold'))
    customlabelQuadratic.grid(row=4, column=0, columnspan=2, sticky='w')
    customlabelQuadraticM = tk.Label(customfitwindow1, text="M = ")
    customlabelQuadraticM.grid(row=5, column=0)
    entryfieldCustomFitQuadraticM = tk.Entry(customfitwindow1, textvariable = customQuadratic1M, bd=2, width=10, bg='white')
    entryfieldCustomFitQuadraticM.grid(row=5, column=1)
    entryfieldCustomFitQuadraticM.insert(0,str(customParameterarray1[3]))
    customlabelQuadraticKa = tk.Label(customfitwindow1, text="Ka = ")
    customlabelQuadraticKa.grid(row=6, column=0)
    entryfieldCustomFitQuadraticKa = tk.Entry(customfitwindow1, textvariable = customQuadratic1Ka, bd=2, width=10, bg='white')
    entryfieldCustomFitQuadraticKa.grid(row=6, column=1)
    entryfieldCustomFitQuadraticKa.insert(0,str(customParameterarray1[4]))
    customlabelQuadraticKb = tk.Label(customfitwindow1, text="Kb = ")
    customlabelQuadraticKb.grid(row=7, column=0)
    entryfieldCustomFitQuadraticKb = tk.Entry(customfitwindow1, textvariable = customQuadratic1Kb, bd=2, width=10, bg='white')
    entryfieldCustomFitQuadraticKb.grid(row=7, column=1)
    entryfieldCustomFitQuadraticKb.insert(0,str(customParameterarray1[5]))
    
    #place BET fitting options
    customlabelBET = tk.Label(customfitwindow1, text="BET")
    customlabelBET.configure(font = ('arial', 10, 'bold'))
    customlabelBET.grid(row=4, column=4, columnspan=2, sticky='w')
    customlabelBETM = tk.Label(customfitwindow1, text="M = ")
    customlabelBETM.grid(row=5, column=4)
    entryfieldCustomFitBETM = tk.Entry(customfitwindow1, textvariable = customBET1M, bd=2, width=10, bg='white')
    entryfieldCustomFitBETM.grid(row=5, column=5)
    entryfieldCustomFitBETM.insert(0,str(customParameterarray1[6]))
    customlabelBETKa = tk.Label(customfitwindow1, text="Ka = ")
    customlabelBETKa.grid(row=6, column=4)
    entryfieldCustomFitBETKa = tk.Entry(customfitwindow1, textvariable = customBET1Ka, bd=2, width=10, bg='white')
    entryfieldCustomFitBETKa.grid(row=6, column=5)
    entryfieldCustomFitBETKa.insert(0,str(customParameterarray1[7]))
    customlabelBETKb = tk.Label(customfitwindow1, text="Kb = ")
    customlabelBETKb.grid(row=7, column=4)
    entryfieldCustomFitBETKb = tk.Entry(customfitwindow1, textvariable = customBET1Kb, bd=2, width=10, bg='white')
    entryfieldCustomFitBETKb.grid(row=7, column=5)
    entryfieldCustomFitBETKb.insert(0,str(customParameterarray1[8]))
    
    #place DSLangmuir fitting options
    customlabelDSLangmuir = tk.Label(customfitwindow1, text="DSLangmuir")
    customlabelDSLangmuir.configure(font = ('arial', 10, 'bold'))
    customlabelDSLangmuir.grid(row=8, column=0, columnspan=2, sticky='w')
    customlabelDSLangmuirM1 = tk.Label(customfitwindow1, text="M1 = ")
    customlabelDSLangmuirM1.grid(row=9, column=0)
    entryfieldCustomFitDSLangmuirM1 = tk.Entry(customfitwindow1, textvariable = customDSLangmuir1M1, bd=2, width=10, bg='white')
    entryfieldCustomFitDSLangmuirM1.grid(row=9, column=1)
    entryfieldCustomFitDSLangmuirM1.insert(0,str(customParameterarray1[9]))
    customlabelDSLangmuirM2 = tk.Label(customfitwindow1, text="M2 = ")
    customlabelDSLangmuirM2.grid(row=10, column=0)
    entryfieldCustomFitDSLangmuirM2 = tk.Entry(customfitwindow1, textvariable = customDSLangmuir1M2, bd=2, width=10, bg='white')
    entryfieldCustomFitDSLangmuirM2.grid(row=10, column=1)
    entryfieldCustomFitDSLangmuirM2.insert(0,str(customParameterarray1[10]))
    customlabelDSLangmuirK1 = tk.Label(customfitwindow1, text="K1 = ")
    customlabelDSLangmuirK1.grid(row=11, column=0)
    entryfieldCustomFitDSLangmuirK1 = tk.Entry(customfitwindow1, textvariable = customDSLangmuir1K1, bd=2, width=10, bg='white')
    entryfieldCustomFitDSLangmuirK1.grid(row=11, column=1)
    entryfieldCustomFitDSLangmuirK1.insert(0,str(customParameterarray1[11]))
    customlabelDSLangmuirK2 = tk.Label(customfitwindow1, text="K2 = ")
    customlabelDSLangmuirK2.grid(row=12, column=0)
    entryfieldCustomFitDSLangmuirK2 = tk.Entry(customfitwindow1, textvariable = customDSLangmuir1K2, bd=2, width=10, bg='white')
    entryfieldCustomFitDSLangmuirK2.grid(row=12, column=1)
    entryfieldCustomFitDSLangmuirK2.insert(0,str(customParameterarray1[12]))
    
    #place Temkin fitting options
    customlabelTemkin = tk.Label(customfitwindow1, text="Temkin")
    customlabelTemkin.configure(font = ('arial', 10, 'bold'))
    customlabelTemkin.grid(row=8, column=4, columnspan=2, sticky='w')
    customlabelTemkinM = tk.Label(customfitwindow1, text="M = ")
    customlabelTemkinM.grid(row=9, column=4)
    entryfieldCustomFitTemkinM = tk.Entry(customfitwindow1, textvariable = customTemkin1M, bd=2, width=10, bg='white')
    entryfieldCustomFitTemkinM.grid(row=9, column=5)
    entryfieldCustomFitTemkinM.insert(0,str(customParameterarray1[13]))
    customlabelTemkinK = tk.Label(customfitwindow1, text="K = ")
    customlabelTemkinK.grid(row=10, column=4)
    entryfieldCustomFitTemkinK = tk.Entry(customfitwindow1, textvariable = customTemkin1K, bd=2, width=10, bg='white')
    entryfieldCustomFitTemkinK.grid(row=10, column=5)
    entryfieldCustomFitTemkinK.insert(0,str(customParameterarray1[14]))
    customlabelTemkinTheta = tk.Label(customfitwindow1, text="Theta = ")
    customlabelTemkinTheta.grid(row=11, column=4)
    entryfieldCustomFitTemkinTheta = tk.Entry(customfitwindow1, textvariable = customTemkin1Theta, bd=2, width=10, bg='white')
    entryfieldCustomFitTemkinTheta.grid(row=11, column=5)
    entryfieldCustomFitTemkinTheta.insert(0,str(customParameterarray1[15]))

def popup_customfit2(): # popup menu for entering initial parameter guesses for model fitting for isotherm 2
    customfitwindow2 = tk.Toplevel()
    customfitwindow2.wm_title("Custom Fit 2")
    customfitwindow2.geometry("275x310")
    
    #connect variables to the main window
    global customLangmuir2M
    customLangmuir2M = tk.StringVar(customfitwindow2)
    global customLangmuir2K
    customLangmuir2K = tk.StringVar(customfitwindow2)
    global customHenry2KH
    customHenry2KH = tk.StringVar(customfitwindow2)
    global customQuadratic2M
    customQuadratic2M = tk.StringVar(customfitwindow2)
    global customQuadratic2Ka
    customQuadratic2Ka = tk.StringVar(customfitwindow2)
    global customQuadratic2Kb
    customQuadratic2Kb = tk.StringVar(customfitwindow2)
    global customBET2M
    customBET2M = tk.StringVar(customfitwindow2)
    global customBET2Ka
    customBET2Ka = tk.StringVar(customfitwindow2)
    global customBET2Kb
    customBET2Kb = tk.StringVar(customfitwindow2)
    global customDSLangmuir2M1
    customDSLangmuir2M1 = tk.StringVar(customfitwindow2)
    global customDSLangmuir2M2
    customDSLangmuir2M2 = tk.StringVar(customfitwindow2)
    global customDSLangmuir2K1
    customDSLangmuir2K1 = tk.StringVar(customfitwindow2)
    global customDSLangmuir2K2
    customDSLangmuir2K2 = tk.StringVar(customfitwindow2)
    global customTemkin2M
    customTemkin2M = tk.StringVar(customfitwindow2)
    global customTemkin2K
    customTemkin2K = tk.StringVar(customfitwindow2)
    global customTemkin2Theta
    customTemkin2Theta = tk.StringVar(customfitwindow2)

    customlabel = tk.Label(customfitwindow2, text="Custom Fit Parameters", font = 'bold')
    customlabel.grid(row=0, column=0, columnspan=6)
    
    placeholder0 = tk.Label(customfitwindow2, text='\n', height=2, width=5)
    placeholder0.grid(row=14, column=3)

    #button to save and export fitting parameters input by the user
    custombutton = tk.Button(customfitwindow2, text="Save", command= lambda: [customfitwindow2.destroy(), customParameterExport2()], width = 15)
    custombutton.grid(row=13, column=0, columnspan = 6)
    
    #place Langmuir fitting options
    customlabelLangmuir = tk.Label(customfitwindow2, text="Langmuir")
    customlabelLangmuir.configure(font = ('arial', 10, 'bold'))
    customlabelLangmuir.grid(row=1, column=0, columnspan=2, sticky='w')
    customlabelLangmuirM = tk.Label(customfitwindow2, text="M = ")
    customlabelLangmuirM.grid(row=2, column=0)
    entryfieldCustomFitLangmuirM = tk.Entry(customfitwindow2, textvariable = customLangmuir2M, bd=2, width=10, bg='white')
    entryfieldCustomFitLangmuirM.grid(row=2, column=1)
    entryfieldCustomFitLangmuirM.insert(0,str(customParameterarray2[0]))
    customlabelLangmuirK = tk.Label(customfitwindow2, text="K = ")
    customlabelLangmuirK.grid(row=3, column=0)
    entryfieldCustomFitLangmuirK = tk.Entry(customfitwindow2, textvariable = customLangmuir2K, bd=2, width=10, bg='white')
    entryfieldCustomFitLangmuirK.grid(row=3, column=1)
    entryfieldCustomFitLangmuirK.insert(0,str(customParameterarray2[1]))
    
    #place Henry fitting options
    customlabelHenry = tk.Label(customfitwindow2, text="Henry")
    customlabelHenry.configure(font = ('arial', 10, 'bold'))
    customlabelHenry.grid(row=1, column=4, columnspan=2, sticky='w')
    customlabelHenryKH = tk.Label(customfitwindow2, text="KH = ")
    customlabelHenryKH.grid(row=2, column=4)
    entryfieldCustomFitHenry = tk.Entry(customfitwindow2, textvariable = customHenry2KH, bd=2, width=10, bg='white')
    entryfieldCustomFitHenry.grid(row=2, column=5)
    entryfieldCustomFitHenry.insert(0,str(customParameterarray2[2]))
    
    #place Quadratic fitting options
    customlabelQuadratic = tk.Label(customfitwindow2, text="Quadratic")
    customlabelQuadratic.configure(font = ('arial', 10, 'bold'))
    customlabelQuadratic.grid(row=4, column=0, columnspan=2, sticky='w')
    customlabelQuadraticM = tk.Label(customfitwindow2, text="M = ")
    customlabelQuadraticM.grid(row=5, column=0)
    entryfieldCustomFitQuadraticM = tk.Entry(customfitwindow2, textvariable = customQuadratic2M, bd=2, width=10, bg='white')
    entryfieldCustomFitQuadraticM.grid(row=5, column=1)
    entryfieldCustomFitQuadraticM.insert(0,str(customParameterarray2[3]))
    customlabelQuadraticKa = tk.Label(customfitwindow2, text="Ka = ")
    customlabelQuadraticKa.grid(row=6, column=0)
    entryfieldCustomFitQuadraticKa = tk.Entry(customfitwindow2, textvariable = customQuadratic2Ka, bd=2, width=10, bg='white')
    entryfieldCustomFitQuadraticKa.grid(row=6, column=1)
    entryfieldCustomFitQuadraticKa.insert(0,str(customParameterarray2[4]))
    customlabelQuadraticKb = tk.Label(customfitwindow2, text="Kb = ")
    customlabelQuadraticKb.grid(row=7, column=0)
    entryfieldCustomFitQuadraticKb = tk.Entry(customfitwindow2, textvariable = customQuadratic2Kb, bd=2, width=10, bg='white')
    entryfieldCustomFitQuadraticKb.grid(row=7, column=1)
    entryfieldCustomFitQuadraticKb.insert(0,str(customParameterarray2[5]))
    
    #place BET fitting options
    customlabelBET = tk.Label(customfitwindow2, text="BET")
    customlabelBET.configure(font = ('arial', 10, 'bold'))
    customlabelBET.grid(row=4, column=4, columnspan=2, sticky='w')
    customlabelBETM = tk.Label(customfitwindow2, text="M = ")
    customlabelBETM.grid(row=5, column=4)
    entryfieldCustomFitBETM = tk.Entry(customfitwindow2, textvariable = customBET2M, bd=2, width=10, bg='white')
    entryfieldCustomFitBETM.grid(row=5, column=5)
    entryfieldCustomFitBETM.insert(0,str(customParameterarray2[6]))
    customlabelBETKa = tk.Label(customfitwindow2, text="Ka = ")
    customlabelBETKa.grid(row=6, column=4)
    entryfieldCustomFitBETKa = tk.Entry(customfitwindow2, textvariable = customBET2Ka, bd=2, width=10, bg='white')
    entryfieldCustomFitBETKa.grid(row=6, column=5)
    entryfieldCustomFitBETKa.insert(0,str(customParameterarray2[7]))
    customlabelBETKb = tk.Label(customfitwindow2, text="Kb = ")
    customlabelBETKb.grid(row=7, column=4)
    entryfieldCustomFitBETKb = tk.Entry(customfitwindow2, textvariable = customBET2Kb, bd=2, width=10, bg='white')
    entryfieldCustomFitBETKb.grid(row=7, column=5)
    entryfieldCustomFitBETKb.insert(0,str(customParameterarray2[8]))
    
    #place DSLangmuir fitting options
    customlabelDSLangmuir = tk.Label(customfitwindow2, text="DSLangmuir")
    customlabelDSLangmuir.configure(font = ('arial', 10, 'bold'))
    customlabelDSLangmuir.grid(row=8, column=0, columnspan=2, sticky='w')
    customlabelDSLangmuirM1 = tk.Label(customfitwindow2, text="M1 = ")
    customlabelDSLangmuirM1.grid(row=9, column=0)
    entryfieldCustomFitDSLangmuirM1 = tk.Entry(customfitwindow2, textvariable = customDSLangmuir2M1, bd=2, width=10, bg='white')
    entryfieldCustomFitDSLangmuirM1.grid(row=9, column=1)
    entryfieldCustomFitDSLangmuirM1.insert(0,str(customParameterarray2[9]))
    customlabelDSLangmuirM2 = tk.Label(customfitwindow2, text="M2 = ")
    customlabelDSLangmuirM2.grid(row=10, column=0)
    entryfieldCustomFitDSLangmuirM2 = tk.Entry(customfitwindow2, textvariable = customDSLangmuir2M2, bd=2, width=10, bg='white')
    entryfieldCustomFitDSLangmuirM2.grid(row=10, column=1)
    entryfieldCustomFitDSLangmuirM2.insert(0,str(customParameterarray2[10]))
    customlabelDSLangmuirK1 = tk.Label(customfitwindow2, text="K1 = ")
    customlabelDSLangmuirK1.grid(row=11, column=0)
    entryfieldCustomFitDSLangmuirK1 = tk.Entry(customfitwindow2, textvariable = customDSLangmuir2K1, bd=2, width=10, bg='white')
    entryfieldCustomFitDSLangmuirK1.grid(row=11, column=1)
    entryfieldCustomFitDSLangmuirK1.insert(0,str(customParameterarray2[11]))
    customlabelDSLangmuirK2 = tk.Label(customfitwindow2, text="K2 = ")
    customlabelDSLangmuirK2.grid(row=12, column=0)
    entryfieldCustomFitDSLangmuirK2 = tk.Entry(customfitwindow2, textvariable = customDSLangmuir2K2, bd=2, width=10, bg='white')
    entryfieldCustomFitDSLangmuirK2.grid(row=12, column=1)
    entryfieldCustomFitDSLangmuirK2.insert(0,str(customParameterarray2[12]))
    
    #place Temkin fitting options
    customlabelTemkin = tk.Label(customfitwindow2, text="Temkin")
    customlabelTemkin.configure(font = ('arial', 10, 'bold'))
    customlabelTemkin.grid(row=8, column=4, columnspan=2, sticky='w')
    customlabelTemkinM = tk.Label(customfitwindow2, text="M = ")
    customlabelTemkinM.grid(row=9, column=4)
    entryfieldCustomFitTemkinM = tk.Entry(customfitwindow2, textvariable = customTemkin2M, bd=2, width=10, bg='white')
    entryfieldCustomFitTemkinM.grid(row=9, column=5)
    entryfieldCustomFitTemkinM.insert(0,str(customParameterarray2[13]))
    customlabelTemkinK = tk.Label(customfitwindow2, text="K = ")
    customlabelTemkinK.grid(row=10, column=4)
    entryfieldCustomFitTemkinK = tk.Entry(customfitwindow2, textvariable = customTemkin2K, bd=2, width=10, bg='white')
    entryfieldCustomFitTemkinK.grid(row=10, column=5)
    entryfieldCustomFitTemkinK.insert(0,str(customParameterarray2[14]))
    customlabelTemkinTheta = tk.Label(customfitwindow2, text="Theta = ")
    customlabelTemkinTheta.grid(row=11, column=4)
    entryfieldCustomFitTemkinTheta = tk.Entry(customfitwindow2, textvariable = customTemkin2Theta, bd=2, width=10, bg='white')
    entryfieldCustomFitTemkinTheta.grid(row=11, column=5)
    entryfieldCustomFitTemkinTheta.insert(0,str(customParameterarray2[15]))

def customParameterExport1(): #exports the custom fitting parameters for fit 1
    global customParameterarray1
    
    inputLangmuirM = customLangmuir1M.get()
    inputLangmuirK = customLangmuir1K.get()
    inputHenryKH = customHenry1KH.get()
    inputQuadraticM = customQuadratic1M.get()
    inputQuadraticKa = customQuadratic1Ka.get()
    inputQuadraticKb = customQuadratic1Kb.get()
    inputBETM = customBET1M.get()
    inputBETKa = customBET1Ka.get()
    inputBETKb = customBET1Kb.get()
    inputDSLangmuirM1 = customDSLangmuir1M1.get()
    inputDSLangmuirM2 = customDSLangmuir1M2.get()
    inputDSLangmuirK1 = customDSLangmuir1K1.get()
    inputDSLangmuirK2 = customDSLangmuir1K2.get()
    inputTemkinM = customTemkin1M.get()
    inputTemkinK = customTemkin1K.get()
    inputTemkinTheta = customTemkin1Theta.get()
    customParameterarray1 = [inputLangmuirM, inputLangmuirK, inputHenryKH, inputQuadraticM, inputQuadraticKa, inputQuadraticKb,
                             inputBETM, inputBETKa, inputBETKb, inputDSLangmuirM1, inputDSLangmuirM2, inputDSLangmuirK1,
                             inputDSLangmuirK2, inputTemkinM, inputTemkinK, inputTemkinTheta]
    for i in range(0, len(customParameterarray1)):
        if customParameterarray1[i] == '':
            i += 1
        else:
            try:
                float(customParameterarray1[i])
            except:
                error = "Wrongful input in custom fit 1 entryfield"
                popup_errormsg(error)
                exit 
            i += 1

def customParameterExport2(): #exports the custom fitting parameters for fit 2
    global customParameterarray2
    
    inputLangmuirM = customLangmuir2M.get()
    inputLangmuirK = customLangmuir2K.get()
    inputHenryKH = customHenry2KH.get()
    inputQuadraticM = customQuadratic2M.get()
    inputQuadraticKa = customQuadratic2Ka.get()
    inputQuadraticKb = customQuadratic2Kb.get()
    inputBETM = customBET2M.get()
    inputBETKa = customBET2Ka.get()
    inputBETKb = customBET2Kb.get()
    inputDSLangmuirM1 = customDSLangmuir2M1.get()
    inputDSLangmuirM2 = customDSLangmuir2M2.get()
    inputDSLangmuirK1 = customDSLangmuir2K1.get()
    inputDSLangmuirK2 = customDSLangmuir2K2.get()
    inputTemkinM = customTemkin2M.get()
    inputTemkinK = customTemkin2K.get()
    inputTemkinTheta = customTemkin2Theta.get()
    customParameterarray2 = [inputLangmuirM, inputLangmuirK, inputHenryKH, inputQuadraticM, inputQuadraticKa, inputQuadraticKb,
                             inputBETM, inputBETKa, inputBETKb, inputDSLangmuirM1, inputDSLangmuirM2, inputDSLangmuirK1,
                             inputDSLangmuirK2, inputTemkinM, inputTemkinK, inputTemkinTheta]
    for i in range(0, len(customParameterarray2)):
        if customParameterarray2[i] == '':
            i += 1
        else:
            try:
                float(customParameterarray2[i])
            except:
                error = "Wrongful input in custom fit 2 entryfield"
                popup_errormsg(error)
                exit 
            i += 1 

def clearCustomFit1(): #clears the custom fit parameter 1 input by settting it to empty strings
    global customParameterarray1
    
    customParameterarray1 =['','','','','','','','','','','','','','','','']

def clearCustomFit2(): #clears the custom fit parameter 2 input by settting it to empty strings
    global customParameterarray2
    
    customParameterarray2 =['','','','','','','','','','','','','','','','']
    
def safe_file(): #opens a window where the user can input a filename and location to save the data. the data in then stored in a table separated by commas and exported as a csv file.
    
    file = filedialog.asksaveasfile(
            defaultextension='.csv',
            filetypes=[
                    ("CSV file", ".csv"),
                    ("Text file", ".txt"),
                    ("All files", ".*"),])
    if file is None:                        
       return
    line = ''
    filetext = str('Total Pressure [bar], Mole Fraction Gas 1, Mole Fraction Gas 2, Partial Pressure Gas 1, Partial Pressure Gas 2, loading Gas 1 [mmol/g], loading Gas 2 [mmol/g], Selectivity \n')
    for i in range(0, len(selectivities)):
        for j in range(0, 8):
            line += str(selectivities[i][j]) + ','
            j += 1
        line += '\n'
        i += 1
        filetext += line
        line = ''
    file.write(filetext)
    file.close()

def popup_errormsg(errormessage): #opens a window with an error message to inform the user
    tk.messagebox.showwarning(title=None, message=errormessage)

def display_help(): #opens a window with the manual. 

    helpwindow = tk.Toplevel()
    image = Image.open("Help_image.png")
    w, h = image.size
    photo = ImageTk.PhotoImage(image)
    helpwindow.photo = photo  # solution for bug in `PhotoImage`
    helpcanvas = tk.Canvas(helpwindow, width=w, height=h)
    helpcanvas.pack()
    helpcanvas.create_image(0, 0, anchor='nw', image=photo)

def close_window():
    window.destroy()
    root.destroy()

#define class for scrollable window_______________________________
class ScrollbarFrame(tk.Frame):
    """
    Extends class tk.Frame to support a scrollable Frame 
    This class is independent from the widgets to be scrolled and 
    can be used to replace a standard tk.Frame
    """
    def __init__(self, parent, **kwargs):
        tk.Frame.__init__(self, parent, **kwargs)

        # The Scrollbar, layout to the right
        vsb = tk.Scrollbar(self, orient="vertical")
        vsb.pack(side="right", fill="y")
        hsb = tk.Scrollbar(self, orient="horizontal")
        hsb.pack(side="bottom", fill="x")

        # The Canvas which supports the Scrollbar Interface, layout to the left
        self.canvas = tk.Canvas(self, borderwidth=0, background="#ffffff")
        self.canvas.pack(side="left", fill="both", expand=True)

        # Bind the Scrollbar to the self.canvas Scrollbar Interface
        self.canvas.configure(yscrollcommand=vsb.set)
        vsb.configure(command=self.canvas.yview)
        self.canvas.configure(xscrollcommand=hsb.set)
        hsb.configure(command=self.canvas.xview)

        # The Frame to be scrolled, layout into the canvas
        # All widgets to be scrolled have to use this Frame as parent
        self.scrolled_frame = tk.Frame(self.canvas, background=self.canvas.cget('bg'))
        self.canvas.create_window((4, 4), window=self.scrolled_frame, anchor="nw")

        # Configures the scrollregion of the Canvas dynamically
        self.scrolled_frame.bind("<Configure>", self.on_configure)

    def on_configure(self, event):
        """Set the scroll region to encompass the scrolled frame"""
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))
#_______________________________________________
#general GUI window settings

root = tk.Tk()
root.geometry("1220x850")  # Size of the window 
root.title('GraphIAST')
sbf = ScrollbarFrame(root)
root.grid_rowconfigure(0, weight=1)
root.grid_columnconfigure(0, weight=1)
sbf.grid(row=0, column=0, sticky='nsew')
window = sbf.scrolled_frame
window.configure(bg='white')
my_font1=('arial', 10, 'bold')

var = tk.IntVar()
var.set(0)

root.protocol("WM_DELETE_WINDOW", close_window)
#______________________________________________________________________________
#import data button 1
labelimport1 = tk.Label(window,text='Import Isotherm 1',width=30,font=my_font1, bg='white')  
labelimport1.grid(row=1,column=3)
buttonimport1 = tk.Button(window, text='Import File', 
   width=20,command = lambda:[import_file1(),plot_data1()])
buttonimport1.grid(row=2,column=3) 

#import data button 2
labelimport2 = tk.Label(window,text='Import Isotherm 2',width=30,font=my_font1, bg='white')  
labelimport2.grid(row=1,column=4)
buttonimport2 = tk.Button(window, text='Import File', 
   width=20,command = lambda:[import_file2(),plot_data2()])
buttonimport2.grid(row=2,column=4) 
    
#Fit button 1
labelFit1 = tk.Label(window,text='Fit Isotherms',width=30,font=my_font1, bg='white')  
labelFit1.grid(row=10,column=4)
buttonFit1 = tk.Button(window, text='Fit Isotherm 1', 
   width=20,command = lambda:fit_isotherm1())
buttonFit1.grid(row=11,column=4) 

#Fit button 2
buttonFit2 = tk.Button(window, text='Fit Isotherm 2', 
   width=20,command = lambda:fit_isotherm2())
buttonFit2.grid(row=12,column=4) 

#Safe button
labelSave = tk.Label(window,text='Save Data in File',width=30,font=my_font1, bg='white')  
labelSave.grid(row=28,column=3)
buttonSave = tk.Button(window, text='Save', 
   width=20,command = lambda:safe_file())
buttonSave.grid(row=29,column=3) 

#Quit button
buttonQuit = tk.Button(window, text='Quit', 
   width=20,command = close_window)
buttonQuit.grid(row=29,column=4) 

#custom fit parameter button 1 & 2
labelFit2 = tk.Label(window,text='Define Fit Parameters',width=30,font=my_font1, bg='white')  
labelFit2.grid(row=7,column=4)
buttonFitParameters1 = tk.Button(window, text='Custom Fit 1', 
   width=20,command = popup_customfit1)
buttonFitParameters1.grid(row=8,column=4)

buttonFitParameters2 = tk.Button(window, text='Custom Fit 2', 
   width=20,command = popup_customfit2)
buttonFitParameters2.grid(row=9,column=4)

#custom fit parameter clear button 1 & 2
buttonClearCustomFit1 = tk.Button(window, text='Clear', 
   width=5,command = clearCustomFit1)
buttonClearCustomFit1.grid(row=8,column=5)

buttonClearCustomFit2 = tk.Button(window, text='Clear', 
   width=5,command = clearCustomFit2)
buttonClearCustomFit2.grid(row=9,column=5)

#help button
buttonHelp = tk.Button(window, text='Help', 
   width=5,command = display_help)
buttonHelp.grid(row=29,column=5)

#Run IAST button
labelIAST = tk.Label(window,text='Calculate IAST Selectivities',width=30,font=my_font1, bg='white')  
labelIAST.grid(row=25,columnspan=2, column=3)
buttonIAST = tk.Button(window, text='Run IAST', 
   width=20,command = lambda:[run_IAST(), plot_data3()])
buttonIAST.grid(row=26,columnspan=2, column=3) 

#Change Axis button
labelAxis = tk.Label(window, text='Switch x Axis', width=30, font=my_font1, bg='white')
labelAxis.grid(row=13, column=4)
buttonAxis = tk.Button(window, text="Logarithmic Axis", width=20, command = lambda:Simpletoggle())
buttonAxis.grid(row=14, column=4)

#______________________________________________________________________________
#coupled radio buttons for model selection

labelModels = tk.Label(window, text='Select Model for Fit', width=30, font=my_font1, bg='white')
labelModels.grid(row=7, column=3)

tk.Radiobutton(window, text="Interpolator", variable=var, value=0, bg='white').grid(row = 8, column = 3, sticky='W',ipadx=100)
tk.Radiobutton(window, text="Langmuir", variable=var, value=1, bg='white').grid(row = 9, column = 3, sticky='W',ipadx=100)
tk.Radiobutton(window, text="Quadratic", variable=var, value=2, bg='white').grid(row = 10, column = 3, sticky='W',ipadx=100)
tk.Radiobutton(window, text="Henry", variable=var, value=3, bg='white').grid(row = 11, column = 3, sticky='W',ipadx=100)
tk.Radiobutton(window, text="BET", variable=var, value=4, bg='white').grid(row = 12, column = 3, sticky='W',ipadx=100)
tk.Radiobutton(window, text="Temkin", variable=var, value=5, bg='white').grid(row = 13, column = 3, sticky='W',ipadx=100)
tk.Radiobutton(window, text="Dual Side Langmuir", variable=var, value=6, bg='white').grid(row = 14, column = 3, sticky='W',ipadx=100)

#______________________________________________________________________________
#entryfield for pressures
labelPressures = tk.Label(window, text='Define Pressures for IAST [bar] \n (i.e. 0.1, 1.0,...)', width=30, font=my_font1, bg='white').grid(row=21, column=3)
entryfieldPressures = tk.Entry(window, bd=2, width=30, bg='white')
entryfieldPressures.grid(row=22, column=3)

#entryfield for mole fractions
labelMoleFractions = tk.Label(window, text='Define Mole Fractions Gas 1 for IAST \n (i.e. 0.05, 0.1,...)', width=30, font=my_font1, bg='white').grid(row=21, column=4)
entryfieldMoleFractions = tk.Entry(window, bd=2, width=30, bg='white')
entryfieldMoleFractions.grid(row=22, column=4)

#______________________________________________________________________________
#fitting parameters textboxes
labelParameters1 = tk.Label(window, text='               Fitting Parameters', bg='white', width=30, font=my_font1).grid(row=13, column=1)
labelParameters2 = tk.Label(window, text='               Fitting Parameters', bg='white', width=30, font=my_font1).grid(row=13, column=2)

textParameters1 = scrolledtext.ScrolledText(window, height=7, width=29, font=my_font1, bg='white', undo=True)
textParameters1.configure(state='disabled')
textParameters1.grid(rowspan=5, row=14, column=1, sticky='E')

textParameters2 = scrolledtext.ScrolledText(window, height=7, width=29, font=my_font1, bg='white', undo=True)
textParameters2.grid(rowspan=5, row=14, column=2, sticky='E')
textParameters2.configure(state='disabled')

#______________________________________________________________________________
#plots

#plot1
labelFigure1 = tk.Label(window, text='               Gas 1 including Fit', bg='white', width=30, font=my_font1).grid(row=1, column=1)
Figure1 = Figure(figsize=(2.75,2.75), dpi=100, constrained_layout=True)
plot1 = Figure1.add_subplot(111)
canvasFigure1 = FigureCanvasTkAgg(Figure1, master=window)
canvasFigure1.draw()
canvasFigure1.get_tk_widget().grid(column=1, rowspan=10, row=2)

#plot2
labelFigure2 = tk.Label(window, text='               Gas 2 including Fit', bg='white', width=30, font=my_font1).grid(row=1, column=2)
Figure2 = Figure(figsize=(2.75,2.75), dpi=100, constrained_layout=True)
plot2 = Figure2.add_subplot(111)
canvasFigure2 = FigureCanvasTkAgg(Figure2, master=window)
canvasFigure2.draw()
canvasFigure2.get_tk_widget().grid(column=2, rowspan=10, row=2)

#plot3
placeholder4 = tk.Label(window, text='\r\n', bg='white', height=2, width=30)
placeholder4.grid(row=19, columnspan=2, column=1)
labelFigure3 = tk.Label(window, text='           IAST Selectivities vs. Pressure', bg='white', width=30, font=my_font1).grid(row=20, columnspan=2, column=1)
Figure3 = Figure(figsize=(2.75,2.75), dpi=100, constrained_layout=True)
plot3 = Figure3.add_subplot(111)
canvasFigure3 = FigureCanvasTkAgg(Figure3, master=window)
canvasFigure3.draw()
canvasFigure3.get_tk_widget().grid(columnspan=2, column=1, rowspan=10, row=21)

#______________________________________________________________________________

root.mainloop()  #Keep the window open