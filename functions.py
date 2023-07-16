import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

def plot_processed(title,conc,analyte,position):
    b = pd.read_csv('b.csv')
    c = pd.read_csv('c.csv')
    d = pd.read_csv('d.csv')
    e = pd.read_csv('e.csv')
    f = pd.read_csv('f.csv')
    g = pd.read_csv('g.csv')
    h = pd.read_csv('h.csv')  
    
    b['time'] = b['Time1'] - b['Time1'][0]
    c['time'] = c['Time1'] - c['Time1'][0]
    d['time'] = d['Time1'] - d['Time1'][0]
    e['time'] = e['Time1'] - e['Time1'][0]
    f['time'] = f['Time1'] - f['Time1'][0]
    g['time'] = g['Time1'] - g['Time1'][0]
    h['time'] = h['Time1'] - h['Time1'][0]
    


    plt.rcParams['figure.figsize'] = [12, 12]
    plt.title(title,fontdict={'fontsize':16})
    plt.xlabel('Time (s)',fontsize=14)
    plt.ylabel('Wavelength Shift (nm)',fontsize=14)
    plt.plot(b['time'],b['Data1'],'#88CCEE',label = '{} nM {}'.format(conc[0],analyte))
    plt.plot(c['time'],c['Data1'],'#44AA99',label = '{} nM {}'.format(conc[1],analyte))
    plt.plot(d['time'],d['Data1'],'#117733',label = '{} nM {}'.format(conc[2],analyte))
    plt.plot(e['time'],e['Data1'],'#999933',label = '{} nM {}'.format(conc[3],analyte))
    plt.plot(f['time'],f['Data1'],'#DDCC77',label = '{} nM {}'.format(conc[4],analyte))
    plt.plot(g['time'],g['Data1'],'#CC6677',label = '{} nM {}'.format(conc[5],analyte))
    plt.plot(h['time'],h['Data1'],'#882256',label = '{} nM {}'.format(conc[6],analyte))
    plt.legend(loc=position,prop={'size':15})
    plt.show()


def sep_assoc_dissoc():
    b = pd.read_csv('b.csv')
    c = pd.read_csv('c.csv')
    d = pd.read_csv('d.csv')
    e = pd.read_csv('e.csv')
    f = pd.read_csv('f.csv')
    g = pd.read_csv('g.csv')
    h = pd.read_csv('h.csv')  
    
    b['time'] = b['Time1'] - b['Time1'][0]
    
    assoc_time = []
    z = 0
    for i in range(2501):
        assoc_time.append(b['time'][z])
        z+=1 
    dis_time = []
    z = 2501
    for i in range(2499):
        dis_time.append(round(b['time'][z]-500.2,2))
        z+=1
    b_assoc_data = []
    z = 0
    for i in range(2501):
        b_assoc_data.append(b['Data1'][z])
        z+=1  
    b_dis_data = []
    z = 2501
    for i in range(2499):
        b_dis_data.append(b['Data1'][z])
        z+=1
    c_assoc_data = []
    z = 0
    for i in range(2501):
        c_assoc_data.append(c['Data1'][z])
        z+=1    
    c_dis_data = []
    z = 2501
    for i in range(2499):
        c_dis_data.append(c['Data1'][z])
        z+=1
    d_assoc_data = []
    z = 0
    for i in range(2501):
        d_assoc_data.append(d['Data1'][z])
        z+=1
    d_dis_data = []
    z = 2501
    for i in range(2499):
        d_dis_data.append(d['Data1'][z])
        z+=1
    e_assoc_data = []
    z = 0
    for i in range(2501):
        e_assoc_data.append(e['Data1'][z])
        z+=1 
    e_dis_data = []
    z = 2501
    for i in range(2499):
        e_dis_data.append(e['Data1'][z])
        z+=1
    f_assoc_data = []
    z = 0
    for i in range(2501):
        f_assoc_data.append(f['Data1'][z])
        z+=1
    f_dis_data = []
    z = 2501
    for i in range(2499):
        f_dis_data.append(f['Data1'][z])
        z+=1
    g_assoc_data = []
    z = 0
    for i in range(2501):
        g_assoc_data.append(g['Data1'][z])
        z+=1
    g_dis_data = []
    z = 2501
    for i in range(2499):
        g_dis_data.append(g['Data1'][z])
        z+=1
    h_assoc_data = []
    z = 0
    for i in range(2501):
        h_assoc_data.append(h['Data1'][z])
        z+=1
    h_dis_data = []
    z = 2501
    for i in range(2499):
        h_dis_data.append(h['Data1'][z])
        z+=1
    return (
        assoc_time,b_assoc_data,c_assoc_data,d_assoc_data,e_assoc_data,f_assoc_data,g_assoc_data,h_assoc_data,
        dis_time,b_dis_data,c_dis_data,d_dis_data,e_dis_data,f_dis_data,g_dis_data,h_dis_data
            )
def make_df(t,b,c,d,e,f,g,h):
    df = pd.DataFrame()
    df['time'] = t
    df['b_data'] = b
    df['c_data'] = c
    df['d_data'] = d
    df['e_data'] = e
    df['f_data'] = f
    df['g_data'] = g
    df['h_data'] = h
    return df

def flip_df(df):
    df_flip = pd.DataFrame()
    df_flip['time'] = df['time']
    df_flip['b_data'] = df['b_data'] * -1
    df_flip['c_data'] = df['c_data'] * -1
    df_flip['d_data'] = df['d_data'] * -1
    df_flip['e_data'] = df['e_data'] * -1
    df_flip['f_data'] = df['f_data'] * -1
    df_flip['g_data'] = df['g_data'] * -1
    df_flip['h_data'] = df['h_data'] * -1
    return df_flip


def plot_df(df,title,conc,analyte,position):
    plt.rcParams['figure.figsize'] = [14, 12]
    plt.title(title,fontdict={'fontsize':16})
    plt.xlabel('Time (s)',fontsize=14)
    plt.ylabel('Wavelength Shift (nm)',fontsize=14)
    plt.plot(df['time'],df['b_data'],'#88CCEE',label = '{} nM {}'.format(conc[0],analyte),marker='o')
    plt.plot(df['time'],df['c_data'],'#44AA99',label = '{} nM {}'.format(conc[1],analyte),marker='o')
    plt.plot(df['time'],df['d_data'],'#117733',label = '{} nM {}'.format(conc[2],analyte),marker='o')
    plt.plot(df['time'],df['e_data'],'#999933',label = '{} nM {}'.format(conc[3],analyte),marker='o')
    plt.plot(df['time'],df['f_data'],'#DDCC77',label = '{} nM {}'.format(conc[4],analyte),marker='o')
    plt.plot(df['time'],df['g_data'],'#CC6677',label = '{} nM {}'.format(conc[5],analyte),marker='o')
    plt.plot(df['time'],df['h_data'],'#882256',label = '{} nM {}'.format(conc[6],analyte),marker='o')

    plt.legend(loc=position,prop={'size':15})


def subtract_first(df):
    df_sub = pd.DataFrame()
    df_sub['time'] = df['time']
    df_sub['b_data'] = df['b_data']
    df_sub['c_data'] = abs(df['c_data'] - df['b_data'])
    df_sub['d_data'] = abs(df['d_data'] - df['b_data'])
    df_sub['e_data'] = abs(df['e_data'] - df['b_data'])
    df_sub['f_data'] = abs(df['f_data'] - df['b_data'])
    df_sub['g_data'] = abs(df['g_data'] - df['b_data'])
    df_sub['h_data'] = abs(df['h_data'] - df['b_data'])
    return df_sub


def plot_raw(title,conc,analyte):
    data = pd.read_csv('RawData.csv')
    plt.rcParams['figure.figsize'] = [10, 8]
    plt.title(title,fontdict={'fontsize':16})
    plt.xlabel('Time (s)',fontsize=14)
    plt.ylabel('Wavelength Shift (nm)',fontsize=14)
    end = 7500

    plt.plot(data['A_time'][0:end:],data['A_signal'][0:end:],'#332288',label='Reference Buffer')
    plt.plot(data['B_time'][0:end:],data['B_signal'][0:end:],'#88CCEE',label = '{}nM {}'.format(conc[0],analyte))
    plt.plot(data['C_time'][0:end:],data['C_signal'][0:end:],'#44AA99',label = '{}nM {}'.format(conc[1],analyte))
    plt.plot(data['D_time'][0:end:],data['D_signal'][0:end:],'#117733',label = '{}nM {}'.format(conc[2],analyte))
    plt.plot(data['E_time'][0:end:],data['E_signal'][0:end:],'#999933',label = '{}nM {}'.format(conc[3],analyte))
    plt.plot(data['F_time'][0:end:],data['F_signal'][0:end:],'#DDCC77',label = '{}nM {}'.format(conc[4],analyte))
    plt.plot(data['G_time'][0:end:],data['G_signal'][0:end:],'#CC6677',label = '{}nM {}'.format(conc[5],analyte))
    plt.plot(data['H_time'][0:end:],data['H_signal'][0:end:],'#882256',label = '{}nM {}'.format(conc[6],analyte))

    plt.legend(loc='upper left',prop={'size':10})
    
def kinetic_analysis(df,type,conc,analyte,lip_type,lip_conc,position):
    if type == 'association':
        def Gauss(x, A, B):
            y = A*(1-np.exp(-B*x/100))
            return y
    elif type == 'dissociation':
        def Gauss(x, A, B):
            y = A*np.exp(-B*x/100)
            return y
    param_b, cov_b = curve_fit(Gauss,df['time'],df['b_data'])
    fit_asymp_b = param_b[0]
    fit_kobs_b = param_b[1]
    fit_y_b = Gauss(df['time'], fit_asymp_b, fit_kobs_b)

    param_c, cov_c = curve_fit(Gauss,df['time'],df['c_data'])
    fit_asymp_c = param_c[0]
    fit_kobs_c = param_c[1]
    fit_y_c = Gauss(df['time'], fit_asymp_c, fit_kobs_c)

    param_d, cov_d = curve_fit(Gauss,df['time'],df['d_data'])
    fit_asymp_d = param_d[0]
    fit_kobs_d = param_d[1]
    fit_y_d = Gauss(df['time'], fit_asymp_d, fit_kobs_d)

    param_e, cov_e = curve_fit(Gauss,df['time'],df['e_data'])
    fit_asymp_e = param_e[0]
    fit_kobs_e = param_e[1]
    fit_y_e = Gauss(df['time'], fit_asymp_e, fit_kobs_e)

    param_f, cov_f = curve_fit(Gauss,df['time'],df['f_data'])
    fit_asymp_f = param_f[0]
    fit_kobs_f = param_f[1]
    fit_y_f = Gauss(df['time'], fit_asymp_f, fit_kobs_f)

    param_g, cov_g = curve_fit(Gauss,df['time'],df['g_data'])
    fit_asymp_g = param_g[0]
    fit_kobs_g = param_g[1]
    fit_y_g = Gauss(df['time'], fit_asymp_g, fit_kobs_g)

    param_h, cov_h = curve_fit(Gauss,df['time'],df['h_data'])
    fit_asymp_h = param_h[0]
    fit_kobs_h = param_h[1]
    fit_y_h = Gauss(df['time'], fit_asymp_h, fit_kobs_h)

    well_kinetics = [fit_kobs_b/100,fit_kobs_c/100,fit_kobs_d/100,fit_kobs_e/100,fit_kobs_f/100,fit_kobs_g/100,fit_kobs_h/100]
    print(well_kinetics)
    
    fig, axs = plt.subplots(4, 2,figsize=(15,15))
    fig.tight_layout(pad=3.5)

    axs[0,0].set_title('{}mM POPC {}s, {}nM {}, fitted'.format(lip_conc,lip_type,conc[0],analyte))
    axs[0,0].plot(df['time'],df['b_data'],'#88CCEE',marker='o',label='{}nM {}'.format(conc[0],analyte))
    axs[0,0].plot(df['time'], fit_y_b,'black',label='fit, K = {}'.format((fit_kobs_b/100).round(5)))
    axs[0,0].legend(prop={'size':15},loc=position)
    axs[0,0].set(xlabel='Time (s)', ylabel='Wavelength Shift (nm)')

    axs[0,1].set_title('{}mM POPC {}s, {}nM {}, fitted'.format(lip_conc,lip_type,conc[1],analyte))
    axs[0,1].plot(df['time'],df['c_data'],'#44AA99',marker='o',label='{}nM {}'.format(conc[1],analyte))
    axs[0,1].plot(df['time'], fit_y_c,'black',label='fit, K = {}'.format((fit_kobs_c/100).round(5)))
    axs[0,1].legend(prop={'size':15},loc=position)
    axs[0,1].set(xlabel='Time (s)', ylabel='Wavelength Shift (nm)')

    axs[1,0].set_title('{}mM POPC {}s, {}nM {}, fitted'.format(lip_conc,lip_type,conc[2],analyte))
    axs[1,0].plot(df['time'],df['d_data'],'#117733',marker='o',label='{}nM {}'.format(conc[2],analyte))
    axs[1,0].plot(df['time'], fit_y_d,'black',label='fit, K = {}'.format((fit_kobs_d/100).round(5)))
    axs[1,0].legend(prop={'size':15},loc=position)
    axs[1,0].set(xlabel='Time (s)', ylabel='Wavelength Shift (nm)')

    axs[1,1].set_title('{}mM POPC {}s, {}nM {}, fitted'.format(lip_conc,lip_type,conc[3],analyte))
    axs[1,1].plot(df['time'],df['e_data'],'#999933',marker='o',label='{}nM {}'.format(conc[3],analyte))
    axs[1,1].plot(df['time'], fit_y_e,'black',label='fit, K = {}'.format((fit_kobs_e/100).round(5)))
    axs[1,1].legend(prop={'size':15},loc=position)
    axs[1,1].set(xlabel='Time (s)', ylabel='Wavelength Shift (nm)')

    axs[2,0].set_title('{}mM POPC {}s, {}nM {}, fitted'.format(lip_conc,lip_type,conc[4],analyte))
    axs[2,0].plot(df['time'],df['f_data'],'#DDCC77',marker='o',label='{}nM {}'.format(conc[4],analyte))
    axs[2,0].plot(df['time'], fit_y_f,'black',label='fit, K = {}'.format((fit_kobs_f/100).round(5)))
    axs[2,0].legend(prop={'size':15},loc=position)
    axs[2,0].set(xlabel='Time (s)', ylabel='Wavelength Shift (nm)')
    
    axs[2,1].set_title('{}mM POPC {}s, {}nM {}, fitted'.format(lip_conc,lip_type,conc[5],analyte))
    axs[2,1].plot(df['time'],df['g_data'],'#CC6677',marker='o',label='{}nM {}'.format(conc[5],analyte))
    axs[2,1].plot(df['time'], fit_y_g,'black',label='fit, K = {}'.format((fit_kobs_g/100).round(5)))
    axs[2,1].legend(prop={'size':15},loc=position)
    axs[2,1].set(xlabel='Time (s)', ylabel='Wavelength Shift (nm)')

    axs[3,0].set_title('{}mM POPC {}s, {}nM {}, fitted'.format(lip_conc,lip_type,conc[6],analyte))
    axs[3,0].plot(df['time'],df['h_data'],'#882256',marker='o',label='{}nM {}'.format(conc[6],analyte))
    axs[3,0].plot(df['time'], fit_y_h,'black',label='fit, K = {}'.format((fit_kobs_h/100).round(5)))
    axs[3,0].legend(prop={'size':15},loc=position)
    axs[3,0].set(xlabel='Time (s)', ylabel='Wavelength Shift (nm)')

    axs[3, 1].set_title('{}mM POPC {}s,K observed Plot '.format(lip_conc,lip_type))
    axs[3,1].plot(conc,well_kinetics,'o')
    axs[3,1].set(xlabel='PDBu []', ylabel='K observed')

    return well_kinetics
