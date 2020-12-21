#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 11:33:24 2020

@author: rgryan
"""

import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob

## Path to ozonesondes
path = '/Volumes/GoogleDrive/My Drive/Documents/Postdoc/rockets/'
osp = path+'ozonesondes/ndacc_2019_ozonesondes/test/'
osp = path+'ozonesondes/ndacc_2019_ozonesondes/'
outpath = path+'ozonesondes/obsPack_input/'
plotAvgs = True
avgdfs = []

stations = ["dumont d'urville", "South Pole", "Belgrano", "Boulder", "Hilo",
             "Santa Cruz", "Natal Brazil" , "Neumayer", "Ny-Aalesund", 
             "Samoa", "Ittoqqortoormiit", "Wallops Flight Facility"]

def calc_DT(seconds):
        return startDT + pd.Timedelta(str(seconds)+' seconds')
    
def calc_improvement(s,wr,nr):
    wr_s, nr_s = wr-s, nr-s
    if abs(wr_s) > abs(nr_s): # things got worse!
        return -1*abs(wr-nr)
    else:
        return abs(wr-nr)
def improvement_cat(s,wr,nr):
    wr_s, nr_s = wr-s, nr-s
    if wr_s and nr_s < 0: # both simulations were less than sonde
        return 'blue'
    elif wr_s and nr_s > 0: # both simulations were greater than sonde
        return 'red'
    else:
        return 'orange'
    

## Using glob, this reads all .csv sonde files in the folder
files = glob.glob(osp+'*.*')
for osf in files[:]:
    print(osf[len(osp):])
    
    ## Use 'Readlines' to find the header of the file and location lat/lon
    f = open(osf, "r")
    c,flag,footer,l =0, 0, 0, [] # Counters for finding information in file
    for line in f:
        l.append(line)  
        c = c+1
        try:
            if list(filter(None, line.split(' ')))[0] == 's':
                header=c-1
                flag=1
            elif list(filter(None, line.split(' ')))[0] == 'hPa':
                header=c-1
                flag=1
        except IndexError:
            continue
        
        if line[:-1] in stations:
            #print( '    Station is '+line[:-1] )
            station, iGeo = line[:-1], c
    
    if flag>0:
        columns = list(filter(None, l[header-1].split(' ')))
        if columns[-1] == '\n':
            columns = columns[:-1]
        else:
            columns = columns
        #print(columns)
        
        top = list(filter(None, l[0].split(' ')))
        if top[-1] == '\n':
            top = top[:-1]
        else:
            top = top
        startDT = pd.to_datetime( top[-3:-2][0] + ' '+ top[-2:-1][0][:8],dayfirst=True) 
        endDT   = pd.to_datetime(top[-2:-1][0][8:]+' '+top[-1:][0][:8] ,dayfirst=True)
        gf = pd.read_csv(osf, delim_whitespace=True, header=header)
        gf.columns = columns
        
        df = pd.DataFrame()
        if 'ElapTime' in columns:
            df['Duration'] = gf['ElapTime']
        else:
            df['Duration'] = gf['Time']
        
        if 'Lat' and 'Lon' in columns:
            df['lat'], df['lon'] = gf['Lat'], gf['Lon']
        elif 'GPSLon' and 'GPSLat' in columns:
            df['lat'], df['lon'] = gf['GPSLat'], gf['GPSLon']
        else:
            geo = list(filter(None, l[iGeo].split(' ')))
            df['lat'], df['lon'] = float(geo[4]), float(geo[3])
            
        if 'Press' in columns:
            df['Pressure'] = gf['Press']
        else:
            print('different pressure label')
            
        if 'PO3' in columns:
            df['O3_sonde'] = 1e9*((gf['PO3']/1000)/(df['Pressure']*100))
        else:
            print('different o3 label')
    
      
        df['Date_Time'] = list(map(calc_DT, df['Duration']))
        df.set_index(pd.DatetimeIndex(df['Date_Time'], dayfirst=True), inplace=True)
        mth = str(df.index.month[0])
        
        #print('<><><><><> '+date+', '+station+', '+country+' <><><><><>')
        
        try:
            gcr1, gcr2, datestr = pd.DataFrame(), pd.DataFrame(), str(startDT)[:10].split('-')
            dateplus1 = str(startDT + pd.Timedelta('1 day'))[:10]
            nr = Dataset(path+'ozonesondes/obsPack_output/no_rockets/GEOSChem.ObsPack.'+datestr[0]+datestr[1]+datestr[2]+'_0000z.nc4')
            wr = Dataset(path+'ozonesondes/obsPack_output/with_rockets/GEOSChem.ObsPack.'+datestr[0]+datestr[1]+datestr[2]+'_0000z.nc4')

            gcr1['P_gchem'] = wr.variables['pressure'][:]
            gcr1['O3_gcnr'] = nr.variables['O3'][:]*1e9
            gcr1['O3_gcwr'] = wr.variables['O3'][:]*1e9
            gcr1['lat'], gcr1['lon'] = wr.variables['lat'][:], wr.variables['lon'][:]
            gcr1 = gcr1.astype('float')
            gcr1['strLat'] = [ str(round(gcr1['lat'][i],1)) for i in range(len(gcr1)) ]
        
            gcr1 = gcr1[gcr1['strLat'] == str(round(float(gcr1['lat'][0]),1)) ]
            gcr1 = gcr1[gcr1['P_gchem']>0]
        
            df.set_index(df['Pressure'], inplace=True)
            df = df[~df.index.duplicated(keep='first')]
        
            gcr1.set_index(gcr1['P_gchem'], inplace=True)
            gcr1 = gcr1[~gcr1.index.duplicated(keep='first')]
        
            ni = df['Pressure'].values.tolist() + gcr1['P_gchem'].values.tolist()
            ni = list( dict.fromkeys(ni) )
            df_ = df.reindex(ni).sort_index().interpolate()
            #df_ = df_.reindex(gcr1['P_gchem'].values.tolist())
            gcr1['sonde'] = df_['O3_sonde']
            gcr1['wrDiff'] = gcr1['O3_gcwr'] - gcr1['sonde']
            gcr1['nrDiff'] = gcr1['O3_gcnr'] - gcr1['sonde']
            gcr1['imp'] = list(map(calc_improvement, gcr1['sonde'],gcr1['O3_gcwr'],gcr1['O3_gcnr']))
            gcr1['cat'] = list(map(improvement_cat, gcr1['sonde'],gcr1['O3_gcwr'],gcr1['O3_gcnr']))
            gcr1['zero'] = 0
            imp, son = gcr1['imp'], gcr1['sonde']
        gcw, gcn = gcr1['O3_gcwr'], gcr1['O3_gcnr']
        
        if plotAvgs == False:
            fig = plt.figure(figsize=(4,3))
            h=fig.add_subplot(111)
            #h.plot(gcr1['sonde'].values, gcr1['P_gchem'].values, color='orange', zorder=0)
            #h.plot(gcr1['O3_gcnr'].values, gcr1['P_gchem'].values,  '--', color='darkred', zorder=1)
            #h.plot(gcr1['O3_gcwr'].values, gcr1['P_gchem'].values, color='dodgerblue', zorder=2)
        
            #h.scatter(gcr1['wrDiff'].values, gcr1['P_gchem'].values, c='blue',marker='o')
            #h.scatter(gcr1['nrDiff'].values, gcr1['P_gchem'].values, c='red', marker='x')
            h.scatter(gcr1['imp'].values, gcr1['P_gchem'].values, color=gcr1['cat'].values, marker='x')
            h.plot(gcr1['zero'].values, gcr1['P_gchem'].values, '--', color='black' )
            
        try:
            wrT =Dataset(path+'ozonesondes/obsPack_output/with_rockets/GEOSChem.ObsPack.'+datestrplus1+'_0000z.nc4')
            nrT = Dataset(path+'ozonesondes/obsPack_output/no_rockets/GEOSChem.ObsPack.'+datestrplus1+'_0000z.nc4')
        
            gcr2['P_gchem'] = wrT.variables['pressure'][:]
            gcr2['O3_gcnr'] = nrT.variables['O3'][:]*1e9
            gcr2['O3_gcwr'] = wrT.variables['O3'][:]*1e9
            gcr2['lat'], gcr2['lon'] = wrT.variables['lat'][:], wrT.variables['lon'][:]
            gcr2 = gcr2.astype('float')
            gcr2['strLat'] = [ str(round(gcr2['lat'][i],1)) for i in range(len(gcr2)) ]
            gcr2 = gcr2[gcr2['strLat'] == str(round(float(lat),1)) ]
            gcr2 = gcr2[gcr2['P_gchem']>0]
            gcr2.set_index(gcr2['P_gchem'], inplace=True)
            gcr2 = gcr2[~gcr2.index.duplicated(keep='first')]
        
            ni = df['Pressure'].values.tolist() + gcr2['P_gchem'].values.tolist()
            ni = list( dict.fromkeys(ni) )
            df_ = df.reindex(ni).sort_index().interpolate()
            #df_ = df_.reindex(gcr2['P_gchem'].values.tolist())
            gcr2['sonde'] = df_['O3_sonde']
            gcr2['wrDiff'] = gcr2['O3_gcwr'] - gcr2['sonde']
            gcr2['nrDiff'] = gcr2['O3_gcnr'] - gcr2['sonde']
            
            gcr2['imp'] = list(map(calc_improvement, gcr2['sonde'],gcr2['O3_gcwr'],gcr2['O3_gcnr']))
            gcr2['cat'] = list(map(improvement_cat, gcr2['sonde'],gcr2['O3_gcwr'],gcr2['O3_gcnr']))
            gcr2['zero'] = 0
            imp, son = pd.concat([gcr1['imp'], gcr2['imp']]), pd.concat([gcr1['sonde'], gcr2['sonde']])
            gcw, gcn = pd.concat([gcr1['O3_gcwr'], gcr2['O3_gcwr']]), pd.concat([gcr1['O3_gcnr'], gcr2['O3_gcnr']])
            
            if plotAvgs == False:
                #h.plot(gcr2['sonde'].values, gcr2['P_gchem'].values, color='orange', zorder=0)
                #h.plot(gcr2['O3_gcwr'].values, gcr2['P_gchem'].values, color='dodgerblue', zorder=2)
                #h.plot(gcr2['O3_gcnr'].values, gcr2['P_gchem'].values, '--', color='darkred', zorder=2)
            
                #h.scatter(gcr2['wrDiff'].values, gcr2['P_gchem'].values, c='blue', marker='o')
                #h.scatter(gcr2['nrDiff'].values, gcr2['P_gchem'].values, c='red', marker='x')
            
                h.scatter(gcr2['imp'].values, gcr2['P_gchem'].values, color=gcr2['cat'].values, marker='x')
                h.plot(gcr2['zero'].values, gcr2['P_gchem'].values, '--', color='black')

            
        except FileNotFoundError:
            print(' No tomorrow data ')
        
        if plotAvgs == False:
            h.set_ylim(5,500)
            h.set_yscale('log')
            h.set_ylabel('Pressure (hPa)')
            h.set_xlabel('O$_3$ diff (ppb)')
            plt.gca().invert_yaxis()
            #h.legend(['Rockets-Sonde', 'No Rockets-Sonde'], loc='upper right')
            h.set_title(station+', '+country+' '+date)
            fig.savefig(path+'/ozonesondes/plots/imp_'+station+'_'+country+'_'+date+'.png', bbox_inches='tight')
        
        if plotAvgs==True:
            # get space-free station code:
            sc = station.replace(' ','')
            sc = sc.replace( ")","" )
            sc = sc.replace( "(","" )
            try:
                exec(sc+'_imp_'+mth+'['+datestr+'] = imp')
                exec(sc+'_son_'+mth+'['+datestr+'] = son')
                exec(sc+'_gcw_'+mth+'['+datestr+'] = gcw')
                exec(sc+'_gcn_'+mth+'['+datestr+'] = gcn')
            except NameError:
                exec(sc+'_imp_'+mth+'=pd.DataFrame()')
                exec(sc+'_son_'+mth+'=pd.DataFrame()')
                exec(sc+'_gcw_'+mth+'=pd.DataFrame()')
                exec(sc+'_gcn_'+mth+'=pd.DataFrame()')
                
                avgdfs_imp.append(sc+'_imp_'+mth)
                avgdfs_son.append(sc+'_son_'+mth)
                avgdfs_gcn.append(sc+'_gcn_'+mth)
                avgdfs_gcw.append(sc+'_gcw_'+mth)
                
                exec(sc+'_imp_'+mth+'['+datestr+'] = imp')
                exec(sc+'_son_'+mth+'['+datestr+'] = son')
                exec(sc+'_gcn_'+mth+'['+datestr+'] = gcn')
                exec(sc+'_gcw_'+mth+'['+datestr+'] = gcw')

    except (KeyError, TypeError, ValueError) as error:
        print(' -------X--->  Data incomplete/cannot be processed at '+station+', '+country)

#%%    
if plotAvgs == True:
    if len(avgdfs_imp)>0:
        for d in range(len(avgdfs_imp)):
            exec('n=len('+avgdfs_imp[d]+'.columns)')
            exec(avgdfs_imp[d]+'["mean"] = '+avgdfs_imp[d]+'.mean(axis=1)')
            exec(avgdfs_imp[d]+'["std"] = '+avgdfs_imp[d]+'.std(axis=1)')
            
            exec(avgdfs_son[d]+'["mean"] = '+avgdfs_son[d]+'.mean(axis=1)')
            exec(avgdfs_son[d]+'["std"] = '+avgdfs_son[d]+'.std(axis=1)')
            
            exec(avgdfs_gcn[d]+'["mean"] = '+avgdfs_gcn[d]+'.mean(axis=1)')
            exec(avgdfs_gcn[d]+'["std"] = '+avgdfs_gcn[d]+'.std(axis=1)')
            
            exec(avgdfs_gcw[d]+'["mean"] = '+avgdfs_gcw[d]+'.mean(axis=1)')
            exec(avgdfs_gcw[d]+'["std"] = '+avgdfs_gcw[d]+'.std(axis=1)')
            
            #fjg = plt.figure(figsize=(4,3))
            fig, j = plt.subplots(1, 2, sharey=True, figsize=(6,4))
            
            exec('j[0].errorbar('+avgdfs_son[d]+'["mean"].values, '+avgdfs_son[d]+'.index.values, xerr='+avgdfs_son[d]+'["std"].values ,\
                             fmt = "o--", color="orange", mfc="yellow")')
            exec('j[0].errorbar('+avgdfs_gcn[d]+'["mean"].values, '+avgdfs_gcn[d]+'.index.values, xerr='+avgdfs_gcn[d]+'["std"].values ,\
                             fmt = "s-", color="red", mfc="white")')
            exec('j[0].errorbar('+avgdfs_gcw[d]+'["mean"].values, '+avgdfs_gcw[d]+'.index.values, xerr='+avgdfs_gcw[d]+'["std"].values ,\
                             fmt = "d-", color="blue", mfc="cyan")')
            j[0].set_ylim(5,500)
            j[0].set_yscale('log')
            j[0].set_xlabel('O$_3$ profile (ppb)')
            #plt.gca().invert_yaxis()
            j[0].set_title(avgdfs_son[d]+', n='+str(n))
            j[0].legend(['Sonde', 'GC_NR', 'GC_WR'], loc='lower right')
            
            exec('j[1].errorbar('+avgdfs_imp[d]+'["mean"].values, '+avgdfs_imp[d]+'.index.values, xerr='+avgdfs_imp[d]+'["std"].values ,\
                             fmt = "o", color="teal", mfc="white")')
            j[1].set_ylim(5,500)
            j[1].set_yscale('log')
            j[1].set_ylabel('Pressure (hPa)')
            j[1].set_xlabel('O$_3$ improvement (ppb)')
            plt.gca().invert_yaxis()
            j[1].set_title(avgdfs_imp[d]+', n='+str(n))
                        
            
            fig.savefig(path+'/ozonesondes/plots/imp_Avgs_andProfiles_'+avgdfs_imp[d]+'.png', bbox_inches='tight')