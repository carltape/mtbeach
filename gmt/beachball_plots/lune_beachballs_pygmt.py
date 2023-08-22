import pygmt

import numpy as np

with pygmt.config(PS_MEDIA='9.5ix21.5i', PS_PAGE_ORIENTATION='landscape', PROJ_LENGTH_UNIT='inch', MAP_FRAME_TYPE='plain', MAP_TICK_LENGTH='0.0c', MAP_FRAME_PEN='2p', FONT_ANNOT_PRIMARY='12p,Helvetica,black', FONT_HEADING='18p,Helvetica,black', FONT_LABEL='10p,Helvetica,black'):
    fig = pygmt.Figure()

    # First Lune
    fig.basemap(projection='H0/2.8i', region='-30/30/-90/90', frame=['xa10f5g10','ya10f5g10','wesn+g200'], xshift='0.5', yshift='0.5')
    
    x,y = np.genfromtxt('../dfiles/sourcetype_patch_01.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, color='120', close='+x0')

    x,y = np.genfromtxt('../dfiles/sourcetype_patch_02.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, color='255', close='+x0')

    fig.basemap(projection='H0/2.8i', region='-30/30/-90/90', frame=['xa10f5g10','ya10f5g10','wesn'])

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_01.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,160/32/240')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_02.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,255/0/0')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_03.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_04.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_06.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    fig.meca(spec='../dfiles/beachballs_ipts1_iref1_lune_psmeca', scale='0.45i', convention='mt', no_clip=True, G='255/0/0', L='0.5p,0/0/0')

    
    # Title text
    fig.text(x=0, y=0, text='Reference sets of moment tensors (input files available in carltape compearth github repository)', font='16p,Helvetica-Bold,black', angle=0, justify='LM', projection='X1i', region='0/1/0/1', xshift='a0.5', yshift='a8.3', no_clip=True)
    fig.text(x=0, y=0, text='Plotted using pygmt', font='12p,Helvetica-Bold,black', angle=0, justify='LM', projection='X1i', region='0/1/0/1', xshift='a0.5', yshift='a8.1', no_clip=True)

    
    # Second Lune
    fig.basemap(projection='H0/2.8i', region='-30/30/-90/90', frame=['xa10f5g10','ya10f5g10','wesn+g200'], xshift='3.5')

    x,y = np.genfromtxt('../dfiles/sourcetype_patch_01.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, color='120', close='+x0')

    x,y = np.genfromtxt('../dfiles/sourcetype_patch_02.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, color='255', close='+x0')

    fig.basemap(projection='H0/2.8i', region='-30/30/-90/90', frame=['xa10f5g10','ya10f5g10','wesn'])

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_01.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,160/32/240')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_02.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,255/0/0')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_03.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_04.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_06.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    fig.meca(spec='../dfiles/beachballs_ipts1_iref2_lune_psmeca', scale='0.45i', convention='mt', no_clip=True, G='255/0/0', L='0.5p,0/0/0')


    # Third Lune
    fig.basemap(projection='H0/2.8i', region='-30/30/-90/90', frame=['xa10f5g10','ya10f5g10','wesn+g200'], xshift='3.5')
    
    x,y = np.genfromtxt('../dfiles/sourcetype_patch_01.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, color='120', close='+x0')

    x,y = np.genfromtxt('../dfiles/sourcetype_patch_02.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, color='255', close='+x0')

    fig.basemap(projection='H0/2.8i', region='-30/30/-90/90', frame=['xa10f5g10','ya10f5g10','wesn'])

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_01.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,160/32/240')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_02.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,255/0/0')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_03.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_04.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_06.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    fig.meca(spec='../dfiles/beachballs_ipts1_iref3_lune_psmeca', scale='0.45i', convention='mt', no_clip=True, G='255/0/0', L='0.5p,0/0/0')


    # Fourth Lune
    fig.basemap(projection='H0/2.8i', region='-30/30/-90/90', frame=['xa10f5g10','ya10f5g10','wesn+g200'], xshift='3.5')
    
    x,y = np.genfromtxt('../dfiles/sourcetype_patch_01.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, color='120', close='+x0')

    x,y = np.genfromtxt('../dfiles/sourcetype_patch_02.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, color='255', close='+x0')

    fig.basemap(projection='H0/2.8i', region='-30/30/-90/90', frame=['xa10f5g10','ya10f5g10','wesn'])

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_01.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,160/32/240')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_02.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,255/0/0')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_03.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_04.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_06.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    fig.meca(spec='../dfiles/beachballs_ipts1_iref4_lune_psmeca', scale='0.45i', convention='mt', no_clip=True, G='255/0/0', L='0.5p,0/0/0')

    # Fifth Lune
    fig.basemap(projection='H0/2.8i', region='-30/30/-90/90', frame=['xa10f5g10','ya10f5g10','wesn+g200'], xshift='3.5')
    
    x,y = np.genfromtxt('../dfiles/sourcetype_patch_01.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, color='120', close='+x0')

    x,y = np.genfromtxt('../dfiles/sourcetype_patch_02.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, color='255', close='+x0')

    fig.basemap(projection='H0/2.8i', region='-30/30/-90/90', frame=['xa10f5g10','ya10f5g10','wesn'])

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_01.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,160/32/240')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_02.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,255/0/0')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_03.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_04.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_06.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    fig.meca(spec='../dfiles/beachballs_ipts1_iref5_lune_psmeca', scale='0.45i', convention='mt', no_clip=True, G='255/0/0', L='0.5p,0/0/0')

    # Sixth Lune
    fig.basemap(projection='H0/2.8i', region='-30/30/-90/90', frame=['xa10f5g10','ya10f5g10','wesn+g200'], xshift='3.5')
    
    x,y = np.genfromtxt('../dfiles/sourcetype_patch_01.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, color='120', close='+x0')

    x,y = np.genfromtxt('../dfiles/sourcetype_patch_02.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, color='255', close='+x0')

    fig.basemap(projection='H0/2.8i', region='-30/30/-90/90', frame=['xa10f5g10','ya10f5g10','wesn'])

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_01.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,160/32/240')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_02.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,255/0/0')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_03.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_04.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    x,y = np.genfromtxt('../dfiles/sourcetype_arc_06.dat', unpack=True, usecols=(0,1))
    fig.plot(x=x, y=y, pen='3p,30/144/255')

    fig.meca(spec='../dfiles/beachballs_ipts2_iref3_lune_psmeca', scale='0.45i', convention='mt', no_clip=True, G='255/0/0', L='0.5p,0/0/0')

    fig.savefig('lune_hammer_iplot2_lplot1_kplot1_pygmt.pdf')
    



