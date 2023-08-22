import pygmt


with pygmt.config(PS_MEDIA='3.2ix7.7i', PS_PAGE_ORIENTATION='portrait', PROJ_LENGTH_UNIT='inch', MAP_FRAME_TYPE='plain', MAP_TICK_LENGTH='0.0c', MAP_FRAME_PEN='2p', MAP_TICK_PEN='2p'):
    fig = pygmt.Figure()
    fig.basemap(projection='X2.25i/6.72i',region='-30/30/-90/90',frame=['xa10f5g10','ya10f5g10','wesn+g200'],xshift='0.5',yshift='0.5')
    fig.meca(spec='beachballs_ipts1_iref1_lune_psmeca', scale='0.45i', convention='mt', no_clip=True, G='255/0/0', L='0.5p,0/0/0', W='0.5p,0/0/0')
    fig.savefig('lune_pygmt.pdf')
