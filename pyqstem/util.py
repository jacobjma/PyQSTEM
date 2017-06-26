import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import warnings
from ase.data.colors import cpk_colors
from ase.data import covalent_radii,chemical_symbols
from matplotlib.patches import Circle
from scipy.cluster.hierarchy import linkage,fcluster
from ase import Atoms
import collections

def energy2wavelength(v0):
    m0=0.5109989461*10**3 # keV / c**2
    h=4.135667662*10**(-15)*10**(-3) # eV * s
    c=2.99792458*10**8 # m / s
    return h*c/np.sqrt(v0*(2*m0+v0))*10**10

def focal_spread(Cc,energy,energy_spread,current_stability):
    return Cc*np.sqrt((dv0/v0)**2+2*current_stability**2)

def scherzer(v0,Cs):
    assert Cs > 0.
    return -1.2*np.sqrt(wavelength(v0)*Cs)

def scherzer_point_resolution(v0,Cs):
    assert Cs > 0.
    return 0.6*wavelength(v0)**(3/4.)*Cs**(1/4.)

def spatial_frequencies(shape,sampling,return_polar=False,return_nyquist=False,wavelength=None):

    if not isinstance(shape, collections.Iterable):
        shape = (shape,)*2

    if not isinstance(sampling, collections.Iterable):
        sampling = (sampling,)*2
    
    dkx=1/(shape[0]*sampling[0])
    dky=1/(shape[1]*sampling[1])

    if shape[0]%2==0:
        kx = np.fft.fftshift(dkx*np.arange(-shape[0]/2,shape[0]/2,1))
    else:
        kx = np.fft.fftshift(dkx*np.arange(-shape[0]/2-.5,shape[0]/2-.5,1))

    if shape[1]%2==0:
        ky = np.fft.fftshift(dky*np.arange(-shape[1]/2,shape[1]/2,1))
    else:
        ky = np.fft.fftshift(dky*np.arange(-shape[1]/2-.5,shape[1]/2-.5,1))

    ky,kx = np.meshgrid(ky,kx)

    k2 = kx**2+ky**2

    ret = (kx,ky,k2)
    if return_nyquist:
        knx = 1/(2*sampling[0])
        kny = 1/(2*sampling[1])
        Kx=kx/knx
        Ky=ky/knx
        K2 = Kx**2+Ky**2
        ret += (Kx,Ky,K2)
    if return_polar:
        theta = np.sqrt(k2*wavelength**2)
        phi = np.arctan2(ky,kx)
        ret += (theta,phi)

    return ret

def atoms_plot(atoms,direction=2,ax=None,scan_range=None,potential_extent=None,probe_extent=None,s=40,boundary=1,scale_atoms=1,legend=False):

    if not np.allclose(atoms.get_cell(), np.diag(np.diag(atoms.get_cell()))):
        warnings.warn("Ignoring non-diagonal components of unit cell")

    if direction == 'x':
        direction = 0
    elif direction == 'y':
        direction = 1
    elif direction == 'z':
        direction = 2

    axes=np.delete([0,1,2],direction)
    labels=np.delete(['x','y','z'],direction)

    direction_sort=np.argsort(-atoms.get_positions()[:,direction])
    atoms=atoms[direction_sort]

    positions=atoms.get_positions()[:,axes]

    cell=atoms.get_cell()[axes,axes]

    if ax is None:
        fig, ax = plt.subplots()

    handles={}
    for atom,pos in zip(atoms,positions):
        handles[atom.symbol]=ax.add_artist(Circle(xy=(pos[0], pos[1]), facecolor=cpk_colors[atom.number], radius=scale_atoms*covalent_radii[atom.number],
                        zorder=0,label=atom.symbol,lw=1,edgecolor='k'))

    handles = list(handles.values())

    if potential_extent is not None:
        handles.append(ax.add_patch(patches.Rectangle((potential_extent[0],potential_extent[2]),
                potential_extent[1]-potential_extent[0],potential_extent[3]-potential_extent[2],
                alpha=0.2,color='g',linewidth=None,zorder=1,label='potential extent')))
        ax.add_patch(patches.Rectangle((potential_extent[0],potential_extent[2]),
                potential_extent[1]-potential_extent[0],potential_extent[3]-potential_extent[2],
                linestyle='dashed',linewidth=1,fill=False,zorder=1))

    if probe_extent is not None:
        handles.append(ax.add_patch(patches.Rectangle((probe_extent[0],probe_extent[2]),
                probe_extent[1]-probe_extent[0],probe_extent[3]-probe_extent[2],
                alpha=0.2,color='b',linewidth=None,zorder=1,label='probe max extent')))
        ax.add_patch(patches.Rectangle((probe_extent[0],probe_extent[2]),
                probe_extent[1]-probe_extent[0],probe_extent[3]-probe_extent[2],
                linestyle='dashed',linewidth=1,fill=False,zorder=1))

    if scan_range is not None:
        handles.append(ax.add_patch(patches.Rectangle((scan_range[0][0],scan_range[1][0]),
                scan_range[0][1]-scan_range[0][0],scan_range[1][1]-scan_range[1][0],
                alpha=0.2,color='r',linewidth=None,zorder=3,label='scan range')))
        ax.add_patch(patches.Rectangle((scan_range[0][0],scan_range[1][0]),
                scan_range[0][1]-scan_range[0][0],scan_range[1][1]-scan_range[1][0],
                linestyle='dashed',linewidth=1,fill=False,zorder=3))
        dx=(scan_range[0][1]-scan_range[0][0])/float(scan_range[0][2])
        dy=(scan_range[1][1]-scan_range[1][0])/float(scan_range[1][2])
        for i in range(1,scan_range[0][2]):
            ax.plot([scan_range[0][0]+i*dx]*2,[scan_range[1][0],scan_range[1][1]],'k-',alpha=.4,zorder=2)
        for i in range(1,scan_range[1][2]):
            ax.plot([scan_range[0][0],scan_range[0][1]],[scan_range[1][0]+i*dy]*2,'k-',alpha=.4,zorder=2)

    if legend:
        ax.legend(handles=handles)

    ax.plot([0,0,cell[0],cell[0],0],[0,cell[1],cell[1],0,0],'k',linewidth=1.5)
    ax.axis('equal')
    ax.set_xlim([-boundary,cell[0]+boundary])
    ax.set_ylim([-boundary,cell[1]+boundary])
    ax.set_xlabel('{0} [Angstrom]'.format(labels[0]))
    ax.set_ylabel('{0} [Angstrom]'.format(labels[1]))

def project_positions(atoms,distance=1,return_counts=False):
    
    positions=atoms.get_positions()[:,:2]

    clusters = fcluster(linkage(positions), distance, criterion='distance')
    unique, indices = np.unique(clusters, return_index=True)
    positions = np.array([np.mean(positions[clusters == u], axis=0) for u in unique])
    
    if return_counts:
        counts=np.array([np.sum(clusters == u) for u in unique])
        return positions, counts
    else:
        return positions

def draw_scalebar(pil_img,scale_length,sampling,units='nm',placement=[5,5],margin=3,bar_height=3,
                  font=None,bar_color=0,bg_color=None,formatting='1f',anchor='top left'):

    from PIL import ImageDraw

    placement=list(placement)

    draw = ImageDraw.Draw(pil_img)
    bar_length=scale_length/sampling

    text='{0:.{formatting}} {1}'.format(scale_length,units,formatting=formatting)
    text_size=draw.textsize(text, font=font)
    text_spacing=0

    bg_height=text_spacing+bar_height+text_size[1]

    if anchor=='top right':
        placement[0]-=margin*2+bar_length
    elif anchor=='bottom left':
        placement[1]-=margin*2+bg_height
    elif anchor=='bottom right':
        placement[0]-=margin*2+bar_length
        placement[1]-=margin*2+bg_height

    if bg_color is not None:
        draw.rectangle([(placement[0],placement[1]),
                (placement[0]+bar_length+2*margin,placement[1]+bg_height+2*margin)], fill=bg_color)

    draw.rectangle([(placement[0]+margin,placement[1]+bg_height-bar_height+margin),
                (placement[0]+bar_length+margin,placement[1]+bg_height+margin)], fill=bar_color)

    draw.text((placement[0]+margin+bar_length//2-text_size[0]//2,placement[1]+margin),text,bar_color,font=font)

    return pil_img
