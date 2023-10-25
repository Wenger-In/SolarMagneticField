from cProfile import label
import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
from matplotlib.lines import Line2D
from matplotlib.colors import SymLogNorm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sunpy.map
import pfsspy
from pfsspy import coords
from pfsspy import tracing

flap_case = 3

###############################################################################
# Fuctions
def set_axes_lims(ax):
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 180)
    return 0
###############################################################################
# Step 1: import magnetogram
if flap_case == 1:
    gong_fname = 'E:/Research/Data/GONG/20200608.fits/202006080114.fits'
elif flap_case == 2:
    gong_fname = 'E:/Research/Data/GONG/20200921.fits/202009210114.fits'
elif flap_case == 3:
    gong_fname = 'E:/Research/Data/GONG/20210117.fits/202101171314.fits'

gong_map = sunpy.map.Map(gong_fname)
# Remove the mean, sothat curl B = 0; set colorbar to be symlog
norm = SymLogNorm(linthresh=5)
gong_map = sunpy.map.Map(gong_map.data - np.mean(gong_map.data), gong_map.meta, plot_settings={'norm': norm})
# Change limit of magnetogram colorbar
lim = 20
for i in range(len(gong_map.data)):
    for j in range(len(gong_map.data[i])):
        if gong_map.data[i][j] > lim:
            gong_map.data[i][j] = lim
        elif gong_map.data[i][j] < -lim:
            gong_map.data[i][j] = -lim
###############################################################################
# Step 2: set grids and calculate
nrho = 30
rss = 2.5  # unit: solar radii
input = pfsspy.Input(gong_map, nrho, rss)
output = pfsspy.pfss(input)
###############################################################################
# Step 3: plot input GONG magnetogram (Figure 1)
m = input.map
# Create the figure and axes
fig = plt.figure()
ax = plt.subplot(projection=m)
# Plot the GONG magnetogram
m.plot()
# Plot formatting
plt.colorbar()
ax.set_title('Input magnetogram field')
set_axes_lims(ax)
###############################################################################
# Step 4: plot output source surface map (Figure 2)
ss_br = output.source_surface_br
# Create the figure and axes
fig = plt.figure()
ax = plt.subplot(projection=ss_br, label='Neutral Line')
# Plot the source surface map
ss_br.plot()

# Plot the polarity inversion line
ax.plot_coord(output.source_surface_pils[0])
# Plot PSP trace
PSP_orb = [(67.923843, 90 - 2.4845026), (97.512878, 90 - 3.8445756)]
(PSP_orb_beg, PSP_orb_end) = zip(*PSP_orb)
ax.add_line(Line2D(PSP_orb_beg, PSP_orb_end, linewidth=2, color='black', label='Trajectory of PSP'))
# PSP_bfr = [(67.923843, 90 - 2.4845026), (83.312325, 90 - 3.4472520)]
# (PSP_bfr_beg, PSP_bfr_end) = zip(*PSP_bfr)
# ax.add_line(Line2D(PSP_bfr_beg, PSP_bfr_end, linewidth=2, color='black'))
# PSP_case = [(83.312325, 90 - 3.4472520), (85.887192, 90 - 3.5546808)]
# (PSP_beg, PSP_end) = zip(*PSP_case)
# ax.add_line(Line2D(PSP_beg, PSP_end, linewidth=2, color='black'))
# PSP_ltr = [(85.887192, 90 - 3.5546808), (97.512878, 90 - 3.8445756)]
# (PSP_ltr_beg, PSP_ltr_end) = zip(*PSP_ltr)
# ax.add_line(Line2D(PSP_ltr_beg, PSP_ltr_end, linewidth=2, color='black'))
# Plot formatting
plt.colorbar(orientation='horizontal')
ax.set_title('Magnetic field on source surface @ 2.5 Rs')
set_axes_lims(ax)

plt.legend(loc='upper right', frameon=True)

plt.show()
###############################################################################
# Step 5: plot PFSS solution in meridian plane (Figure 3)
fig, ax = plt.subplots()
ax.set_aspect('equal')

# Take 100 start points spaced equally in theta
r = 1.01 * const.R_sun
lon = np.pi / 2 * u.rad
lat = np.linspace(-np.pi / 2, np.pi / 2, 20) * u.rad
seeds = SkyCoord(lon, lat, r, frame=output.coordinate_frame)

tracer = pfsspy.tracing.PythonTracer()
field_lines = tracer.trace(seeds, output)

for field_line in field_lines:
    coords = field_line.coords
    coords.representation_type = 'cartesian'
    color = {0: 'black', -1: 'tab:red', 1: 'tab:blue'}.get(field_line.polarity)
    ax.plot(coords.y / const.R_sun,
            coords.z / const.R_sun, color=color)

# Add inner and outer boundary circles
ax.add_patch(mpatch.Circle((0, 0), 1, color='k', fill=False))
ax.add_patch(mpatch.Circle((0, 0), input.grid.rss, color='k', linestyle='--',
                           fill=False))
ax.set_title('PFSS solution')
###############################################################################
# Step 6: plot PFSS solution in 3D (Figure 4)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
###############################################################################
r = 1.2 * const.R_sun
lat = np.linspace(-np.pi / 2, np.pi / 2, 10, endpoint=False)
lon = np.linspace(0, 2 * np.pi, 10, endpoint=False)
lat, lon = np.meshgrid(lat, lon, indexing='ij')
lat, lon = lat.ravel() * u.rad, lon.ravel() * u.rad

seeds = SkyCoord(lon, lat, r, frame=output.coordinate_frame)

field_lines = tracer.trace(seeds, output)

for field_line in field_lines:
    color = {0: 'black', -1: 'tab:red', 1: 'tab:blue'}.get(field_line.polarity)
    coords = field_line.coords
    coords.representation_type = 'cartesian'
    ax.plot(coords.x / const.R_sun,
            coords.y / const.R_sun,
            coords.z / const.R_sun,
            color=color, linewidth=1)

ax.set_title('PFSS solution')
plt.show()
