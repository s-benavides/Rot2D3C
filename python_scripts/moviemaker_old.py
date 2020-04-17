import imageio

writer = imageio.get_writer('./movies/N4F1d375K22_ww_freeze.mp4', fps = 24)

qmax = 60*props['Q0']*np.max(Q)
Fmin = np.min(F)
Fmax = np.max(F)
umin = np.min(u)
umax = np.max(u)
wmin = np.min(w)
wmax = np.max(w)
cycles = 5

for t in np.linspace(0, cycles*2*np.pi/om, cycles*100):

     e = np.exp(-1j*om*t)
     fig = plt.figure(figsize = (25,8))

     plt.imshow(np.real(e)*60*props['Q0']*Q,
                 extent = [min(x)/1e3, max(x)/1e3, max(z)/1e3, 
min(z)/1e3], interpolation = 'none', cmap = 'RdBu_r')
     plt.colorbar(label = r'Prescribed buoyancy source $(m/s^2/min)$')
     plt.clim([-qmax, qmax])
     quiv = plt.quiver(
         xx[0::20,0::20]/1e3, zz[0::20,0::20]/1e3, 
np.real(e*u[0::20,0::20]), np.real(e*w[0::20,0::20]),
         scale = 50.)
     qk = plt.quiverkey(quiv, 0.65, 0.32, 1, '1 m/s', labelpos='W',
                    coordinates='figure')
     plt.gca().set_aspect(1)
     plt.ylim([0, 35])
     plt.xlabel('x (km)')
     plt.ylabel('y (km)')
     plt.title('Buoyancy source and wave motion')

     fig.canvas.draw()
     img = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
     img = img.reshape(fig.canvas.get_width_height()[::-1] + (3,))
     writer.append_data(img)
     plt.close()

writer.close()
