import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors
import pickle
import os

import scipy.stats
import scipy.interpolate

plt.close('all')

z0     = np.asarray(
    [np.complex(np.cos(-0.25*np.pi),np.sin(-0.25*np.pi)),
     np.complex(np.cos(0.25*np.pi),np.sin(0.25*np.pi)),
     np.complex(np.cos(0.75*np.pi),np.sin(0.75*np.pi))])

def angle_to_unit_circle(r):
    
    import numpy as np
    
    # Angle must be provided in radians, counter-clockwise from 3 o'clock
    return np.cos(r)+1j*np.sin(r)

def find_moebius_coefficients(r,z0):
    
    import numpy as np
    
    # Find the images of the z0 control points
    w0  = angle_to_unit_circle(r)
    
    # Then calculate the four parameters for the corresponding Möbius map
    a = np.linalg.det(np.asarray(
        [[z0[0]*w0[0],     w0[0],          1],
         [z0[1]*w0[1],     w0[1],          1],
         [z0[2]*w0[2],     w0[2],          1]]))
    
    b = np.linalg.det(np.asarray(
        [[z0[0]*w0[0],     z0[0],     w0[0]],
         [z0[1]*w0[1],     z0[1],     w0[1]],
         [z0[2]*w0[2],     z0[2],     w0[2]]]))
    
    c = np.linalg.det(np.asarray(
        [[z0[0],           w0[0],          1],
         [z0[1],           w0[1],          1],
         [z0[2],           w0[2],          1]]))
    
    d = np.linalg.det(np.asarray(
        [[z0[0]*w0[0],     z0[0],     1],
         [z0[1]*w0[1],     z0[1],     1],
         [z0[2]*w0[2],     z0[2],     1]]))
    
    return a,b,c,d

def moebius(z,a,b,c,d,inverse=False):
    
    if not inverse:
        z = (a*z+b)/(c*z+d)
    else:
        z = (-d*z+b)/(c*z-a)
    
    return z

# =============================================================================
# Check for all possible rotations
resolution  = 360
rots        = np.linspace(-np.pi,np.pi,resolution+1)[:-1]

# Create a matrix checking for inversion
inversion_matrix    = np.ones((resolution,resolution,resolution))

x,y = np.meshgrid(
    np.linspace(-1,1,101),
    np.linspace(-1,1,101))

xy = np.ndarray.flatten(x + 1j*y)
xy  = xy[np.where(np.abs(xy) <= 1)]

"""

counter1    = 0
counter2    = 0

for idx1 in range(resolution):
    print(idx1)
    for idx2 in range(resolution):
        for idx3 in range(resolution):
            
            a,b,c,d = find_moebius_coefficients(
                r   = [rots[idx1],rots[idx2],rots[idx3]],
                z0  = z0)
            
            if abs(b/d) < 1:
                inversion_matrix[idx1,idx2,idx3]= 0
                
            # if np.random.uniform() < 0.0001 and counter1 < 10 and inversion_matrix[idx1,idx2,idx3] == 0:
            #     z = moebius(xy,a,b,c,d,inverse=False)
            #     counter1 += 1
            #     plt.figure()
            #     plt.title('no inversion')
            #     plt.tricontour(np.real(xy),np.imag(xy),np.real(z))
            #     plt.tricontour(np.real(xy),np.imag(xy),np.imag(z))
            #     plt.colorbar()

            # if np.random.uniform() < 0.0001 and counter2 < 10 and inversion_matrix[idx1,idx2,idx3] == 1:
            #     z = moebius(xy,a,b,c,d,inverse=False)
            #     counter2 += 1
            #     plt.figure()
            #     plt.title('inversion')
            #     plt.tricontour(np.real(xy),np.imag(xy),np.real(z))
            #     plt.tricontour(np.real(xy),np.imag(xy),np.imag(z))
            #     plt.colorbar()

            # if counter1 == 10 and counter2 == 10:
            #     raise Exception

pickle.dump(inversion_matrix,open('moebius_support_inversion_matrix.p','wb'))

"""

inversion_matrix    = pickle.load(open('moebius_support_inversion_matrix.p','rb'))




# fig = plt.figure()
# ax = fig.gca(projection='3d')

# cmap = plt.get_cmap("viridis")
# norm= plt.Normalize(0, 1)
# ax.voxels(inversion_matrix[:50,:50,:50], facecolors=cmap(norm(inversion_matrix[:50,:50,:50])), edgecolor=None, alpha = 0.25)

# plt.show()

# raise Exception



inversion_matrix[np.where(inversion_matrix > 0)] = np.nan

# Then we prepare a new colormap, mainly for cosmetic purposes.
# cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", 
#     ["xkcd:grey",
#      "xkcd:cerulean"])
# cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", 
#     ["xkcd:grey",
#      "xkcd:kermit green"])
# cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", 
#     ["xkcd:grey",
#      "xkcd:orangish red"])
cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", 
    ["xkcd:cerulean",
     "xkcd:cerulean"])
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", 
    ["xkcd:kermit green",
     "xkcd:kermit green"])
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", 
    ["xkcd:orangish red",
     "xkcd:orangish red"])

root_directory = os.path.dirname(os.path.realpath(__file__))

if not os.path.exists(root_directory+'\\'+'moebius support figure'):
    os.makedirs(root_directory+'\\'+'moebius support figure')
    
X,Y     = np.meshgrid(
    rots,rots)


os.chdir(root_directory + '\\' + 'moebius support figure')

for idx in range(resolution):
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.contourf(
        inversion_matrix[idx,:,:],X,Y,
        levels  = 2,
        alpha   = 0.5,
        zdir    = 'x',
        offset  = rots[idx],
        cmap    = cmap1)
    plt.contourf(
        X,inversion_matrix[:,idx,:],Y,
        levels  = 2,
        alpha   = 0.5,
        zdir    = 'y',
        offset  = rots[idx],
        cmap    = cmap2)
    plt.contourf(
        X,Y,inversion_matrix[:,:,idx],
        levels  = 2,
        alpha   = 0.5,
        zdir    = 'z',
        offset  = rots[idx],
        cmap    = cmap3)
    
    plt.gca().set_xlabel('rotation for control point $A$')
    plt.gca().set_ylabel('rotation for control point $B$')
    plt.gca().set_zlabel('rotation for control point $C$')
    
    plt.gca().set_xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
    plt.gca().set_xticklabels(['$-\pi$','$-\pi/2$','0','$\pi/2$','$\pi$'])
    plt.gca().set_yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
    plt.gca().set_yticklabels(['$-\pi$','$-\pi/2$','0','$\pi/2$','$\pi$'])
    plt.gca().set_zticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
    plt.gca().set_zticklabels(['$-\pi$','$-\pi/2$','0','$\pi/2$','$\pi$'])
    
    plt.savefig('img_'+str(idx).zfill(3)+'.png')
    plt.close('all')

