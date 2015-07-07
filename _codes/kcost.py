# Import NUMPY and MATPLOTLIB libraries
    import polycentroid as centroid
    import numpy as np
    import scipy.optimize as opt
    import matplotlib.pyplot as plt
    from matplotlib.path import Path
    import matplotlib.patches as patches
    from gkwmod import gkwmod as gkw
    from mpl_toolkits.axes_grid.axislines import Subplot
#from matplotlib2tikz import save as tikz_save


def kcost():
    
    """
    Created on Tue Apr 27 14:18:29 2014

    @author: skywalker (tcy.kam@gmail.com)
    """



# Parameters
    siggma = sigma # Matching probability
    betta = beta  # Discount factor
    thetta = theta    # CRRA parameter
    alfa = alpha      # Convexity: cost function 
    A = 1.0         # Scaling: cost function
    qbar = 0.001  # shift: domain for u(q), so we have u(q + qbar)

    Inflation = Pi    # Money supply growth (Home)
    InflationF = PiF      # Money supply growth (Foreign)

    qmin = 0.01
    qmax = 30.0

# Plot font sizes
    ftsmall = 10
    ftmed = 11
    ftlarge = 12

# Instantiate Model Class

    model = gkw( siggma, betta, thetta, alfa, A, qbar )


#
# Case 3: 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1. Solve for qhat (SME allocation of q) by Brent's algorithm

    def excessmargin_sme(q, model, Inflation):
        """Gap between MU and MC"""
        wedge = ( Inflation - model.betta*(1-model.siggma) ) \
                                    / (model.siggma * model.betta)
        cprime = model.costdiff(q)
        uprime = model.utildiff(q)
        xcess = wedge * cprime - uprime
        return xcess

    def excessmargin_firstbest(q, model, Inflation):
        """Gap between MU and MC"""
        cprime = model.costdiff(q)
        uprime = model.utildiff(q)
        xcess = cprime - uprime
        return xcess

# Case 3's SME q traded
    qhat = opt.brentq(excessmargin_sme, qmin, qmax,\
                        args=(model, Inflation), xtol=1e-12,    \
                        rtol=4.4408920985006262e-16, maxiter=100, \
                        full_output=True, disp=True)

# First best q traded
    qstar = opt.brentq(excessmargin_firstbest, qmin, qmax, \
                        args=(model, Inflation), xtol=1e-12,    \
                        rtol=4.4408920985006262e-16, maxiter=100, \
                        full_output=True, disp=True)

# Case 3: 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2. Vertices of polygons from Case 3 linear inequalities

    kmax= 2.0  # arbitrary upper bound for plotting figure box
    kfmax = kmax*0.25

    k_threshold = model.cost(qhat[0])*(Inflation - betta*(1-siggma)) 
    kf_threshold = model.cost(qhat[0])*(InflationF - betta*(1-siggma)) 

    v1 = (kmax, 0.)
    v2 = (kmax, kf_threshold)
    v3 = (0., kf_threshold)
    v4 = (k_threshold, 0.)

    verts3 = [ v1,   # R-B
               v2,   # R-T
               v3,   # L-T
               v4,   # L-B
               v1,    # complete simple walk: ignored in patch
            ]

    graph3 = [ Path.MOVETO,
               Path.LINETO,
               Path.LINETO,
               Path.LINETO,
               Path.CLOSEPOLY,
            ]

    path3 = Path(verts3, graph3)


    c3 = centroid.calculate_polygon_centroid(verts3)


# Case 4:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    origin = (0., 0.)

    verts4 = [  origin,
                v3,
                v4,
                origin,
             ]

    graph4 = [ Path.MOVETO,
               Path.LINETO,
               Path.LINETO,
               Path.CLOSEPOLY,
            ]

    path4 = Path(verts4, graph4)


    c4 = centroid.calculate_polygon_centroid(verts4)


# Case 1 and 2 are obvious
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    verts12 = [ v3,
                (0., kfmax),
                (kmax,kfmax),
                (kmax, kf_threshold),
                v3,
              ]

    graph12 = [ Path.MOVETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.LINETO,
                 Path.CLOSEPOLY,
              ]

    path12 = Path(verts12, graph12)

    c12 = centroid.calculate_polygon_centroid(verts12)


# PLOT 'EM ALL
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fig1 = plt.figure(facecolor='white')
    ax1 = Subplot(fig1, 111)
    fig1.add_subplot(ax1)
    ax1.axis["right"].set_visible(False)
    ax1.axis["top"].set_visible(False)
    ax1.axis["left"].set_visible(False)
    ax1.axis["bottom"].set_visible(False)

# Turn off tick labels and tick marks
    ax1.set_xticks([]) 
    ax1.set_yticks([])

# Case 3 region ...
    patch3 = patches.PathPatch(path3, \
                facecolor='0.6', edgecolor='red', lw=0)
    plt.gca().add_patch(patch3)

# Case 4 region ...
    patch4 = patches.PathPatch(path4, \
                facecolor='0.8', edgecolor='red', lw=0)
    plt.gca().add_patch(patch4)

# Case 1 and Case 2 region ...
    patch12 = patches.PathPatch(path12, \
                facecolor='0.9', edgecolor='none', lw=0)

#patch12.set_label('Case 1 - 2')

    plt.gca().add_patch(patch12)

# Intercepts ...
    intercept_y = (0., kf_threshold)

    ax1.annotate(r"$k^{f}(\hat{q},\Pi^f, \beta, \sigma)$", xy= intercept_y,\
                xytext=intercept_y,
                horizontalalignment='right', verticalalignment='top',
                )

    intercept_x = (k_threshold, 0.0)

    ax1.annotate(r"$k(\hat{q},\Pi, \beta, \sigma)$", xy= intercept_x, \
                 xytext=(k_threshold, -0.015),
                horizontalalignment='center', verticalalignment='top',
                )

# Label regions/ patches:
    ax1.annotate(r"Case 1 and Case 2", xy = c12, \
                 xytext = c12,
                horizontalalignment='left', verticalalignment='top',
                )

    ax1.annotate(r"Case 3", xy = c3, \
                 xytext = c3,
                horizontalalignment='left', verticalalignment='top',
                )

    ax1.annotate(r"Case 4", xy = c4*0.25, \
                 xytext = c4*0.25,
                horizontalalignment='left', verticalalignment='bottom',
                )


# CUSTOM: Draw axes with arrow tips:

    ax1.set_xlim(0,kmax*1.1)
    ax1.set_ylim(0,kfmax*1.1)

    plt.arrow(0, 0, kmax*1.02, 0, width=0.0005, color="k", \
                clip_on=False, head_width=0.007, head_length=0.05) # x-axis

    plt.arrow(0, 0, 0, kfmax*1.02, width=0.0005, color="k", \
                clip_on=False, head_width=0.02, head_length=0.015) # y-axis

# CUSTOM: Set axes labels:
#plt.xlabel(r"$\kappa$", fontsize=ftmed)
#plt.ylabel(r"$\kappa^{f}$", fontsize=ftmed)

    xlabel_pos = (kmax*1.05, 0)
    ylabel_pos = (-0.0015,kfmax)

    ax1.annotate(r"$\kappa$", xy = xlabel_pos, \
                 xytext = xlabel_pos,
                horizontalalignment='left', verticalalignment='top',
                )


    ax1.annotate(r"$\kappa^{f}$  ", xy = ylabel_pos, \
                 xytext = ylabel_pos,
                horizontalalignment='right', verticalalignment='bottom',
                )

    ax1.annotate(r"$O$", xy = (0,0), \
                 xytext = (-0.08,-0.03),
                horizontalalignment='left', verticalalignment='bottom',
                )



# Save figures in EPS and PNG formats:
    plt.savefig('regimes1.eps')
    plt.savefig('regimes1.png')

# Display all figures at the end:
    plt.ion() # Tell pylab to run UI in separate thread

    plt.show()

