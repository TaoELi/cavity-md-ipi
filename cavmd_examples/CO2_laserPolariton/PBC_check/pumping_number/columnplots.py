'''
One part of sciplot package, plot everthing in one figure.
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
#mpl.use('pdf')
def initialize(col=1, row=1, width=4,
               height=None,
               sharex=False,
               sharey=False,
               commonX=None,
               commonY=None,
               commonYs=None,
               labelthem=None,
               labelthemPosition=None,
               fontsize=8,
               return_fig_args=False,
               LaTeX=False
               ):
    plt.rc('font', family='sans-serif')#, serif='Times')
    plt.rc('text', usetex=LaTeX)
    plt.rc('xtick', labelsize=fontsize)
    plt.rc('ytick', labelsize=fontsize)
    plt.rc('axes', labelsize=fontsize)
    if height is None:
        height = width / 1.618
    fig, axes = plt.subplots(col, row, sharex=sharex, sharey=sharey)
    fig.set_size_inches(width, height)
    if commonX is not None:
        fig.text(commonX[0], commonX[1], commonX[2], ha='center', fontsize=fontsize)
    if commonY is not None:
        fig.text(commonY[0], commonY[1], commonY[2],  va='center', rotation='vertical', fontsize=fontsize)
    if commonYs is not None:
        for i in range(row):
            fig.text(commonYs[i][0], commonYs[i][1], commonYs[i][2],  va='center', rotation='vertical')
    x0, y0 = -0.10, 1.15
    if labelthemPosition is not None:
        x0, y0 = labelthemPosition
    if labelthem is True:
        if row is 1:
            label1List = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)"]
            for i in range(col):
                axes[i].text(x0, y0, label1List[i], transform=axes[i].transAxes, fontsize=fontsize+5, fontweight='bold', va='top', ha='right')
        elif row is 2:
            label2List = [ ["(a)", "(b)"], ["(c)", "(d)"], ["(e)", "(f)"], ["(g)", "(h)"], ["(i)", "(j)"], ["(k)", "(l)"]]
            #if LaTeX:
            #    label2List = [ [r"\textbf{a}", r"\textbf{b}"], [r"\textbf{c}", r"\textbf{d}"], [r"\textbf{e}", r"\textbf{f}"], [r"\textbf{g}", r"\textbf{h}"] ]
            if col is 1:
                i = 0
                for j in range(row):
                    axes[j].text(x0, y0, label2List[i][j], transform=axes[j].transAxes, fontsize=fontsize+5, fontweight='bold', va='top', ha='right')
            else:
                for i in range(col):
                    for j in range(row):
                        axes[i,j].text(x0, y0, label2List[i][j], transform=axes[i, j].transAxes, fontsize=fontsize+5, fontweight='bold', va='top', ha='right')
        elif row is 3:
            label3List = [ ["(a)", "(b)", "(c)"], ["(d)", "(e)", "(f)"], ["(g)", "(h)", "(i)"]]
            if col is 1:
                for i in range(row):
                    axes[i].text(x0, y0, label3List[0][i], transform=axes[i].transAxes, fontsize=fontsize+5, fontweight='bold', va='top', ha='right')
            else:
                for i in range(col):
                    for j in range(row):
                        axes[i,j].text(x0, y0, label3List[i][j], transform=axes[i, j].transAxes, fontsize=fontsize+5, fontweight='bold', va='top', ha='right')
        elif row is 4:
            label2List = [ ["(a)", "(b)", "(c)", "(d)"], ["(e)", "(f)", "(g)", "(h)"], ["(i)", "(j)", "(k)", "(l)"] ]
            if col == 1:
                i = 0
                for j in range(row):
                    axes[j].text(x0, y0, label2List[i][j], transform=axes[j].transAxes, fontsize=fontsize+5, fontweight='bold', va='top', ha='right')
            else:
                for i in range(col):
                    for j in range(row):
                        axes[i,j].text(x0, y0, label2List[i][j], transform=axes[i, j].transAxes, fontsize=fontsize+5, fontweight='bold', va='top', ha='right')
    if return_fig_args:
        return fig, axes
    else:
        return axes

def plotone(
    xs,
    ys,
    ax,
    colors=None,
    labels=None,
    lw=2,
    markersize=None,
    xlabel=None,
    ylabel=None,
    xlog=False,
    ylog=False,
    xlim=None,
    ylim=None,
    showlegend=True,
    alphaspacing=0.0,
    alpha=1.0,
    mfc="none",
    bothyticks=True,
    yscientific=False,
    yscientificAtLabel=False,
    sharewhichx=None,
    sharewhichy=None,
    legendFontSize=None,
    rainbowColor=False,
        ):
    # Find the scienfic order for y axis
    if yscientificAtLabel is True:
        if ylim is not None:
            yorder = np.floor(np.log10(np.max(np.abs(ylim[1]))))
        else:
            yorder = np.floor(np.log10(np.max(np.abs(ys[0]))))
    else:
        yorder = 0
    lines = []
    if colors is None:
        for i in range(len(xs)):
            if labels is None:
                line, = ax.plot(xs[i], ys[i]/10**yorder, lw=lw, markersize=markersize, mfc=mfc, alpha=1.0 - i*alphaspacing if alphaspacing else alpha)
            else:
                line, = ax.plot(xs[i], ys[i]/10**yorder, lw=lw, label=labels[i], markersize=markersize, mfc=mfc, alpha=1.0 - i*alphaspacing if alphaspacing else alpha)
            lines.append(line)
    else:
        for i in range(len(xs)):
            if labels is None:
                line, = ax.plot(xs[i], ys[i]/10**yorder, colors[i], lw=lw, markersize=markersize, mfc=mfc, alpha=1.0 - i*alphaspacing if alphaspacing else alpha)
            else:
                line, = ax.plot(xs[i], ys[i]/10**yorder, colors[i], lw=lw, label=labels[i], markersize=markersize, mfc=mfc, alpha=1.0 - i*alphaspacing if alphaspacing else alpha)
            lines.append(line)
    if yscientific is True:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    # Add xlabel and ylabel
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        if yorder is not 0:
            ax.set_ylabel(ylabel + " [$\\times$ 10$^{%d}$]" %yorder)
        else:
            ax.set_ylabel(ylabel)
    # Set xlim or ylim
    if xlim is not None:
        ax.set_xlim(xlim[0], xlim[1])
    if ylim is not None:
        ax.set_ylim(ylim[0]/10**yorder, ylim[1]/10**yorder)
    # Choose the log or normal sacling for axes
    if xlog:
        ax.set_xscale('log')
    if ylog:
        ax.set_yscale('log')
    # Add twin axes
    if bothyticks:
        ax.tick_params(direction='in')
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
    # Save file
    if showlegend:
        ax.legend(fontsize=legendFontSize, markerscale=legendFontSize)
    if sharewhichx is not None:
        ax.get_shared_x_axes().join(sharewhichx, ax)
        sharewhichx.set_xticklabels([])
    if sharewhichy is not None:
        ax.get_shared_y_axes().join(sharewhichy, ax)
        ax.set_yticklabels([])
    if rainbowColor:
        colormap = plt.cm.gist_rainbow #nipy_spectral, Set1,Paired
        colors = [colormap(i) for i in np.linspace(0, 1,len(ax.lines))]
        for i,j in enumerate(ax.lines):
            j.set_color(colors[i])
    return lines

def broken_y(axes, d=0.015):
    axes[0].spines['right'].set_visible(False)
    axes[1].spines['left'].set_visible(False)
    axes[0].yaxis.tick_left()
    axes[0].tick_params(labelright='off')
    axes[1].yaxis.tick_right()

    d = d # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=axes[0].transAxes, color='k', clip_on=False)
    axes[0].plot((1-d,1+d), (-d,+d), **kwargs)
    axes[0].plot((1-d,1+d),(1-d,1+d), **kwargs)

    kwargs.update(transform=axes[1].transAxes)  # switch to the bottom axes
    axes[1].plot((-d,+d), (1-d,1+d), **kwargs)
    axes[1].plot((-d,+d), (-d,+d), **kwargs)

def adjust(
    wspace=0,
    hspace=0,
    tight_layout=False,
    savefile=None,
    includelegend=None
    ):
    #plt.subplots_adjust(wspace=wspace, hspace=hspace)
    if tight_layout is True:
        plt.tight_layout()
    if savefile is not None and includelegend is None:
        plt.savefig(savefile, dpi=300, bbox_inches="tight")
    if savefile is not None and includelegend is not None:
        plt.savefig(savefile, dpi=300, bbox_extra_artists=(includelegend,), bbox_inches='tight')
    # show figure
    plt.show()

if __name__ == '__main__':
    x = np.linspace(0, 1e4, 100)
    y = np.exp(-1e-4*x)
    y2 = np.cos(x)
    axes = initialize(2,1)
    plotone([x, x], [y, y2], axes[0], labels=['data 1', 'data 2'],  ylabel='$P_2$')
    plotone([x, x], [y*0.5, y2*2], axes[1], labels=['data 1', 'data 2'], xlabel='time', ylabel='$P_2$')
    adjust(tight_layout=True)
