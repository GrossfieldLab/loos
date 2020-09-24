import loos
# from loos import pyloos

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
# from matploltib.pyplot docs
def heatmap(data, row_labels, col_labels, ax=None, cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right", rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(
    im,
    data=None,
    valfmt="{x:.2f}",
    textcolors=["black", "white"],
    threshold=None,
    **textkw
):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max()) / 2.0

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center", verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


pre = 'clusters/'

roc_exemplar_names = [
    "ROC-kgs-pool-protons-100-0.6626-20.pdb",
    "ROC-kgs-pool-protons-100-0.1851-13.pdb",
    "ROC-kgs-pool-protons-100-0.077-11.pdb",
]

sel = '!hydrogen'
roc_exemplars = [loos.selectAtoms(loos.createSystem(pre+fn), sel)
                 for fn in roc_exemplar_names]
roc_labels = ["ROC 1 (0.6626)", 'ROC 2 (0.1851)', "ROC 3 (0.077)"]
ol3_exemplar_names = [
    "OL3-kgs-pool-protons-100-0.6529-13.pdb",
    "OL3-kgs-pool-protons-100-0.1246-5.pdb",
    "OL3-kgs-pool-protons-100-0.0901-16.pdb"
]
ol3_exemplars = [loos.selectAtoms(loos.createSystem(pre+fn), sel)
                 for fn in ol3_exemplar_names]
ol3_labels = ["OL3 1 (0.6529)", "OL3 2 (0.1246)", "OL3 3 (0.0901)"]

catted_labels = roc_labels + ol3_labels
catted_structures = roc_exemplars + ol3_exemplars
rmsds = np.zeros(36).reshape(6,6)
for i, r in enumerate(catted_structures):
    for j, o in enumerate(catted_structures):
        r.alignOnto(o)
        rmsds[i][j] = r.rmsd(o)

fig, ax = plt.subplots()

im, cbar = heatmap(
    rmsds,
    catted_labels,
    catted_labels,
    ax=ax,
    cmap="Greens",
    cbarlabel="Heavy atom RMSD (Å)"
)

annotations = annotate_heatmap(im, size=MEDIUM_SIZE)
fig.tight_layout()
fig_fn = "slag_plots/clusters-fl.rmsds"
plt.savefig(fig_fn+'.svg', format="svg", transparent=True, bbox_inches="tight")
plt.savefig(fig_fn+'.pdf', format="pdf", transparent=True, bbox_inches="tight")
plt.show()
