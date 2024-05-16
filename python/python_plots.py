import pandas as pd
import numpy as np
import os, os.path
import matplotlib
import matplotlib.pyplot as plt
import support_functions as sf

##################################
#    SOME FILE NAME FUNCTIONS    #
##################################

def fnt_img_paired_outcome_histogram_plots(
    grouping: str, 
    outcome_field: str
):
    return f"overlayed_outcome_histogram_{grouping}-{outcome_field}"

def fnt_img_expected_outcome_matrix_plots(
    scenario: int
):
    return f"expected_outcome_matrix_scenario-{scenario}"


############################
#    PLOTTING FUNCTIONS    #
############################

# overlayed histogram plots
def generate_overlayed_outcome_histogram_plots(
        df_outcome_values: pd.DataFrame,
        attr_scenarios: sf.AttributeTable,
        dict_scenario_groups: dict,
        fields_outcome: list,
        dir_export: str,
        n_bins: int = 20,
        line_weight: int = 4,
        filename_template = fnt_img_paired_outcome_histogram_plots,
        extension: str = "jpg"
    ):

    # check output directory
    if not os.path.exists(dir_export):
        os.makedirs(dir_export, exist_ok = True)
    
    
    scenario_groups = list(dict_scenario_groups.keys())
    scenario_groups.sort()
    

    for grouping in scenario_groups:
        
        df_ov_filt = df_outcome_values[df_outcome_values[attr_scenarios.key].isin(dict_scenario_groups[grouping])]
        scens = list(set(df_ov_filt[attr_scenarios.key]))
        scens.sort()
        
        for field_outcome in fields_outcome:
            
            fig, ax = plt.subplots(figsize = (15, 10))
            ax.set_xlabel(field_outcome, size = 15)
            ax.set_ylabel("count", size = 15)
            ax.set_title(f"Histogram Traces: {field_outcome}, {grouping} Scenarios", size = 20)
            for scen in scens:
                
                scen_name = attr_scenarios.field_maps[f"{attr_scenarios.key}_to_scenario_name"].get(scen)
                outcomes = df_ov_filt[df_ov_filt[attr_scenarios.key] == scen][field_outcome]
                
                h = np.histogram(outcomes, bins = n_bins)
                
                vec_bin_edges = h[1]
                x = (vec_bin_edges[0:-1] + vec_bin_edges[1:])/2
                y = h[0]

                ax.plot(x, y, lw = line_weight, label = scen_name)
            
            ax.legend(fontsize = 14)
            #plt.show()
            fn_exp = filename_template(grouping, field_outcome)
            fp_exp = os.path.join(dir_export, f"{fn_exp}.{extension}")
            plt.savefig(fp_exp, dpi = 300, bbox_inches = "tight")

    
    return 0



###   
#  heatmap & annotate_heatmap from: https://matplotlib.org/3.5.0/gallery/images_contours_and_fields/image_annotated_heatmap.html
#  

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
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

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_yticklabels(labels=row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
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
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
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
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
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



# get data for heat maps
def get_transition_outcomes(
    df_in: pd.DataFrame, 
    scen: int, 
    attr_states: sf.AttributeTable, 
    field_scenkey: str = "scenario_id", 
    field_i: str = "i", 
    field_j: str = "j", 
    field_value: str = "expected_outcome_across_edge"
) -> np.ndarray:
    
    sf.check_keys(df_in, [field_i, field_j, field_value])
    m = attr_states.n_key_values
    dict_labels = attr_states.field_maps["state_to_name"]
    
    # initialize output array and get indices for overwrite
    array_out = np.zeros((m , m)).astype(float)
    array_inds = np.array(df_in[(df_in[field_scenkey] == scen) & (~df_in[field_value].isna())][[field_i, field_j, field_value]])
    inds = np.array(array_inds[:, 0] - 1)*m + (array_inds[:, 1]) - 1
    inds = inds.astype(int)
    
    # get values and overwrite in array
    vals = np.nan_to_num(array_inds[:, 2])
    np.put(array_out, inds, vals)
    
    # find states to keep
    states_keep = list(set(array_inds[:, 0].astype(int) - 1) | set(array_inds[:, 1].astype(int) - 1))
    states_keep.sort()
    
    labels_keep = [dict_labels.get(x + 1) for x in states_keep]
    
    #return array_out
    return array_out[states_keep, :][:, states_keep], labels_keep



# heat maps
def generate_expected_value_matrix_plots(
        df_expected_outcomes: pd.DataFrame,
        attr_scenarios: sf.AttributeTable,
        attr_states: sf.AttributeTable,        
        dir_export: str,
        filename_template = fnt_img_expected_outcome_matrix_plots,
        extension: str = "jpg",
        c_map: str = "RdYlGn",
        **kwargs
    ):

    # check output directory
    if not os.path.exists(dir_export):
        os.makedirs(dir_export, exist_ok = True)
    
    # ranges for color bars (fixed)
    rng_val = np.max(np.abs(np.nan_to_num(np.array(df_expected_outcomes["expected_outcome_across_edge"]))))

    # scenarios to include
    scens = list(set(df_expected_outcomes[attr_scenarios.key]))
    scens.sort()
        
    for scen in scens:
        # get relevant information for this scenario
        gto = get_transition_outcomes(
            df_expected_outcomes, 
            scen,
            attr_states
        )
        
        # scenario name
        scen_name = attr_scenarios.field_maps[f"{attr_scenarios.key}_to_scenario_name"].get(scen)
        
        # set up plot 
        fig, ax = plt.subplots(figsize = (18, 18))
        ax.set_title(f"Expected Outcomes by Edge Traversal - {scen_name}\n", size = 18)
        
        im, cbar = heatmap(
            gto[0], gto[1], gto[1], 
            ax=ax,
            cmap=c_map, cbarlabel="Expected Outcomes for Pathways that Traverse Edge i<->j ", 
            vmin = -rng_val, vmax = rng_val
        )
        texts = annotate_heatmap(im, valfmt="{x:.1f}")
        fig.tight_layout()
        
        # set up and export
        fn_exp = filename_template(scen)
        fp_exp = os.path.join(dir_export, f"{fn_exp}.{extension}")
        plt.savefig(fp_exp, dpi = 300, bbox_inches = "tight")

    
    return 0
