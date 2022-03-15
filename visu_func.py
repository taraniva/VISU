import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib as mpl
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Classes


class Node:
    def __init__(self, nID, nx, ny):
        self.nID = nID
        self.nx = nx
        self.ny = ny


class Cell:
    def __init__(self, cID, nvert, itr_nodes, ctr_vals):
        self.cID = cID
        self.nvert = nvert
        self.itr_nodes = itr_nodes
        self.ctr_vals = ctr_vals

# Functions


def truncate_colourmap(cmap, minval, maxval, n=100):
    """Modifies the span of inserted colourmap.
    
    :param cmap: Matplotlib `colormap` instance. 
	Example colourmaps eg. 'jet', 'viridis' etc.
    :param minval: `float` in interval [0,1] representing lower boundary colour
    :param maxval: `float` in interval [0,1] representing upper boundary colour
    
    :return: Colourmap with updated colour span.
    :rtype: Matplotlib `colourmap` instance.
    """
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def get_cvert_coords_x(cell_ID, nodes, cells):
    """Gets x-coordinates of all nodes in a given cell.
    
    :param cell_ID: `integer` representing cell ID
    :param nodes: List containing objects of `node` class
    :param cells: List containing objects of `cell` class
    
    :return: `list` containing x-coordinates of nodes in cell nr. `cell_ID`
    """
    coord_x = [nodes[i-1].nx for i in cells[cell_ID].itr_nodes]

    return coord_x


def get_cvert_coords_y(cell_ID, nodes, cells):
    """Gets y-coordinates of all nodes in a given cell.
    
    :param cell_ID: `integer` representing cell ID
    :param nodes: List containing objects of `node` class
    :param cells: List containing objects of `cell` class
    
    :return: `list` containing y-coordinates of nodes in cell nr. `cell_ID`
    """
    coord_y = [nodes[j-1].ny for j in cells[cell_ID].itr_nodes]

    return coord_y


# Unused at the moment within the draw routine
def connect_nodes(figure_nodes, first_node_ID, second_node_ID, nodes):
    """Plots a line between two selected nodes.

    :param figure_nodes: Matplotlib `figure` instance
    :param first_node_ID: `integer` ID of first node
    :param first_node_ID: `integer` ID of second node
    :param nodes: `list` of nodes
    """
    figure_nodes.plot([nodes[first_node_ID].nx, nodes[second_node_ID].nx],
                      [nodes[first_node_ID].ny, nodes[second_node_ID].ny],
                      'ro-')


def insert_cell_variables(nc,
                          cells,
                          data_variables,
                          data_cell_vertices):
    """Transfers cell variables from imported data tables into the cells.
    
    :param nc: `integer` representing the total number of cells
    :param cells: `list` containing all `cell` objects
    :param data_variables: `pandas` table containing all cell variables ordered by cell ID
    :param data_cell_vertices: `pandas` table containing list of vertices in all cells
    """

    for k in range(nc):

        if cells[k].nvert == 3 or 4 or 5 or 6 or 7 or 8 or 9:
            cells[k].itr_nodes = [data_cell_vertices[f"v{i}"][k]
                                  for i in range(cells[k].nvert, 0, -1)]
            cells[k].itr_nodes = [int(l) for l in cells[k].itr_nodes]

        else:
            raise ValueError('Not implemented for ',
                             cells[k].nvert, ' vertices of cell nr.',
                             cells[k].cID, '.')

    # Insert values into cell-centres
    for k in range(nc):
        cells[k].ctr_vals = [data_variables.cx[k],
                             data_variables.cy[k],
                             data_variables.cu[k],
                             data_variables.cv[k],
                             data_variables.cell_density[k],
                             data_variables.cell_pressure[k],
                             data_variables.cell_internal_energy[k]]


##############################################################################
#                  M A I N   P L O T T I N G   R O U T I N E
##############################################################################

def draw(filename_mesh,
         filename_cellsid,
         filename_connectivity,
         filename_vals,
         val_name,
         x_min,
         x_max,
         y_min,
         y_max,
         c_min,
         c_max):
    """Main plotting routine. Generates nodal plot, coloured plot and mesh plot
    and saves them in `.eps` format.

    :param filename_mesh: Path to file containing nodal coordinates
    :param filename_cellsid: Path to file containing cell characteristics
    :param filename_connectivity: Path to file containing cell connectivity
    :param filename_vals: Path to file containing cell-centred variables
    :param val_name: Varible for coloured plot, either `density`, `pressure` or `internal_energy`
    :param x_min: "left" x-axis boundary for plotting
    :param x_max: "right" x-axis boundary for plotting
    :param y_min: "lower" y-axis boundary for plotting
    :param y_max: "upper" y-axis boundary for plotting
    """

    ##########################################################################
    #                   D  A  T  A      I  M  P  O  R  T
    ##########################################################################
    names_nodes = ['nx', 'ny']
    data_nodes = pd.read_csv(filename_mesh,
                             delimiter=' ',
                             skipinitialspace=True,
                             names=names_nodes)
    nn = data_nodes.shape[0]
    nodes = [Node(i, data_nodes.nx[i-1], data_nodes.ny[i-1])
             for i in range(1, nn+1)]

    # Cell ID, boundary identificator, no. of vertices, no. of neigh. cells
    names_cell_connectivity = ['cID',
                               'IDbdry',
                               'nvert',
                               'nbdrcells']
    data_cell_connectivity = pd.read_csv(filename_cellsid,
                                         delimiter=' ',
                                         skipinitialspace=True,
                                         names=names_cell_connectivity)
    nc = data_cell_connectivity.shape[0]
    cells = [Cell(data_cell_connectivity.cID[k-1],
                  data_cell_connectivity.nvert[k-1],
                  [],
                  []) for k in range(1, nc+1)]

    # Internal nodes and neigh. cells by ID for each cell
    names_cell_vertices = ['cID',
                           'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7', 'v8', 'v9',
                           'ng1', 'ng2', 'ng3', 'ng4', 'ng5', 'ng6', 'ng7', 'ng8', 'ng9']

    data_cell_vertices = pd.read_csv(filename_connectivity,
                                     delimiter=' ',
                                     skipinitialspace=True,
                                     names=names_cell_vertices)

    # Cell-centred values
    names_variables = ['cx',
                       'cy',
                       'cu',
                       'cv',
                       'cell_density',
                       'cell_pressure',
                       'cell_internal_energy']

    data_variables = pd.read_csv(filename_vals,
                                 delimiter=' ',
                                 skipinitialspace=True,
                                 names=names_variables)
    ##########################################################################
    #                          U  T  I  L  I  T  Y
    ##########################################################################
    # Fills cells with lists of their vertices by ID
    # !!!!!! Vertices MUST BE INSERTED CLOCKWISE !!!!!!

    insert_cell_variables(nc,
                          cells,
                          data_variables,
                          data_cell_vertices)
    #########################################################################
    # Plot nodes to a separate figure
    x = [nodes[i].nx for i in range(nn)]
    y = [nodes[j].ny for j in range(nn)]
    node_plot_onoff = 1
    if node_plot_onoff == 1:
        figure_nodes = plt.figure()
        node_plot = figure_nodes.add_subplot(111)
        node_plot.scatter(x, y)
        node_plot.set_title('VISU Nodes')
        node_plot.set_xlabel('x')
        node_plot.set_ylabel('y')
        node_plot.set_xlim([x_min, x_max])
        node_plot.set_ylim([y_min, y_max])
        for k in range(nn):
            print_string = str(nodes[k].nID)
            node_plot.annotate(print_string, (x[k], y[k]))
        figure_nodes.savefig('nodes_plot.pdf', format='pdf')
        plt.show()

    #######################################################################
    #                             P  L  O  T
    #######################################################################
    figure_colour_plot, colour_plot = plt.subplots()
    divider = make_axes_locatable(colour_plot)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    figure_mesh_plot, mesh_plot = plt.subplots()
    figure_contour_plot, contour_plot = plt.subplots()

    colour_plot.set_xlim([x_min, x_max])
    colour_plot.set_ylim([y_min, y_max])

    mesh_plot.set_xlim([x_min, x_max])
    mesh_plot.set_ylim([y_min, y_max])

    contour_plot.set_xlim([x_min, x_max])
    contour_plot.set_ylim([y_min, y_max])

    patches = []

    # Determine which variable to plot
    if val_name == 'density':
        selector = 4
    elif val_name == 'pressure':
        selector = 5
    elif val_name == 'internal_energy':
        selector = 6
    else:
        raise ValueError('Not a variable: ', val_name,)

    for k in range(nc):
        cx = get_cvert_coords_x(k, nodes, cells)
        cy = get_cvert_coords_y(k, nodes, cells)
        points = np.column_stack((cx, cy))

        polygon = mpl.patches.Polygon(points)
        patches.append(polygon)

        colors = [cells[k].ctr_vals[selector] for k in range(nc)]
        colors = np.array(colors)

    cell_center_x = [cells[k].ctr_vals[1] for k in range(nc)]
    cell_center_y = [cells[k].ctr_vals[2] for k in range(nc)]

    # For debug
    # print(colors)

    # General plotting settings
    plt.rcParams["figure.figsize"] = [12*2, (12-(0.05*12))*2]
    # plt.rcParams["figure.figsize"] = [7, 3-(0.05*3)]
    font = {'family': 'Sans',
            'weight': 'normal',
            'size': 14}

    mpl.rc('font', **font)

    """
    #T I C K S
    mpl.rcParams['xtick.major.pad'] = 60
    mpl.rcParams['ytick.major.pad'] = 60
    """

    ######################### Coloured plot #########################
    custom_cmap = truncate_colourmap(plt.get_cmap('jet'), 0.20, 0.95)
    p = PatchCollection(patches, cmap=custom_cmap, alpha=0.8)
    p.set_array(colors)
    p.set_edgecolor("k")
    colour_plot.add_collection(p)
    figure_colour_plot.colorbar(p, cax=cax)
    colour_plot.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))
    colour_plot.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))
    colour_plot.xaxis.set_ticks_position('bottom')
    colour_plot.yaxis.set_ticks_position('left')

    # Span of colourbar values
    p.set_clim([c_min, c_max])

    #plt.rcParams.update({'font.size': 14})
    #colour_plot.autoscale_view()
    #colour_plot.set_title('', **font)
    #colour_plot.set_xlabel('x [cm]', **font)
    #colour_plot.set_ylabel('y [cm]', **font)
    figure_colour_plot.savefig(
		'colour_plot.pdf', dpi=800, bbox_inches='tight', pad_inches=0)
	# figure_colour_plot.savefig(
	# 'colour_plot.eps',dpi=300,bbox_inches='tight', pad_inches=0)

    ######################### Mesh plot #########################
    q = PatchCollection(patches, alpha=0.8)
    q.set_edgecolor("k")
    q.set_facecolor("cyan")
    mesh_plot.add_collection(q)

    mesh_plot.autoscale_view()
    mesh_plot.set_xlabel('x')
    mesh_plot.set_ylabel('y')

    figure_mesh_plot.savefig(
		'grid_plot.pdf', dpi=800, bbox_inches='tight', pad_inches=0)

    ######################### Contour plot #########################
    contour_onoff = 0
    if contour_onoff == 1:
        contour_plot.tricontour(
            cell_center_x, cell_center_y, colors, levels=14, linewidths=0.5, colors='k')
        cntr2 = contour_plot.tricontourf(
            cell_center_y, cell_center_x, colors, levels=16, cmap="RdBu_r")

        figure_contour_plot.colorbar(cntr2, ax=contour_plot)
        #contour_plot.plot(x, y, 'ko', ms=3)
        contour_plot.set(xlim=(x_min, x_max), ylim=(y_min, y_max))
        figure_contour_plot.savefig('contour_plot.png')

    plt.show()

##########################################################################
#                      E N D   O F   M O D U L E
##########################################################################


"""
W I P
"""

#def import_mesh_data(filename_mesh,
#                     filename_cellsid,
#                     filename_connectivity,
#                     filename_vals,
#                     nn,
#                     nodes,
#                     nc,
#                     cells,
#                     data_cell_vertices,
#                     data_variables):
#
#    names_nodes = ['nx', 'ny']
#    data_nodes = pd.read_csv(filename_mesh,
#                       delimiter=' ',
#                       skipinitialspace=True,
#                       names=names_nodes)
#    nn = data_nodes.shape[0]
#    print(nn)
#    nodes = [Node(i, data_nodes.nx[i-1], data_nodes.ny[i-1]) for i in range(1, nn+1)]
#
#    # Cell ID, boundary identificator, no. of vertices, no. of neigh. cells
#    names_cell_connectivity = ['cID',
#                               'IDbdry',
#                               'nvert',
#                               'nbdrcells']
#    data_cell_connectivity = pd.read_csv(filename_cellsid,
#                                         delimiter=' ',
#                                         skipinitialspace=True,
#                                         names=names_cell_connectivity)
#    nc = data_cell_connectivity.shape[0]
#    cells = [Cell(data_cell_connectivity.cID[k-1],
#                  data_cell_connectivity.nvert[k-1],
#                  [],
#                  []) for k in range(1, nc+1)]
#
#
#    # Internal nodes and neigh. cells by ID for each cell
#    names_cell_vertices = ['cID',
#                           'v1', 'v2', 'v3', 'v4', 'v5', 'v6',
#                           'ng1', 'ng2', 'ng3', 'ng4', 'ng5', 'ng6']
#
#    data_cell_vertices = pd.read_csv(filename_connectivity,
#                                     delimiter=' ',
#                                     skipinitialspace=True,
#                                     names=names_cell_vertices)
#
#    # Cell-centred values
#    names_variables = ['cx',
#                       'cy',
#                       'cu',
#                       'cv',
#                       'cell_density',
#                       'cell_pressure',
#                       'cell_internal_energy']
#
#    data_variables = pd.read_csv(filename_vals,
#                                 delimiter=' ',
#                                 skipinitialspace=True,
#                                 names=names_variables)
