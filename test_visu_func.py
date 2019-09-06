import visu_func
from visu_func import Node
from visu_func import Cell
import pandas as pd


def test_cell_variables():
    
    # Testing on mesh 4x4 cells for Noh problem
    names_nodes = ['nx', 'ny']
    data_nodes = pd.read_csv(r"/Users/Ivan/Documents/FJFI/BC_Arbeit/CCLuS/mesh_000.dat",
                       delimiter=' ',
                       skipinitialspace=True,
                       names=names_nodes)    
    nn = data_nodes.shape[0]
    nodes = [Node(i, data_nodes.nx[i-1], data_nodes.ny[i-1]) for i in range(1, nn+1)]
    
    # Cell ID, boundary identificator, no. of vertices, no. of neigh. cells
    names_cell_connectivity = ['cID',
                               'IDbdry',
                               'nvert',
                               'nbdrcells']
    data_cell_connectivity = pd.read_csv(r"/Users/Ivan/Documents/FJFI/BC_Arbeit/CCLuS/input_cellsid.dat",
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
                           'v1', 'v2', 'v3', 'v4', 'v5', 'v6',
                           'ng1', 'ng2', 'ng3', 'ng4', 'ng5', 'ng6']
    
    data_cell_vertices = pd.read_csv(r"/Users/Ivan/Documents/FJFI/BC_Arbeit/CCLuS/input_cnconn.dat",
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
    
    data_variables = pd.read_csv(r"/Users/Ivan/Documents/FJFI/BC_Arbeit/CCLuS/vals_000.dat",
                                 delimiter=' ',
                                 skipinitialspace=True,
                                 names=names_variables)
    
    visu_func.insert_cell_variables(nc,
                                    cells,
                                    data_variables,
                                    data_cell_vertices)
    assert nn == 25
    assert nc == 16
    
    #Counterclockwise node placement test
    assert nodes[0].nID == 1
    assert nodes[0].nx == -1.0
    assert nodes[3].nID == 4
    assert nodes[7].nID == 8
    
    assert cells[0].cID == 1