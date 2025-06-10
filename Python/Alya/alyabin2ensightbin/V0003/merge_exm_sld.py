"""
pass all exmedi arrays to solidz mesh
!!!ONLY TRIANGLE MESHES!!!!
!!!Expects that in the file of every variable in ensi, the comment lime start with t=time, time - time that corresponds to this variable
"""

import vtk
import sys
import os

mpi_on = True
try:
    mpi_on = int(os.environ["MPI_SINGLETON"]) != 1
except:
    mpi_on = True

if mpi_on:
    from mpi4py import MPI  # needs install


import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from progressbar import progressbar, ProgressBar
import shutil
import glob
import re
import pandas as pd
import itertools

if mpi_on:
    comm = MPI.COMM_WORLD
    my_rank = comm.Get_rank()
    my_name = MPI.Get_processor_name()
    num_procs = comm.Get_size()
else:
    comm = None
    my_rank = 0
    my_name = "nompi"
    num_procs = 1

variables_to_ignore = ["ISOCH"]

casenames = None
coarse_mesh = None
out_path = None
out_casename = None


if my_rank == 0:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("coarse_mesh", help="Coarse mesh in VTP")
    parser.add_argument("out_path", help="Ouptut path")
    parser.add_argument("out_casename", help="Ouptut case name (prefix)")
    parser.add_argument(
        "--case", action="append", help="Path and case file", required=True
    )
    parser.add_argument(
        "-i",
        "--interval",
        help="Postprocess only timesteps in the specified time interval",
        type=float,
        nargs=2,
        metavar=("time_start", "time_end"),
    )
    args = parser.parse_args()

    casenames = args.case
    coarse_mesh = args.coarse_mesh
    out_path = args.out_path
    out_casename = args.out_casename
    time_interval = args.interval


if num_procs > 1:
    casenames = comm.bcast(casenames, root=0)
    coarse_mesh = comm.bcast(coarse_mesh, root=0)
    out_path = comm.bcast(out_path, root=0)
    out_casename = comm.bcast(out_casename, root=0)
    time_interval = comm.bcast(time_interval, root=0)


displ_varname = "DISPL"  # variable name for displacement
time_check_tolerance = 1e-7  # tolerance when time values are compared
time_ndigits = 6  # number of digits to keep after comma for time
iterationid_number_of_digits = (
    6  # numnber of digits in ensi variable names used for number of iteration
)


# ==============================================================
#
#   FUNCTIONS
#


def calculate_interpolation_weights(mesh, mesh_coarse):
    """Calculate intepolation matrices, !!!!!!ONLY TRIANGLES!!!!!! AND ONLY POINT DATA !!!!

    mesh, mesh_coarse -- vtkPolyData

    For each point of mesh_coarse use vtkCellLocator to find closest point on mesh, interpolation weights and point ids for interpolation

    Return interp_weights - npts x 3 matrix of weights, interp_ptids - npts x 3 matrix of point ids
    so that interpolated value of array X is:
    X[i] = interp_weights[i,0] * a[interp_ptids[i,0]] + interp_weights[i,1] * a[interp_ptids[i,1]] + interp_weights[i,2] * a[interp_ptids[i,2]]
    or
    X = np.sum(interp_weights * a[interp_ptids], axis=1)

    """

    npts = mesh_coarse.GetNumberOfPoints()
    interp_weights = np.zeros((npts, 3), dtype=np.float32)
    interp_ptids = np.zeros((npts, 3), dtype=np.int64)

    loc = vtk.vtkCellLocator()
    loc.SetDataSet(mesh)
    loc.BuildLocator()

    for i in progressbar(range(npts)):
        # for i in range(npts):
        closest_point = [0, 0, 0]
        cellId = vtk.mutable(0)
        subId = vtk.mutable(0)
        dist2 = vtk.mutable(0)
        loc.FindClosestPoint(
            mesh_coarse.GetPoint(i), closest_point, cellId, subId, dist2
        )

        # get interpolation weights
        cell = mesh.GetCell(cellId)
        closestPoint1 = [0, 0, 0]
        pcoords = [0, 0, 0]
        weights = [0] * cell.GetNumberOfPoints()
        reval = cell.EvaluatePosition(
            closest_point, closestPoint1, subId, pcoords, dist2, weights
        )
        assert reval >= 0, "mesh.GetCell() failed for unknown reason"

        interp_weights[i, :] = weights

        # get point ids
        ptids = cell.GetPointIds()
        for j in range(ptids.GetNumberOfIds()):
            interp_ptids[i, j] = ptids.GetId(j)

    return interp_weights, interp_ptids


def interpolate_array(array, interp_weights, interp_ptids):
    """Interpolate array "array"

    array          : np.ndarray(K) or np.ndarray((K,N))  vector or columns of values to interpolate from. Must contain at least max(interp_ptids.ravel()) values
    interp_weights : np.ndarray((npts,3)) matrix of weights
    interp_ptids   : np.ndarray((npts,3)) matrix of point ids

    returns np.ndarray( npts ) of interpolate values
    """

    result = None

    if array.ndim == 1:
        result = np.sum(interp_weights * array[interp_ptids], axis=1)
    elif array.ndim == 2:
        result = np.zeros((interp_weights.shape[0], array.shape[1]))
        for i in range(array.shape[1]):
            result[:, i] = np.sum(interp_weights * array[interp_ptids, i], axis=1)

    assert (
        not result is None
    ), f"interpolate_array: array dimension has to be 1 or 2, passed {array.ndim}, shape: {array.shape}"

    return result


#
#  END OF FUNTIONS
#
# ===============================================================

# ==============================================================
#
#   ENSI Writer. !!!ONLY TRIANGLES!!!!
#


def write_geometry(surf_mesh: vtk.vtkPolyData, path: str, casename: str) -> None:
    point_coordinates = vtk_to_numpy(surf_mesh.GetPoints().GetData())
    connectivity = np.reshape(
        vtk_to_numpy(surf_mesh.GetPolys().GetData()), (-1, 4)
    )  # N x 4, columnn 0 == 3, i.e. 3 points, columns 1:4 are point ids
    assert all(
        connectivity[:, 0] == 3
    ), "The coarse mesh has elements that are not triangles"
    connectivity = connectivity[:, 1:]  # discard the 1st column

    # Ensight groups elements by type. Create grouping
    element_alya2ensi = {
        37: {"Name": b"hexa8", "Vertices": 8},
        30: {"Name": b"tetra4", "Vertices": 4},
        32: {"Name": b"pyramid5", "Vertices": 5},
        34: {"Name": b"penta6", "Vertices": 6},
        10: {"Name": b"tria3", "Vertices": 3},
        12: {"Name": b"quad4", "Vertices": 4},
    }

    # geometry ensight
    with open(os.path.join(path, f"{casename}.ensi.geo"), "wb") as f:
        f.write(b"C Binary".ljust(80))
        f.write(b"description line 1".ljust(80))
        f.write(b"description line 2".ljust(80))
        f.write(b"node id given".ljust(80))
        f.write(b"element id given".ljust(80))
        f.write(b"part".ljust(80))
        f.write(np.array([1], dtype=np.int32))  # int
        f.write(b"description line 1".ljust(80))
        f.write(b"coordinates".ljust(80))

        number_of_points = point_coordinates.shape[0]
        f.write(np.array([number_of_points], dtype=np.int32))  # int
        f.write(np.arange(1, number_of_points + 1, dtype=np.int32))

        ndim = point_coordinates.shape[1]

        for iii in range(3):
            f.write(point_coordinates[:, iii].ravel().astype(np.float32))  # x coord

        number_of_elements = connectivity.shape[0]
        # print("Number of elements = ",number_of_elements)

        f.write(element_alya2ensi[10]["Name"].ljust(80))  # name of the elements
        f.write(np.array([number_of_elements], dtype=np.int32))  # int
        f.write(np.arange(1, number_of_elements + 1, dtype=np.int32))
        f.write(connectivity.ravel().astype(np.int32) + 1)


def write_ensi_case(df, timesets):
    """
    variables : dict of {varname : timeset}

    df : dataframe with ['Variable','FileNo','Time','ndim']
    """

    case_file = os.path.join(out_path, f"{out_casename}.ensi.case")
    print(f"Writing {case_file}")

    with open(case_file, "w") as f:
        f.write("# Converted from Alya\n")
        f.write("# Ensight Gold Format\n")
        f.write("#\n")
        f.write(f"# Problem name: {out_casename}\n")
        f.write("FORMAT\n")
        f.write("type: ensight gold\n")
        f.write("\n")
        f.write("GEOMETRY\n")
        f.write(f"model: 1 {out_casename}.ensi.geo\n")
        f.write("\n")
        f.write("VARIABLE\n")

        variable_tuples = [[varname, timeset] for varname, timeset in timesets.items()]
        timeset_df = pd.DataFrame(
            variable_tuples, columns=["Variable", "Timeset"]
        ).sort_values(by="Timeset", ascending=True, axis=0)

        for varname, timeset in timesets.items():
            ndim = df[df["Variable"] == varname]["ndim"].to_numpy()[0]
            vartype = "scalar" if ndim == 1 else "vector"
            f.write(
                f"{vartype} per node: {timeset} {varname} {out_casename}.ensi.{varname}-"
                + "*" * iterationid_number_of_digits
                + "\n"
            )

        # to make sure the whole array gets printed
        # otherwise numpy converts to string a summary, e.g. (1,2,3,...,5,6,7)
        np.set_printoptions(threshold=np.inf)

        f.write("\n")
        f.write("TIME\n")

        for timeset in timeset_df["Timeset"].unique():
            varname = timeset_df[timeset_df["Timeset"] == timeset][
                "Variable"
            ].to_numpy()[
                0
            ]  # get one variable with the timeset

            var_df = df[df["Variable"] == varname].sort_values(
                by="Time", ascending=True, axis=0
            )

            f.write(f"time set: {timeset}\n")
            f.write(f"number of steps: {var_df.shape[0]}\n")
            f.write(f"filename numbers: \n")
            f.write(str(var_df["FileNo"].to_numpy())[1:-1] + "\n")
            f.write("time values:\n")
            f.write(str(var_df["Time"].to_numpy())[1:-1] + "\n")


def write_variable_pernode(data, outputfolder, project_name, varname, iteration):
    # variable ensight
    iterationid_number_of_digits = 6
    fmt = "%s.ensi.%s-" + f"%0{iterationid_number_of_digits}d"
    filename = os.path.join(outputfolder, fmt % (project_name, varname, iteration))

    write_variable_pernode(data, filename)


def write_variable_pernode(data, filename):
    # variable ensight

    with open(filename, "wb") as f:
        f.write(b"description line 1".ljust(80))
        f.write(b"part".ljust(80))
        f.write(np.array([1], dtype=np.int32))  # int
        f.write(b"coordinates".ljust(80))

        if data.ndim == 1:
            f.write(data)
        elif data.ndim == 2:
            # data has coordinates in the order [[x,y,z],[x,y,z],...]
            # expected order of coordinates
            # vx_n1 vx_n2 ... vx_nn nn floats
            # vy_n1 vy_n2 ... vy_nn nn floats
            # vz_n1 vz_n2 ... vz_nn nn floats
            # Rearrange the  matrix
            f.write(data.ravel(order="F").astype(np.float32))
        else:
            assert (
                False
            ), f"write_variable_pernode: array must have dimension 1 or 2, got {data.ndim}, shape {data.shape}"


def get_case_name(filename):
    assert filename.endswith(
        ".ensi.case"
    ), f"get_case_name was passed not ensi.case, passed: {filename}"
    return filename.replace(".ensi.case", "")


def read_variable_pernode(filename) -> (float, np.ndarray):
    """Returns time and array values stored in this file
    Variable file has to have t=time and type=scalar|vector in the header
    """

    # variable ensight

    with open(filename, "rb") as f:
        description = f.read(80).strip()
        part_word = f.read(80).strip()
        part_number = np.fromfile(f, count=1, dtype=np.int32)[0]
        coord_keyword = f.read(80).strip()

        time_part, type_part = description.split(b",")

        assert (
            time_part[:2] == b"t="
        ), f"File {filename} has no time value in the description line"
        assert (
            type_part[:5] == b"type="
        ), f"File {filename} has no type in the description line"
        assert part_word == b"part", f"File {filename} is not a valid ensi variable"
        assert (
            coord_keyword == b"coordinates"
        ), f"File {filename} is not a valid ensi variable"
        assert part_number == 1, f"File {filename} is not a valid ensi variable"

        time = float(time_part.split(b"=")[1])
        vartype = type_part.split(b"=")[1]
        array = None
        if vartype == b"scalar":
            array = np.fromfile(f, dtype=np.float32)
        elif vartype == b"vector":
            array = np.reshape(np.fromfile(f, dtype=np.float32), (3, -1)).T

        assert not array is None, f"Unknown variable type in {filename}"

        return time, array  # to reverse .ravel(order='F'), see the writer


#
#   ENSI Writer
#
# ==============================================================


# ==============================================================
#
#   Other functions
#
def save_geometry_and_get_iterp_weights(
    casename: str, coarse_mesh_vtp_filename: str, out_casename: str = None
):
    """
    Arguments:
        casename : str
            Path and filename to the ensi.case
        coarse_mesh_vtp_filename : str
            Path and filename to the mesh that will be used to save the result
        out_casename : str
            Path and filename (ensi.case) where to save the problem. If None, do not save anything

    Return: dictionary:
        {
            'weights'  : Nx3 ndarray of floats
            'point_ids': Nx3 ndarray of ints
            'variables': set of variables in the case
        }
    """

    print(f"Reading {casename}")
    assert os.path.isfile(casename)
    rdr = vtk.vtkEnSightGoldBinaryReader()
    rdr.SetCaseFileName(casename)
    rdr.ReadAllVariablesOn()
    rdr.Update()

    print(f"Reading coarse mesh {coarse_mesh_vtp_filename}")
    assert os.path.isfile(
        coarse_mesh_vtp_filename
    ), f"Missing {coarse_mesh_vtp_filename}"
    rdr_coarse = vtk.vtkXMLPolyDataReader()
    rdr_coarse.SetFileName(coarse_mesh_vtp_filename)
    rdr_coarse.Update()

    # list variables
    varnames = set()
    for i in range(rdr.GetNumberOfVariables()):
        varnames.add(rdr.GetDescription(i))

    assert len(varnames) > 0, "No variables to save"

    # ---------------------------------------------------------
    #
    #   precalculate interpolation weights
    #

    # extract the points of the exmedi mesh for interpolation later on
    print("Loading geometry")
    rdr.SetTimeValue(rdr.GetTimeSets().GetItem(0).GetTuple1(0))
    rdr.Update()

    print("Precalcualting interpolation weights")
    interp_weights, interp_ptids = calculate_interpolation_weights(
        rdr.GetOutput().GetBlock(0), rdr_coarse.GetOutput()
    )

    if not out_casename is None:
        write_geometry(
            surf_mesh=rdr_coarse.GetOutput(),
            path=os.path.dirname(out_casename),
            casename=get_case_name(os.path.basename(out_casename)),
        )

    return {"weights": interp_weights, "point_ids": interp_ptids, "variables": varnames}


def process_one_variable(
    filenamein, filenameout, interp_weights, interp_ptids
) -> (float, int):
    """Return time and dimensionality of the variable to decide for ensi if it's a scalar or vector"""

    time, array = read_variable_pernode(filenamein)
    data = interpolate_array(array, interp_weights, interp_ptids)
    write_variable_pernode(data, filenameout)

    ncomponents = 3 if data.ndim > 1 else 1

    return time, ncomponents


def process_filelist(files, interp_weights, interp_ptids):
    """Files need to be a list of fiiles only on master"""

    files2scatter = None

    if my_rank == 0:
        print("Scattering exm variables")
        files2scatter = [tuple(files[i::num_procs]) for i in range(num_procs)]

    if num_procs > 1:
        files2scatter = comm.scatter(files2scatter, root=0)
    else:
        files2scatter = files2scatter[0]

    results = []

    if my_rank == 0:
        bar = ProgressBar(maxval=len(files2scatter))

    kk = []
    for i, filename in enumerate(files2scatter):
        if my_rank == 0:
            bar.update(i)

        if not filename is None:
            (prefix, varname, fileno) = re.match(
                "^(.+)\.ensi\.(.+)-(\d+)$", filename
            ).groups()
            if not varname in variables_to_ignore:
                time, ndim = process_one_variable(
                    filename,
                    os.path.join(out_path, f"{out_casename}.ensi.{varname}-{fileno}"),
                    interp_weights,
                    interp_ptids,
                )
                if (np.round(time, time_ndigits) >= time_interval[0]) and (
                    np.round(time, time_ndigits) <= time_interval[1]
                ):
                    results.append(
                        (varname, int(fileno), np.round(time, time_ndigits), ndim)
                    )

        else:
            results.append((None, None, None, None))

    if my_rank == 0:
        print("Gathering exm variables")

    if num_procs > 1:
        results = comm.gather(results, root=0)
    else:
        results = [results]

    return results


def variables_equal_timesets(df_sorted, var1, var2):
    """sorted_df -- must be solrted by time or filenumber"""

    var1_df = df_sorted[df_sorted["Variable"] == var1]
    var2_df = df_sorted[df_sorted["Variable"] == var2]

    equal = False

    print(f"Comparing {var1} and {var2}")
    if var1_df.shape[0] == var2_df.shape[0]:  # definitely not the same timeset
        timediff = (
            np.abs(var1_df["Time"].to_numpy() - var2_df["Time"].to_numpy())
            < time_check_tolerance
        ).all()
        fileno_diff = (
            np.abs(var1_df["FileNo"].to_numpy() - var2_df["FileNo"].to_numpy()) == 0
        ).all()
        print(f"Comparing times for {var1} and {var2}")
        print(f"Times are equal? {timediff}")
        print(f"Filenumbers are equal? {fileno_diff}")

        equal = timediff & fileno_diff

    return equal


def identify_timesets(df):
    """Return dict {varname : timeset,...}"""

    df_sorted = df.sort_values(by="FileNo", ascending=True, axis=0)

    # need to be careful, since the timesteps will be different for solidz and exmedi
    variables = df_sorted["Variable"].unique()

    variable_timeset = {}

    # start with one
    for varname in variables:
        variable_timeset[varname] = None

    variable_timeset[variables[0]] = 1

    timeset = 2
    # go over all pairs of variables
    for i in range(len(variables)):
        # find a similar variable with a set timeset
        for j in range(len(variables)):
            if i != j:
                if not (variable_timeset[variables[j]] is None):
                    if variables_equal_timesets(df_sorted, variables[i], variables[j]):
                        variable_timeset[variables[i]] = variable_timeset[variables[j]]

        # if not found, assign a new timeset
        if variable_timeset[variables[i]] is None:
            variable_timeset[variables[i]] = timeset
            timeset += 1

    return variable_timeset


#
#   Other functions
#
# ==============================================================


# ---------------------------------------------------------
#
#   Main

results = []
for casename in casenames:
    if my_rank == 0:
        print("===========================================================\n")
        print(f"    Processing {casename}")
        print("\n===========================================================")

    interp_weights = None
    interp_ptids = None

    files = []

    if my_rank == 0:
        os.makedirs(out_path, exist_ok=True)

        data4saving = save_geometry_and_get_iterp_weights(
            casename, coarse_mesh, os.path.join(out_path, f"{out_casename}.ensi.case")
        )

        print("Create a list of variable - timestep pairs to parallelize")
        for varname in data4saving["variables"]:
            path = os.path.dirname(casename)
            file_prefix = get_case_name(os.path.basename(casename))
            files.extend(
                glob.glob(os.path.join(path, f"{file_prefix}.ensi.{varname}-*"))
            )

        interp_weights = data4saving["weights"]
        interp_ptids = data4saving["point_ids"]

    if num_procs > 1:
        interp_weights = comm.bcast(interp_weights, root=0)
        interp_ptids = comm.bcast(interp_ptids, root=0)

    # =============================================================
    #
    #   Append arrays for scattering
    #
    if my_rank == 0:
        while len(files) < num_procs:
            files.append(None)

    # =============================================================
    #
    #   Scatter variables
    #
    rrr = process_filelist(files, interp_weights, interp_ptids)
    if my_rank == 0:
        results.extend(rrr)
    #
    #   Scatter variables
    #
    # =============================================================


# =============================================================
#
# process the results and write ensi.case
#
if my_rank == 0:
    print("===========================================================\n")
    print(f"    Writing ensi case")
    print("\n===========================================================")

    results_df = pd.DataFrame(
        itertools.chain(*results), columns=["Variable", "FileNo", "Time", "ndim"]
    )
    results = None  # deallocate

    results_df.dropna(axis=0, inplace=True)

    print("Identify timesets")
    timesets = identify_timesets(results_df)
    print(timesets)

    print(f"Writing casefile {out_casename} to {out_path}")
    write_ensi_case(results_df, timesets)
