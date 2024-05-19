# Compare the output produced by different runs of the A1 script.
# This output is the "MicrohabitatMatrix*.mat" files and they must
# be in version <7.3 (No HDF5-based format).
# Before running this script, edit the `path1` and `path2` variables,
# to provide the paths to the corresponding output and, if needed,
# edit the `array_name`, `start` and `end` variables too.
import math
from pathlib import Path

import numpy as np
import scipy.io as sio


def get_ndarrays(path1, path2, filename, array_name):
    filepath1 = path1 / filename
    filepath2 = path2 / filename

    f1 = sio.loadmat(filepath1)
    f2 = sio.loadmat(filepath2)

    Microhabitat1 = f1[array_name]
    Microhabitat2 = f2[array_name]

    return Microhabitat1, Microhabitat2


def compare(arr1, arr2):
    # Check size & shape
    if arr1.size != arr2.size:
        raise Exception("Different size")
    if arr1.shape != arr2.shape:
        raise Exception("Different shape")

    # At this point the shapes are the same
    shape = arr1.shape

    all_same = True
    for index in np.ndindex(shape):
        val1 = arr1[index]
        val2 = arr2[index]

        if np.isnan(val1) and np.isnan(val2):
            # Both are NaN
            continue
        else:
            # In math.isclose() NaN is not considered close to any other value, including NaN
            if not math.isclose(arr1[index], arr2[index], rel_tol=1e-9):
                print(f"{index}: {val1} != {val2}")
                all_same = False

    return all_same

def main():
    path1 = Path("")
    path2 = Path("")
    array_name = "Microhabitat"
    start = 1
    end = 40

    for N in range(start, end):
        filename = f"MicrohabitatMatrix{N}.mat"
        Microhabitat1, Microhabitat2 = get_ndarrays(path1, path2, filename, array_name)
        allsame = compare(Microhabitat1, Microhabitat2)

        print(f"{filename} All same" if allsame else f"{filename} NOT SAME!!!")


if __name__ == "__main__":
    main()
