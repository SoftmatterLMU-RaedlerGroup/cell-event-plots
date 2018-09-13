#! /usr/bin/env python3
from datetime import datetime
import inspect
import numpy as np
import os
import pandas as pd
import pickle
import scipy.io as scio
import sys
import time
"""
Read a mat-file of event times and convert it to CSV and pickle.

Call this script with the path of the mat-file as commandline argument
and the converted files are written to the `out` directory.
Examples for reading the CSV and pickle files can be found in the
functions `read_csv` and `read_pickled_pandas`.

Copyright © 2018 Daniel Woschée <daniel.woschee@physik.lmu.de>
Faculty of Physics / Ludwig-Maximilians-Universität München

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 2
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""


def getTimeStamp():
    """Returns a human-readable string representation of the current time"""
    return datetime.now().strftime("%Y-%m-%d–%H%M%S")


def getOutpath(filename='', timestamp=None):
    """Returns (and creates, if necessary) the path to a directory
    called “out” inside the parent directory.
    If `filename` is given, the filename is appended to the output directory.
    A timestamp will be added to the filename if `timestamp != ''`.
    If timestamp is `None`, the current timestamp is used.
    """
    # Create output directory
    try:
        main_dir = os.path.dirname(os.path.abspath(
            inspect.getsourcefile(sys.modules["__main__"])))
        outpath = os.path.abspath(os.path.join(main_dir, '..', 'out'))
    except TypeError:
        outpath = os.path.join(os.getcwd(), 'out')

    if not os.path.isdir(outpath) and not os.path.lexists(outpath):
        os.mkdir(outpath)

    # If requested, build filename
    if len(filename) > 0:
        if timestamp == None:
            timestamp = getTimeStamp()
        outpath = os.path.join(outpath, ((timestamp + '_') if len(timestamp) > 0 else '') + filename)
    return outpath


def to_pickled_pandas(mfilename):
    """Read matfile from path `mfilename` and pickle it as DataFrame"""
    write_pickled_pandas(matimport(mfilename))


def write_pickled_pandas(df):
    """Pickle DataFrame and save it to “out” directory"""
    picklepath = getOutpath("event_times.pickle.xz")
    df.to_pickle(picklepath, compression='xz')
    print(f"Event times pickled to: {picklepath}")


def read_pickled_pandas(path):
    """Read pickled DataFrame with event times"""
    df = pd.read_pickle(path)
    return df


def to_csv(mfilename):
    """Read matfile from path `mfilename` and save it as CSV"""
    write_csv(matimport(mfilename))


def write_csv(df):
    """Save DataFrame to CSV in “out” directory"""
    csvpath = getOutpath("event_times.csv")
    df.to_csv(csvpath, na_rep="NaN")
    print(f"Event times written to: {csvpath}")


def read_csv(path):
    """Read event times from CSV into DataFrame"""
    df = pd.read_csv(path, index_col=0, dtype={
            "id": np.uint32,
            "condition": np.object,
            "measurement": np.object,
            "position": np.object,
            "trace": np.uint16,
            "marker0": np.object,
            "marker1": np.object,
            "t_event0": np.float64,
            "t_event1": np.float64,
            "delay": np.float64,
            })
    return df


def matimport(mfilename):
    """Import event times as `pandas.DataFrame`.
    
    The argument `mfilename` must be a string of a path at which a
    Matlab matfile with a version earlier than 7.3 is located.
    The matfile at this path is read and returned as `pandas.DataFrame`.

    The following variables are expected in the matfile:
    * `id`: an integer array of unique trace identifiers used as
          table index
    * `condition`: a cellstring holding, for each trace, the condition
    * `measurement`: a cellstring holding, for each trace, the measurement
    * `position`: a cellstring, holding, for each trace, the position
    * `trace`: an integer array indicating the (one-based) index of the
          trace in its file
    * `markers`: a cellstring with two columns holding
          the names of the two markers used for the traces in the
          respective rows
    * `t_event`: a double array with two columns holding in each field
          the event time detected for the corresponding trace and marker
    * `delay`: the difference between first and second event. Unless
          both events are detected, the delay is NaN.

    All variables must have as many rows as there are traces. Missing
    variables result in an error. Other variables are ignored.
    """
    print(f"Importing {mfilename} …")
    mat = scio.loadmat(mfilename)

    mat_id = mat['id'].astype(np.uint32)
    mat_trace = mat['trace'].astype(np.uint16)
    mat_t_event = mat['t_event'].astype(np.float64)
    mat_delay = mat['delay'].astype(np.float64)

    mat_condition = mat['condition']
    for i_mc, mc in enumerate(mat_condition.flat):
        mat_condition.flat[i_mc] = str(mc[0])

    mat_measurement = mat['measurement']
    for i_mm, mm in enumerate(mat_measurement.flat):
        mat_measurement.flat[i_mm] = str(mm[0])

    mat_position = mat['position']
    for i_mp, mp in enumerate(mat_position.flat):
        mat_position.flat[i_mp] = str(mp[0])

    mat_markers = mat['markers']
    for i_mm, mm in enumerate(mat_markers.flat):
        mat_markers.flat[i_mm] = str(mm[0])

    # Create DataFrame
    index = pd.Index(data=mat_id.flat, name="id", dtype=np.uint32)
    df = pd.DataFrame(mat_condition, index=index, columns=["condition"])
    df.loc[:, 'measurement'] = mat_measurement
    df.loc[:, 'position'] = mat_position
    df.loc[:, 'trace'] = mat_trace
    df.loc[:, 'marker0'] = mat_markers[:,0]
    df.loc[:, 'marker1'] = mat_markers[:,1]
    df.loc[:, 't_event0'] = mat_t_event[:,0]
    df.loc[:, 't_event1'] = mat_t_event[:,1]
    df.loc[:, 'delay'] = mat_delay

    return df


if __name__ == "__main__":
    df = matimport(sys.argv[1])
    write_pickled_pandas(df)
    write_csv(df)
