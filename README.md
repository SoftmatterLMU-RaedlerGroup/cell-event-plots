# Plotting cell death data
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1418463.svg)](https://doi.org/10.5281/zenodo.1418463)

This repository contains the Matlab code for plotting cell death data.

## License
The software in this repository may be used under the terms of the
General Public License version 2, see the `LICENSE` file for details.
However, the file `tight_subplot.m` may be used under the license terms it contains.

## Preparation
### Installation
Clone out this repository and add its base directory to your Matlab path (or enter the directory within Matlab).
You may need Matlab R2017a or newer for running these scripts.

The scripts use the UTF-8 encoding.
If your Matlab installation is configured to use another encoding,
some characters may be displayed incorrectly,
but the basic functionality of the scripts should not be affected.

### Data
Your data needs to be formatted as described below in the section “Data format”.

There is [example data](https://doi.org/10.5281/zenodo.1418377) available.
To use it, download one of the ZIP files and extract it.
Then copy the directories `Raw` and `Fitted` into the base directory of this repository and rename them
to `datadir_raw` and `datadir_state`, respectively.

## Usage
At first, the data have to be loaded. This can be done by calling `loadCorrelationData`.

You will be asked whether to use the default data directories (see above) or custom ones.
If you decide to use custom base directories, you are prompted for a raw and a state base directory.

Then, you are asked to select, from a list ouf found conditions, the conditions you want to load.
Note that the data from all loaded conditions will be evaluated together.
If you want to plot data from different conditions separately, a separate run is needed for each condition.

After collecting the data, which may take some time, you are asked to select the marker combinations to be loaded.
The loaded marker combinations will then be plotted in a script-dependent way.

Before actual plotting, you may want to apply corrections for measurement-dependent time delays.
This can be done by calling `applyTimeDelays` after loading the data.
The measurement-dependent time delays are read from the file `delays.json`.

Now, you can call the scripts for plotting.
They all start with `plotCorrelation`.
In most cases, the figures and tables created by those scripts are automatically saved in the directory `out`
inside this repository.

The following scripts can be used:
* `plotCorrelation_measurement_single` plots, for each marker combination separately, a scatter plot of event times
    in which events are drawn in different colors depending on the corresponding measurements.

* `plotCorrelation_single` plots, for each marker combination separately, a scatter plot of event times.

* `plotCorrelation_cluster_single` plots, for each marker combination separately, a scatter plot of event times.
    If the variable `makeEllipses` exists and is true, an error ellipse will be plotted in the scatter plot.
    If the variable `makeClusters` exists and is true, event time clusters are searched by k-means clustering.

* `plotCorrelation_cluster_meanshift` works similarly to  `plotCorrelation_cluster_single`.
    However, meanshift clustering is performed instead of k-means clustering,
    and a 3D plot is shown for the clustering result
    (right-click on the gray edge of the 3D plot for saving it or changing its appearance).
    Optionally, ellipses can be asymmetric.

* `plotCorrelation_info` prints information about clusters in the command window and in a text file.

* `plotCorrelation_ROS` plots, for each marker combination involving the ROS marker, a scatter plot of the
    slope of the ROS signal and the event time found by the other marker.

* `plotCorrelation_histograms` plots, per marker, a histogram of event times.
    An interactive window appears, allowing selection of histogram properties:
    * Choose whether to normalize the histogram.
    * Choose whether to set the histogram heights manually.
        If yes, you will be prompted to adjust the upper axis limit for each histogram.
    * Choose which data to include in the histograms: all events, only events of cells showing events for both markers,
        only events of cells showing events for only the first or the second marker.
    * Choose the criteria for sub-histograms:
        `None` (no sub-histograms), `Combinations` (a sub-histogram for each marker combination),
        `Measurements` (a sub-histogram for each measurement),
        `Experiments` (a sub-histogram for each measurement-marker combination).

## Data format
You also need the data you want to plot in a defined format.

In general, the information acquired by this software is contained in three types of files.
All of these files are CSV files in the default format written by `csvwrite`, which is: no header,
comma as field delimiter, only numeric data with `.` as decimal separator, no quotation marks.

* The first type of files are the raw files. They contain the raw fluorescence time traces of the cells.

    The first column is the vector of times at which the fluorescence values have been evaluated.
    It is a strictly monotonically increasing series of non-negative values.
    The recommended time unit is hours.

    Each other column contains the fluorescence values of one cell, where the value in a given row
    corresponds to the time in that row.

* The second type of files are the parameter files. They contain the best estimates of the parameters
    of the corresponding model function for fitting the fluorescence time traces.
    
    Each row corresponds to one cell, and each column to a parameter.
    
* The third type of files are the state files. They contain information about the fluorescence time trace
    that was extracted in postprocessing.
    
    Also here, each row corresponds to one cell. The columns correspond to given properties, which are
    (in correct order):
    1. The index of the trace in the file (in all three types of files)
    2. The event time recognized. Non-finite values may be used to indicate that no event was found;
        negative values are currently not used.
    3. The absolute amplitude of the trace (maximum - minimum)
    4. The relative amplitude of the trace (absulute amplitude divided by the largest non-outlier amplitude of the file)
    5. The logarithmic likelihood of the best fit.
    6. The fit type, indicated by a number. The mapping depends on the fitting software and is currently not used.
    7. The slope of the trace at the event

These files are stored in a well-defined directory structure described below.
Since the names of the directories (and files) are evaluated by this software,
the naming convention specified below must be adhered to.
For historical reasons, the raw files are stored separately from the parameter and state files.
However, it should technically be possible to store them together in the same directory structure.

* The base directories of the directory structures will be called “raw base directory” and “state base directory”
    in the following. As mentioned above, the base directories may also be identical.

* The base directories may be specified when prompted by the software.
    Alternatively, a default raw base directory and a default state base directory may be created as directories
    called `datadir_raw` and `datadir_state`, respectively, in this repository.
    Also symlinks pointing to the corresponding directories are accepted.

* Inside the base directories, there must be directories named after the measurements from which the data originates.
    A typical measurement name (and thus directory name) is the date on which the measurement was performed,
    formatted as `YYYYMMDD`.

* Inside the measurement directories, there must be directories named after the conditions that have been measured.
    Typical condition names used in our project are `ctrl`, `NP25`, `NP100` and `sts`.

* Inside the condition directories, there must be directories named after the fluorescent markers evaluated.
    Typical names of fluorescent markers used in our project are:
    `lyso`, `ros`, `tmrm`, `casp`, `psiva`, `pi` and `toto`.
    Some of these names may be used in special ways
    (for example, the script `plotCorrelation_ROS` plots `ros` data in a special way).

* Inside the marker directories, the data files are located.
    Their names are parsed and, thus, must comply with this convention.
    Special care must be taken for the correct number and position of underscores, which serve as separators
    for parts of the file names.

    All file names start with an initial part indicating the measurement and additional information,
    which is, in our case, the position of the sample to which the file corresponds.
    The position information typically consists of two integer numbers separated by an underscore.
    The meaning of the position information is not relevant.
    The position information is at the beginning of the file name, and is followed by an underscore
    and the measurement name.
    
    For the raw files, only the filename extension `.csv` (or, for historical reasons, `.txt`) is appended
    to the initial part.
    
    For the parameter and state files, the string `_ALL_`, the keyword `PARAMS` and `STATE`, respectively,
    an underscore, a (random) identifier typically consisting of four digits, and the extension `.csv`
    are appended to the initial part.
    
    For example, for a measurement performed on 17 September 2015 at a position `1_4`, there may be
    * a raw file called `1_4_20150917.csv`,
    * a parameter file called `1_4_20150917_ALL_PARAMS_1234.csv` and
    * a state file called `1_4_20150917_ALL_STATE_1234.csv`.
    
    The position and measurement information are used to match pairs of markers measured within the same cell.
    Therefore, each position-measurement combination should be present in exactly two marker directories.
    The behaviour for more than two markers with the same position-measurement combination is undefined.

## Exporting event times
A list of all event times saved in `corrT` can be exported by calling `export_times`.
This will create a mat-file of all event times in the `out` directory.

The data can be converted into a CSV file or a pickled pandas table for use with Python
with the python script `convert_matfile.py` in the directory `py`.
It requires Python 3.6 or later as well as numpy, scipy and pandas.
Simply call the script with the path to the mat-file of event times as commandline argument,
and it will create an CSV file and a compressed pickle file in the `out` directory.

## For hackers
This section provides information that is helpful when changing the code or adding functionality.

### Basic data structures
After calling `loadCorrelationData`, your workspace contains some new variables:
* `corrT` is a table holding information about the marker combinations.
* `combT` should not be used; it only contains information also present in `corrT` and is kept for historical reasons.
* `fileT` is a table holding information about the files.
* `restrict_to_conditions` is a string array of the conditions loaded
    (this information is needed for adding the correct condition labels to the plots).
* `out_dir` is the path to which the files are written; by default it is `out` in this repository.
