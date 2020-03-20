# DrosoRT

DrosoRT contains the MATLAB code to create an agent based model of drosophila locomotion with increasing complexity. Current version of this code assumes that the environment to model is an circular arena 4 cm in radius with an concentric stimulus zone of some known radius.

## Software
This code was run on MATLAB 2019b. However, it should be compatable with MATLAB 2017 and higher.

## Data location

Please download the dataset from [figshare](https://doi.org/10.6084/m9.figshare.11356952.v1)

If you are running the sharp turn analysis portion of the code, please download the Katsov dataset from [Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.854j2)

## Usage

Check to make sure that the data is located in a Data folder under the Data Full folder.

To run all of the code, simply run: `Run_Pipeline`

The code creates postscript analysis files. To convert postscript to pdf, you can use adobe or [ghostscript](https://www.ghostscript.com/)
