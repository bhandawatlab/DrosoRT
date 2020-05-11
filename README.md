# DrosoRT

DrosoRT contains the MATLAB code to create an agent based model of drosophila locomotion with increasing complexity. Current version of this code assumes that the environment to model is an circular arena 4 cm in radius with an concentric stimulus zone of some known radius.

This code accompanies the paper: Tao, L., S. Ozarkar, and V. Bhandawat, "Mechanisms underlying attraction to odors in walking Drosophila. PLOS Computational Biology", 2020. 16(3): p. e1007718.

## Software
This code was run on MATLAB 2019b. However, it should be compatable with MATLAB 2017 and higher.

## Data location

Please download the dataset from [figshare](https://figshare.com/articles/Mechanisms_underlying_attraction_to_odors_in_walking_Drosophila/11356952/2)

If you are running the sharp turn analysis portion of the code, please download the Katsov dataset from [Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.854j2)

## Usage

Check to make sure that the data is located in a Data folder under the Data Full folder.

To run all of the code, simply run: `Run_Pipeline`

The code creates postscript analysis files. To convert postscript to pdf, you can use adobe or [ghostscript](https://www.ghostscript.com/)
