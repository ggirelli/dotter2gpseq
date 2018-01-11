dotter2gpseq
=============

A collection of scripts that can be used in combination with pygpseq to perform advanced operations.

Each script needs to be able to load the `pygpseq` library. Either install it systemwide or place a link to the `pygpseq` folder in the same directory as the script (e.g., `ln -s REPOFOLDER/pygpseq ./pygpseq`).

### dotter2gpseq.py

**dotter2gpseq.py** requires the table output of DOTTER, with the FISH dots coordinates in the `x`, `y` and `z` columns, and the NM address in the `File` column. Then, prepare a folder containing the DNA staining channels (deconvolved) with file name in DOTTER notation (e.g., `dapi_001_cmle.tif`).

The script will add 7 columns: `angle`, `Allele`, `cellID`, `lamin_dist`, `lamin_dist_norm`, `centr_dist` and `centr_dist_norm`.

The `Allele` column contains allele labeling:

- No value: dot outside of cells
- -1: more than 2 dots per cell/channel
- 0: less than 2 dots per cell/channel
- 1: more central dot
- 2: more peripheral dot

The `angle` column contains the angle between each allele pair and the nucleus center of mass.

Lamina distance and center distance are calculated using an anisotropic euclidean transform based on the specified aspect (`-a`). Lamina is defined as the first black (background) pixels outside of the nucleus. Center is defined as the nuclear pixels most distant from the lamina. Previous versions of the script extrapolated the center distance from the normalized lamin distance instead.

Distances from center and lamina are normalized on their sum, for each dot. Previous versions of the script normalized over the maximum radius of the nucleus, *assuming a spherical shape*.

* Also, an option (`-a Z Y X`) is available to specify the voxel aspect ratio.
* Use `-t` to specify the number of threads to be used for paralelization.
* Use `--dilate` to specify the number of dilation operations to perform. The dilations are performed nucleus-wise, while the nuclear masked is saved for a general dilation (for simplicity).
* Use `./dotter2gpseq.py -h` for more details.

### dotter2gpseq_merge.R

Useful to merge the output generated by multiple runs of `dotter2gpseq.py`, possibly on different cell lines, etc...

```
usage: merge_data.R [--] [--help] [--opts OPTS] [--meta META] [--indir INDIR] [--outdir OUTDIR] [--aspect ASPECT] 

Description: merge dotter2gpseq.py output, add dataset and
cell type information. The software looks for dotter2gpseq.py output in
subfolders of the specified input directories. These subfolders must be named as
dataset_series, with series being in XXX format with leading zeros.

Example 1: output in current directory.
./merge_data.R -m meta_HAP1.tsv meta_IMR90.tsv -i HAP1/dots_auto IMR90/dots_auto

Example 2: output to "/home/user/out" directory.
./merge_data.R -m meta_HAP1.tsv meta_IMR90.tsv -i HAP1/dots_auto IMR90/dots_auto
    -o /home/user/out


flags:
  -h, --help            show this help message and exit

optional arguments:
  -x, --opts OPTS           RDS file containing argument values
  -m, --meta META           List of metadata tables. Needed columns: dataset, series, cell_line, set_label, probe_label.
  -i, --indir INDIR         List of input folders, same order as metadata.
  -o, --outdir OUTDIR           Output folder, created if missing. Default to current one. [default: .]
  -a, --aspect ASPECT           Physical size of Z, Y and X voxel sides. Default: 300.0 130.0 130.0 [default: (300,130,130)]
```
