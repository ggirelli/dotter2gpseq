**!!!NOTE!!! This repo is being merged with [pygpseq](http://github.com/ggirelli/gpseq-img-py). Check the [pre-merge](https://github.com/ggirelli/dotter2gpseq/tree/pre-merge) branch for the pre-merge code.**

## How to switch to the pre-merge version (branch not maintained anymore)

From within the repository folder run the following:

```
git fetch
git checkout pre-merge
```

I would strongly advise to uninstall `dotter2gpseq` and switch to `pygpseq` instead:

```
rm -rf dotter2gpseq
git clone http://github.com/ggirelli/gpseq-seq-py
cd gpseq-seq-py
sudo -H pip3 install -e .
```