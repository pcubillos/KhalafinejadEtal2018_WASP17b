# From the directory where this file is located, execute (this will make things easier):
topdir=`pwd`

# Clone (download) the PB code:
git clone --recursive https://github.com/pcubillos/pyratbay
# Set the version used for this article:
cd $topdir/pyratbay
git checkout fb03915
cd $topdir/pyratbay/modules/MCcubed
git checkout e978f64
cd $topdir/pyratbay/modules/pytips
git checkout f181dc8
cd $topdir/pyratbay/modules/TEA
git checkout 324394c


# Patch pyratbay code for high-res data:
cp $topdir/code/patch/driver.py   $topdir/pyratbay/pyratbay/pbay/driver.py
cp $topdir/code/patch/pyratfit.py $topdir/pyratbay/pyratbay/pbay/pyratfit.py
cp $topdir/code/patch/argum.py    $topdir/pyratbay/pyratbay/pyrat/argum.py
cp $topdir/code/patch/objects.py  $topdir/pyratbay/pyratbay/pyrat/objects.py

# Compile the PB code:
cd $topdir/pyratbay
make

# High-resolution atmospheric runs:
cd $topdir/run/
python $topdir/pyratbay/pbay.py -c upper_atm_solar.cfg

# Retrieval (fixed Na abundance):
cd $topdir/run/
python $topdir/pyratbay/pbay.py -c mcmc_WASP17b.cfg
# Retrieval (free Na abundance):
python $topdir/pyratbay/pbay.py -c mcmc_WASP17b_abundance.cfg

# Best-fit spectrum plot:
cd $topdir/run/
python $topdir/code/fig_hires.py

# Retrieval plots:
cd $topdir/run/
python $topdir/code/fig_post.py
