# From the directory where this file is located, execute (this will make things easier):
topdir=`pwd`

# Clone (download) the PB code:
git clone --recursive https://github.com/pcubillos/pyratbay
cd $topdir/pyratbay
git checkout fb03915

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

# Retrieval:
cd $topdir/run/
python $topdir/pyratbay/pbay.py -c mcmcm_WASP17b.cfg

# Best-fit spectrum plot:
cd $topdir/run/
python $topdir/code/fig_hires.py

# Retrieval plots:
cd $topdir/run/
python $topdir/code/fig_post.py
