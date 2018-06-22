# From the directory where this file is located, execute (this will make things easier):
topdir=`pwd`

# Clone (download) the PB code:
git clone --recursive https://github.com/pcubillos/pyratbay
cd $topdir/pyratbay
git checkout fb03915

# Compile the PB code:
cd $topdir/pyratbay
make

# High-resolution atmospheric runs:
cd $topdir/run/
python $topdir/pyratbay/pbay.py -c upper_atm_solar.cfg

# Retrieval:
cd $topdir/run/
python $topdir/pyratbay/pbay.py -c WASP17b_retrieval.cfg

# Forward model plot:
cd $topdir/run/
python $topdir/fig_hires.py

# Retrieval plots:
cd $topdir/run/
python $topdir/fig_post.py

