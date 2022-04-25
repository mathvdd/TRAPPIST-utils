# TRAPPIST-utils
trapccd, subsets and login.cl are to be put in ~/iraf
Need to change path in 2 shell commands in prograp3
IRAF installation with conda on Fedora (community package did not work with the TRAPPIST pipeline):

https://faculty1.coloradocollege.edu/~sburns/courses/18-19/pc362/Anaconda_IRAF_install.html

(may be necessary: sudo dnf install glib2.i686 libgcc.i686 libstdc++.i686 gtk3.i686)
conda config --add channels http://ssb.stsci.edu/astroconda
conda create -n iraf27 python=2.7 iraf-all pyraf-all stsci
conda activate iraf27
mkdir iraf
cd iraf
mkiraf

can modify iraf/login.cl to include shortcut to TRAPDAT and scripts (see exemple)

progtrap2: trapccd (and subsets?) needs to be in the main repertoire and needs to be launched from there
