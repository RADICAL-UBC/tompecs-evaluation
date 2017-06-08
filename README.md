# mc

To install ransampl, you need gsl. To install gsl on Ubuntu, install the following packages:

gsl-bin: GNU Scientific Library (GSL) -- binary package
libgsl0-dbg: GNU Scientific Library (GSL) -- debug symbols package
libgsl0-dev: GNU Scientific Library (GSL) -- development package
libgsl0ldbl: GNU Scientific Library (GSL) -- library package

Command: sudo apt-get install gsl-bin libgsl0-dbg libgsl0-dev libgsl0ldbl


After that, simply run

      ./configure
      make
      sudo make install

If make install gives you trouble about /libtool not being found, just
copy "libtools" script from ransampl's directory to /

Note: You also need autoconf and automake to build ransampl


Gurobi:
======
Gurobi might not find libgurobi70.so when running the simulation binary

You need to update LD_LIBRARY_PATH with /opt/gurobi702/linux64/lib.
Under Ubuntu, the right way to do it is to add a custom .conf file to /etc/ld.so.conf.d; for example,

      sudo emacs /etc/ld.so.conf.d/baderLibs.conf
      
Inside this file you are supposed to write the complete path to the directory that contains all the libraries that you wish to add to the system, for example,

       /opt/gurobi702/linux64/lib
       
Remember to add only the path to the dir, not the full path to the file; all
the libs inside that path will be automatically indexed.

Now save and run

     sudo ldconfig
     
to update the system with this libs.

Gurobi Lisense: After creating the gurobi license file with the grbgetkey tool, a gurobi.lic file will be generated. Store this file in a location that is in the PATH variable. To update the PATH variable on Ubuntu, modify ~/.bashrc and include the line

       export PATH=$PATH:../

for instance to add ../ to the path. Log out the shell and login again for the changes to take effect.

