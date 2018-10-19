unzip htslib-1.7.zip
cd htslib-1.7
autoheader
autoconf
./configure --disable-lmza --prefix=`pwd`
make
make install
