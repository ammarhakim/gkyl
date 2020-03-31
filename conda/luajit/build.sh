git clone http://luajit.org/git/luajit-2.0.git
cd luajit-2.0
git checkout v2.1
make PREFIX=$PREFIX CFLAGS=-fPIC
make XCFLAGS=-DLUAJIT_ENABLE_GC64 install PREFIX=$PREFIX
ln -sf $PREFIX/bin/luajit-2.1.0-beta3 $PREFIX/bin/lua
