# Simplee instructions to build LuaJIT 2.1
git clone http://luajit.org/git/luajit-2.0.git
cd luajit-2.0
git checkout v2.1
make
make install PREFIX=$HOME/gkylsoft/luajit
