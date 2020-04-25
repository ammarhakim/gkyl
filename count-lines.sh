wc -l *.cxx *.h App/*.lua App/*/*.lua Basis/*.lua Comm/*.lua DataStruct/*.lua Eq/*.lua Grid/*.lua Io/*.lua Lib/*.lua Regression/*.lua Updater/*.lua Proto/*/*.lua Unit/*.lua Tool/*.lua */*.h */*.cpp

# total hand-written loc 
echo "Total LOC, including unit tests"
wc -l *.cxx *.h App/*.lua App/*/*.lua Basis/*.lua Comm/*.lua DataStruct/*.lua Eq/*.lua Grid/*.lua Io/*.lua Lib/*.lua Regression/*.lua Updater/*.lua Proto/*/*.lua Unit/*.lua Tool/*.lua */*.h */*.cpp | grep total

echo "Unit tests LOC"
wc -l Unit/*.lua | grep total

echo "Sqlite3 LOC"
wc -l sqlite3/*.h
