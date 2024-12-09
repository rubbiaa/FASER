FASER_BASE_PATH=$1

if [ -z $FASER_BASE_PATH ]; then
	echo "Use as ./generateLib.sh /full/path/to/Faser/main/dir"
	exit 1
fi

rootcint -f $FASER_BASE_PATH/Python_io/lib/ClassesDict.cxx -c $FASER_BASE_PATH/CoreUtils/*.hh $FASER_BASE_PATH/Python_io/lib/LinkDef.h
g++ -o $FASER_BASE_PATH/Python_io/lib/ClassesDict.so -shared -fPIC `root-config --cflags --libs` $FASER_BASE_PATH/CoreUtils/*.cc $FASER_BASE_PATH/Python_io/lib/ClassesDict.cxx
