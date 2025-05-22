// Python & Java convenience

%{
#include <sstream>
%}

%define STRINGHELPER(classname)
%extend classname {
    const char * toString() {
        std::ostringstream os;
        os << (*self);
        return os.str().c_str();
    }

    std::string __str__() {
        std::ostringstream os;
        os << (*self);
        return os.str();
    }
};
%enddef // STRINGHELPER

%include "rave.i"

%pragma(java) jniclasscode = %{
    static {
        String libname = System.mapLibraryName( "JavaRave" );

        try {
            try {   
                System.load ( "/home/walten/install/lib/jni/" + libname );
                System.out.println( libname + " successfully loaded!" );
            } catch(UnsatisfiedLinkError e) {
                libname = "libJavaRave.so";
                System.load ( "/home/walten/install/lib/jni/" + libname );
                System.out.println( libname + " successfully loaded!" );
            }
        } catch(SecurityException e) {
            System.out.println( libname + " not loaded!" );
            e.printStackTrace();
        } catch(UnsatisfiedLinkError e) {
            System.out.println( libname + " not loaded!" );
            e.printStackTrace();
        }
    }
%}
