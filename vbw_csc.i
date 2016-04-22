/* vbw.i */
%module vbwCSC
%include "std_string.i"
%{
  #define SWIG_FILE_WITH_INIT
  #include "VBW_csc.hh"
%}

 %include "VBW_csc.hh"
