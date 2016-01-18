# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.11
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_vbw', [dirname(__file__)])
        except ImportError:
            import _vbw
            return _vbw
        if fp is not None:
            try:
                _mod = imp.load_module('_vbw', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _vbw = swig_import_helper()
    del swig_import_helper
else:
    import _vbw
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class block(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, block, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, block, name)
    __repr__ = _swig_repr
    __swig_setmethods__["size"] = _vbw.block_size_set
    __swig_getmethods__["size"] = _vbw.block_size_get
    if _newclass:size = _swig_property(_vbw.block_size_get, _vbw.block_size_set)
    __swig_setmethods__["alphas"] = _vbw.block_alphas_set
    __swig_getmethods__["alphas"] = _vbw.block_alphas_get
    if _newclass:alphas = _swig_property(_vbw.block_alphas_get, _vbw.block_alphas_set)
    __swig_setmethods__["saxsExpPtr"] = _vbw.block_saxsExpPtr_set
    __swig_getmethods__["saxsExpPtr"] = _vbw.block_saxsExpPtr_get
    if _newclass:saxsExpPtr = _swig_property(_vbw.block_saxsExpPtr_get, _vbw.block_saxsExpPtr_set)
    __swig_setmethods__["saxsErrPtr"] = _vbw.block_saxsErrPtr_set
    __swig_getmethods__["saxsErrPtr"] = _vbw.block_saxsErrPtr_get
    if _newclass:saxsErrPtr = _swig_property(_vbw.block_saxsErrPtr_get, _vbw.block_saxsErrPtr_set)
    __swig_setmethods__["saxsEnsPtr"] = _vbw.block_saxsEnsPtr_set
    __swig_getmethods__["saxsEnsPtr"] = _vbw.block_saxsEnsPtr_get
    if _newclass:saxsEnsPtr = _swig_property(_vbw.block_saxsEnsPtr_get, _vbw.block_saxsEnsPtr_set)
    __swig_setmethods__["saxsPrePtr"] = _vbw.block_saxsPrePtr_set
    __swig_getmethods__["saxsPrePtr"] = _vbw.block_saxsPrePtr_get
    if _newclass:saxsPrePtr = _swig_property(_vbw.block_saxsPrePtr_get, _vbw.block_saxsPrePtr_set)
    __swig_setmethods__["saxsMixPtr"] = _vbw.block_saxsMixPtr_set
    __swig_getmethods__["saxsMixPtr"] = _vbw.block_saxsMixPtr_get
    if _newclass:saxsMixPtr = _swig_property(_vbw.block_saxsMixPtr_get, _vbw.block_saxsMixPtr_set)
    __swig_setmethods__["saxsScale"] = _vbw.block_saxsScale_set
    __swig_getmethods__["saxsScale"] = _vbw.block_saxsScale_get
    if _newclass:saxsScale = _swig_property(_vbw.block_saxsScale_get, _vbw.block_saxsScale_set)
    __swig_setmethods__["numberProcs"] = _vbw.block_numberProcs_set
    __swig_getmethods__["numberProcs"] = _vbw.block_numberProcs_get
    if _newclass:numberProcs = _swig_property(_vbw.block_numberProcs_get, _vbw.block_numberProcs_set)
    def __init__(self): 
        this = _vbw.new_block()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _vbw.delete_block
    __del__ = lambda self : None;
block_swigregister = _vbw.block_swigregister
block_swigregister(block)


def block_alloc(*args):
  return _vbw.block_alloc(*args)
block_alloc = _vbw.block_alloc

def block_free(*args):
  return _vbw.block_free(*args)
block_free = _vbw.block_free

def block_copy(*args):
  return _vbw.block_copy(*args)
block_copy = _vbw.block_copy

def block_copy_construct(*args):
  return _vbw.block_copy_construct(*args)
block_copy_construct = _vbw.block_copy_construct

def block_destroy(*args):
  return _vbw.block_destroy(*args)
block_destroy = _vbw.block_destroy

def ientropy(*args):
  return _vbw.ientropy(*args)
ientropy = _vbw.ientropy

def jensen_shannon_div(*args):
  return _vbw.jensen_shannon_div(*args)
jensen_shannon_div = _vbw.jensen_shannon_div

def polySolver(*args):
  return _vbw.polySolver(*args)
polySolver = _vbw.polySolver

def find_poly_root(*args):
  return _vbw.find_poly_root(*args)
find_poly_root = _vbw.find_poly_root

def SaxsScaleMean(*args):
  return _vbw.SaxsScaleMean(*args)
SaxsScaleMean = _vbw.SaxsScaleMean

def SaxsScaleStandardDeviation(*args):
  return _vbw.SaxsScaleStandardDeviation(*args)
SaxsScaleStandardDeviation = _vbw.SaxsScaleStandardDeviation

def L_function(*args):
  return _vbw.L_function(*args)
L_function = _vbw.L_function

def L_distance(*args):
  return _vbw.L_distance(*args)
L_distance = _vbw.L_distance

def L_print(*args):
  return _vbw.L_print(*args)
L_print = _vbw.L_print

def L_take_step(*args):
  return _vbw.L_take_step(*args)
L_take_step = _vbw.L_take_step

def run_vbw(*args):
  return _vbw.run_vbw(*args)
run_vbw = _vbw.run_vbw
# This file is compatible with both classic and new-style classes.


