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
            fp, pathname, description = imp.find_module('_vbwCSC', [dirname(__file__)])
        except ImportError:
            import _vbwCSC
            return _vbwCSC
        if fp is not None:
            try:
                _mod = imp.load_module('_vbwCSC', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _vbwCSC = swig_import_helper()
    del swig_import_helper
else:
    import _vbwCSC
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
    __swig_setmethods__["size"] = _vbwCSC.block_size_set
    __swig_getmethods__["size"] = _vbwCSC.block_size_get
    if _newclass:size = _swig_property(_vbwCSC.block_size_get, _vbwCSC.block_size_set)
    __swig_setmethods__["alphas"] = _vbwCSC.block_alphas_set
    __swig_getmethods__["alphas"] = _vbwCSC.block_alphas_get
    if _newclass:alphas = _swig_property(_vbwCSC.block_alphas_get, _vbwCSC.block_alphas_set)
    __swig_setmethods__["saxsExpPtr"] = _vbwCSC.block_saxsExpPtr_set
    __swig_getmethods__["saxsExpPtr"] = _vbwCSC.block_saxsExpPtr_get
    if _newclass:saxsExpPtr = _swig_property(_vbwCSC.block_saxsExpPtr_get, _vbwCSC.block_saxsExpPtr_set)
    __swig_setmethods__["saxsErrPtr"] = _vbwCSC.block_saxsErrPtr_set
    __swig_getmethods__["saxsErrPtr"] = _vbwCSC.block_saxsErrPtr_get
    if _newclass:saxsErrPtr = _swig_property(_vbwCSC.block_saxsErrPtr_get, _vbwCSC.block_saxsErrPtr_set)
    __swig_setmethods__["saxsEnsPtr"] = _vbwCSC.block_saxsEnsPtr_set
    __swig_getmethods__["saxsEnsPtr"] = _vbwCSC.block_saxsEnsPtr_get
    if _newclass:saxsEnsPtr = _swig_property(_vbwCSC.block_saxsEnsPtr_get, _vbwCSC.block_saxsEnsPtr_set)
    __swig_setmethods__["saxsPrePtr"] = _vbwCSC.block_saxsPrePtr_set
    __swig_getmethods__["saxsPrePtr"] = _vbwCSC.block_saxsPrePtr_get
    if _newclass:saxsPrePtr = _swig_property(_vbwCSC.block_saxsPrePtr_get, _vbwCSC.block_saxsPrePtr_set)
    __swig_setmethods__["saxsMixPtr"] = _vbwCSC.block_saxsMixPtr_set
    __swig_getmethods__["saxsMixPtr"] = _vbwCSC.block_saxsMixPtr_get
    if _newclass:saxsMixPtr = _swig_property(_vbwCSC.block_saxsMixPtr_get, _vbwCSC.block_saxsMixPtr_set)
    __swig_setmethods__["saxsScale"] = _vbwCSC.block_saxsScale_set
    __swig_getmethods__["saxsScale"] = _vbwCSC.block_saxsScale_get
    if _newclass:saxsScale = _swig_property(_vbwCSC.block_saxsScale_get, _vbwCSC.block_saxsScale_set)
    __swig_setmethods__["csExpPtr"] = _vbwCSC.block_csExpPtr_set
    __swig_getmethods__["csExpPtr"] = _vbwCSC.block_csExpPtr_get
    if _newclass:csExpPtr = _swig_property(_vbwCSC.block_csExpPtr_get, _vbwCSC.block_csExpPtr_set)
    __swig_setmethods__["csErrPtr"] = _vbwCSC.block_csErrPtr_set
    __swig_getmethods__["csErrPtr"] = _vbwCSC.block_csErrPtr_get
    if _newclass:csErrPtr = _swig_property(_vbwCSC.block_csErrPtr_get, _vbwCSC.block_csErrPtr_set)
    __swig_setmethods__["csRmsPtr"] = _vbwCSC.block_csRmsPtr_set
    __swig_getmethods__["csRmsPtr"] = _vbwCSC.block_csRmsPtr_get
    if _newclass:csRmsPtr = _swig_property(_vbwCSC.block_csRmsPtr_get, _vbwCSC.block_csRmsPtr_set)
    __swig_setmethods__["csEnsPtr"] = _vbwCSC.block_csEnsPtr_set
    __swig_getmethods__["csEnsPtr"] = _vbwCSC.block_csEnsPtr_get
    if _newclass:csEnsPtr = _swig_property(_vbwCSC.block_csEnsPtr_get, _vbwCSC.block_csEnsPtr_set)
    __swig_setmethods__["csPrePtr"] = _vbwCSC.block_csPrePtr_set
    __swig_getmethods__["csPrePtr"] = _vbwCSC.block_csPrePtr_get
    if _newclass:csPrePtr = _swig_property(_vbwCSC.block_csPrePtr_get, _vbwCSC.block_csPrePtr_set)
    __swig_setmethods__["csMixPtr"] = _vbwCSC.block_csMixPtr_set
    __swig_getmethods__["csMixPtr"] = _vbwCSC.block_csMixPtr_get
    if _newclass:csMixPtr = _swig_property(_vbwCSC.block_csMixPtr_get, _vbwCSC.block_csMixPtr_set)
    __swig_setmethods__["csScale"] = _vbwCSC.block_csScale_set
    __swig_getmethods__["csScale"] = _vbwCSC.block_csScale_get
    if _newclass:csScale = _swig_property(_vbwCSC.block_csScale_get, _vbwCSC.block_csScale_set)
    __swig_setmethods__["numberProcs"] = _vbwCSC.block_numberProcs_set
    __swig_getmethods__["numberProcs"] = _vbwCSC.block_numberProcs_get
    if _newclass:numberProcs = _swig_property(_vbwCSC.block_numberProcs_get, _vbwCSC.block_numberProcs_set)
    def __init__(self): 
        this = _vbwCSC.new_block()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _vbwCSC.delete_block
    __del__ = lambda self : None;
block_swigregister = _vbwCSC.block_swigregister
block_swigregister(block)


def block_alloc(*args):
  return _vbwCSC.block_alloc(*args)
block_alloc = _vbwCSC.block_alloc

def block_free(*args):
  return _vbwCSC.block_free(*args)
block_free = _vbwCSC.block_free

def block_copy(*args):
  return _vbwCSC.block_copy(*args)
block_copy = _vbwCSC.block_copy

def block_copy_construct(*args):
  return _vbwCSC.block_copy_construct(*args)
block_copy_construct = _vbwCSC.block_copy_construct

def block_destroy(*args):
  return _vbwCSC.block_destroy(*args)
block_destroy = _vbwCSC.block_destroy

def ientropy(*args):
  return _vbwCSC.ientropy(*args)
ientropy = _vbwCSC.ientropy

def jensen_shannon_div(*args):
  return _vbwCSC.jensen_shannon_div(*args)
jensen_shannon_div = _vbwCSC.jensen_shannon_div

def find_square_root(*args):
  return _vbwCSC.find_square_root(*args)
find_square_root = _vbwCSC.find_square_root

def SaxsScaleMean(*args):
  return _vbwCSC.SaxsScaleMean(*args)
SaxsScaleMean = _vbwCSC.SaxsScaleMean

def SaxsScaleStandardDeviation(*args):
  return _vbwCSC.SaxsScaleStandardDeviation(*args)
SaxsScaleStandardDeviation = _vbwCSC.SaxsScaleStandardDeviation

def L_function(*args):
  return _vbwCSC.L_function(*args)
L_function = _vbwCSC.L_function

def L_distance(*args):
  return _vbwCSC.L_distance(*args)
L_distance = _vbwCSC.L_distance

def L_print(*args):
  return _vbwCSC.L_print(*args)
L_print = _vbwCSC.L_print

def L_take_step(*args):
  return _vbwCSC.L_take_step(*args)
L_take_step = _vbwCSC.L_take_step

def run_vbw(*args):
  return _vbwCSC.run_vbw(*args)
run_vbw = _vbwCSC.run_vbw
# This file is compatible with both classic and new-style classes.


