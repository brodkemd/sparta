#ifndef PYTHON_H
#define PYTHON_H

#include <limits.h>
#include <string>

#include "util.h"

#ifdef SPARTA_LONGLONG_TO_LONG
#define PYLONGTOLONG(val) PyLong_AsLong(val)
#define PYLONGTOULONG(val) PyLong_AsUnsignedLong(val)
#define PYLONGFROMLONG(val) PyLong_FromLong(val)
#define PYLONGFROMULONG(val) PyLong_FromUnsignedLong(val)
#else
#define PYLONGTOLONG(val) PyLong_AsLongLong(val)
#define PYLONGTOULONG(val) PyLong_AsUnsignedLongLong(val)
#define PYLONGFROMLONG(val) PyLong_FromLongLong(val)
#define PYLONGFROMULONG(val) PyLong_FromUnsignedLongLong(val)
#endif

#ifdef SPARTA_BIG
#define PyLong_AsSpartaIndex_t(val) PYLONGTOULONG(val)
#define PyLong_FromSpartaIndex_t(val) PYLONGFROMULONG(val)
#endif

// for problems with more than 2B grid cells
// 32-bit smallint, 64-bit bigint, 64-bit cellint
#ifdef SPARTA_BIGBIG
#define PyLong_AsSpartaIndex_t(val) PYLONGTOULONG(val)
#define PyLong_FromSpartaIndex_t(val) PYLONGFROMULONG(val)
#endif

// for machines that do not support 64-bit ints
// 32-bit smallint and bigint and cellint

#ifdef SPARTA_SMALL
#define PyLong_AsSpartaIndex_t(val) _PyLong_AsInt(val)
#define PyLong_FromSpartaIndex_t(val) PyLong_FromLong(val)
#endif

#define PERR(_msg) END UERR(_msg);
#define PYCHECK python::_err();

// dumby stuff
struct _object; /* Incomplete forward declaration */
typedef struct _object PyObject;


namespace elmer{
    class Section;
}

namespace python {
    enum{PYINT,PYFLOAT,PYSTRING,PYBOOL,PYLIST,PYSPARTINDEX_T};

    void _err();

    bool isNone(PyObject*& obj);

    void convertObjToLong(     PyObject*& obj, long& val);
    void convertObjToInt(      PyObject*& obj, int& val);
    void convertObjToBool(     PyObject*& obj, bool& val);
    void convertObjToDouble(   PyObject*& obj, double& val);
    void convertObjectToString(PyObject*& obj, char*& buf);
    void convertObjToSpartaIndex_t(PyObject*& obj, SPARTA_NS::index_t& val);

    PyObject* loadAttrFromObject(PyObject*& obj, const char* name, int type, bool allow_none = false);
    void loadAttrFromObjectAndConvert(PyObject*& obj, const char* name, long&   val);
    void loadAttrFromObjectAndConvert(PyObject*& obj, const char* name, int&    val);
    void loadAttrFromObjectAndConvert(PyObject*& obj, const char* name, double& val);
    void loadAttrFromObjectAndConvert(PyObject*& obj, const char* name, char*&  val);
    void loadAttrFromObjectAndConvert(PyObject*& obj, const char* name, bool&   val);
    void loadAttrFromObjectAndConvert(PyObject*& obj, const char* name, SPARTA_NS::index_t&   val);

    void setAttrOfObject(PyObject*& obj, const char* name, long val);
    void setAttrOfObject(PyObject*& obj, const char* name, bool val);
    void setAttrOfObject(PyObject*& obj, const char* name, double val);
    void setAttrOfObject(PyObject*& obj, const char* name, char* val);
    void setAttrOfObject(PyObject*& obj, const char* name, SPARTA_NS::index_t val);
    void setAttrOfObjectToNone(PyObject* obj, const char* name);

    class handler {
        private:
            PyObject *_main, *path;

        public:
            handler() {}
            handler(std::string _file);
            PyObject* loadObjectWithSetupFromMain(const char* name);
            void setupObject(PyObject*& obj);
            PyObject* loadObjectWithSetupFromObject(PyObject*& obj, const char* name);
            ~handler();       
    };
}

#endif