#ifndef PYTHON_H
#define PYTHON_H

#include <Python.h>
#include <string>

#define END Py_Finalize();
#define START Py_Initialize();
#define PERR(_msg) END UERR(_msg);
#define PYCHECK python::_err();

namespace elmer{
    class Section;
}

namespace python {
    const long noInt = LONG_MIN;
    enum{PYINT,PYFLOAT,PYSTRING,PYBOOL,PYLIST};

    void _err();

    bool isNone(PyObject*& obj);

    void convertObjToLong(     PyObject*& obj, long& val);
    void convertObjToInt(      PyObject*& obj, int& val);
    void convertObjToBool(     PyObject*& obj, bool& val);
    void convertObjToDouble(   PyObject*& obj, double& val);
    void convertObjectToString(PyObject*& obj, char*& buf);

    PyObject* loadAttrFromObject(PyObject*& obj, const char* name, int type, bool allow_none = false);
    void loadAttrFromObjectAndConvert(PyObject*& obj, const char* name, long&   val);
    void loadAttrFromObjectAndConvert(PyObject*& obj, const char* name, int&    val);
    void loadAttrFromObjectAndConvert(PyObject*& obj, const char* name, double& val);
    void loadAttrFromObjectAndConvert(PyObject*& obj, const char* name, char*&  val);
    void loadAttrFromObjectAndConvert(PyObject*& obj, const char* name, bool&   val);

    void setAttrOfObject(PyObject*& obj, const char* name, long val);
    void setAttrOfObject(PyObject*& obj, const char* name, bool val);
    void setAttrOfObject(PyObject*& obj, const char* name, double val);
    void setAttrOfObject(PyObject*& obj, const char* name, char* val);
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
            ~handler() { END }           
    };
}

#endif