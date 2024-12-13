#include "python.h"

#include <Python.h>

#define END Py_Finalize();
#define START Py_Initialize();

namespace python {
    void _err() {
        if (PyErr_Occurred()) {
            PyObject *type, *value, *traceback;
            PyErr_Fetch(&type, &value, &traceback);
            std::string _value = _PyUnicode_AsString(PyObject_Repr(value));
            PERR(_value);
        }
    }

    bool isNone(PyObject*& obj) { return Py_IsNone(obj); }

    PyObject* loadAttrFromObject(PyObject*& obj, const char* name, int type, bool allow_none) {
        PyObject* attr = PyObject_GetAttrString(obj, name);
        if (isNone(attr) && allow_none) {
            return attr;
        }
        PYCHECK
        switch (type) {       
            case PYINT:
                if (!PyLong_Check(attr)) {
                    UERR(std::string(name) + " is not an int");
                }
                break;
            
            case PYFLOAT:
                if (!PyFloat_Check(attr)) {
                    UERR(std::string(name) + " is not a float");
                }
                break;
            
            case PYSTRING:
                if (!PyUnicode_Check(attr)) {
                    UERR(std::string(name) + " is not a string");
                }
                break;
            
            case PYBOOL:
                if (!PyBool_Check(attr)) {
                    UERR(std::string(name) + " is not an bool");
                }
                break;

            case PYLIST:
                if (!PyList_Check(attr)) {
                    UERR(std::string(name) + " is not an list");
                }
                break;
            
            case PYSPARTINDEX_T:
                if (!PyLong_Check(attr)) {
                    UERR(std::string(name) + " is not an int");
                }
                break;
            
            default:
                UERR("Invalid type provided for: " + std::string(name));
        }
        PYCHECK
        return attr;
    }

    // void convertObjToCString(PyObject* obj, char*& s) {
    //     if (!PyUnicode_Check(obj)) {
    //         UERR("can not be converted to a string, it not a string");
    //     }
    //     s = (char*)_PyUnicode_AsString(obj);
    // }

    void convertObjToSpartaIndex_t(PyObject*& obj, SPARTA_NS::index_t& val) {
        if (!PyLong_Check(obj)) {
            UERR("can not be converted to a int, it is not an int");
        }
        val = PyLong_AsSpartaIndex_t(obj);
    }

    void convertObjToLong(PyObject*& obj, long& val) {
        if (!PyLong_Check(obj)) {
            UERR("can not be converted to a int, it is not an int");
        }
        val = PYLONGTOLONG(obj);
    }

    void convertObjToInt(PyObject*& obj, int& val) {
        if (!PyLong_Check(obj)) {
            UERR("can not be converted to a int, it is not an int");
        }
        val = (int)PyLong_AsLong(obj);
    }

    void convertObjToBool(PyObject*& obj, bool& val) {
        if (!PyBool_Check(obj)) {
            UERR("can not be converted to a bool, it is not a bool");
        }
        val = (bool)PyObject_IsTrue(obj);
    }
    
    void convertObjToDouble(PyObject*& obj, double& val) {
        if (!PyFloat_Check(obj)) {
            UERR("can not be converted to a double, it is not a double");
        }
        val = PyFloat_AsDouble(obj);
    }

    void convertObjectToString(PyObject*& obj, char*& buf) {
        PyObject *temp;
        if (!PyUnicode_Check(obj)) {
            if (!PyObject_HasAttrString(obj, "__str__")) {
                UERR("object has no __str__ method that is required");
            }
            temp = PyObject_Str(obj);
        } else {
            temp = obj;
        }
        buf = (char*)_PyUnicode_AsString(temp);
    }

    void loadAttrFromObjectAndConvert(PyObject*& obj, const char* name, long& val) {
        PyObject* attr = loadAttrFromObject(obj, name, PYINT);
        convertObjToLong(attr, val);
    }

    void loadAttrFromObjectAndConvert(PyObject*& obj, const char* name, SPARTA_NS::index_t& val) {
        PyObject* attr = loadAttrFromObject(obj, name, PYSPARTINDEX_T);
        convertObjToSpartaIndex_t(attr, val);
    }

    void loadAttrFromObjectAndConvert(PyObject*& obj, const char* name, int& val) {
        PyObject* attr = loadAttrFromObject(obj, name, PYINT);
        convertObjToInt(attr, val);
    }

    void loadAttrFromObjectAndConvert(PyObject*& obj, const char* name, double& val) {
        PyObject* attr = loadAttrFromObject(obj, name, PYFLOAT);
        convertObjToDouble(attr, val);
    }

    void loadAttrFromObjectAndConvert(PyObject*& obj, const char* name, char*& val) {
        PyObject* attr = loadAttrFromObject(obj, name, PYSTRING);
        convertObjectToString(attr, val);
    }

    void loadAttrFromObjectAndConvert(PyObject*& obj, const char* name, bool& val) {
        PyObject* attr = loadAttrFromObject(obj, name, PYBOOL);
        convertObjToBool(attr, val);
    }

    void setAttrOfObject(PyObject*& obj, const char* name, long val) {
        PyObject* pyAttrValue = PyLong_FromLong(val);
        if (pyAttrValue == NULL) { PYCHECK }

        int setResult = PyObject_SetAttrString(obj, name, pyAttrValue);
        Py_DECREF(pyAttrValue);
        if (setResult == -1) { PYCHECK }
    }

    void setAttrOfObject(PyObject*& obj, const char* name, SPARTA_NS::index_t val) {
        PyObject* pyAttrValue = PyLong_FromSpartaIndex_t(val);
        if (pyAttrValue == NULL) { PYCHECK }

        int setResult = PyObject_SetAttrString(obj, name, pyAttrValue);
        Py_DECREF(pyAttrValue);
        if (setResult == -1) { PYCHECK }
    }

    void setAttrOfObject(PyObject*& obj, const char* name, bool val) {
        PyObject* pyAttrValue = PyBool_FromLong(val);
        if (pyAttrValue == NULL) { PYCHECK }

        int setResult = PyObject_SetAttrString(obj, name, pyAttrValue);
        Py_DECREF(pyAttrValue);
        if (setResult == -1) { PYCHECK }
    }

    void setAttrOfObject(PyObject*& obj, const char* name, double val) {
        PyObject* pyAttrValue = PyFloat_FromDouble(val);
        if (pyAttrValue == NULL) { PYCHECK }

        int setResult = PyObject_SetAttrString(obj, name, pyAttrValue);
        Py_DECREF(pyAttrValue);
        if (setResult == -1) { PYCHECK }
    }

    void setAttrOfObject(PyObject*& obj, const char* name, char* val) {
        PyObject* pyAttrValue = PyUnicode_FromString(val);
        if (pyAttrValue == NULL) { PYCHECK }

        int setResult = PyObject_SetAttrString(obj, name, pyAttrValue);
        Py_DECREF(pyAttrValue);
        if (setResult == -1) { PYCHECK }
    }

    void setAttrOfObjectToNone(PyObject*& obj, const char* name) {
        int setResult = PyObject_SetAttrString(obj, name, Py_None);
        if (setResult == -1) { PYCHECK }
    }

    handler::handler(std::string _file) {
        // ULOG("Running with Python version: " + std::string(Py_GetVersion()));
        // if (util::_me == 0)
        //     fprintf(util::_screen, "Running with Python Embedded:\n  version:%s\n  path: %ls\n", Py_GetVersion(), Py_GetProgramFullPath());
        
        START
        FILE* file = fopen(_file.c_str(), "r");
        if (file != NULL) {
            PyRun_SimpleFile(file, _file.c_str());
            PYCHECK
            fclose(file);
        } else {
            PERR("Could not open: " + _file);  // Print error if script file couldn't be opened
        }

        this->path  = PyUnicode_FromString("__main__");
        this->_main = PyImport_GetModule(this->path);
        PYCHECK
    }

    handler::~handler() { END }

    PyObject* handler::loadObjectWithSetupFromMain(const char* name) {
        return loadObjectWithSetupFromObject(this->_main, name);
    }

    void handler::setupObject(PyObject*& obj) {
        PYCHECK
        PyObject_CallMethod(obj, "_setup", NULL);
        PYCHECK
    }

    PyObject* handler::loadObjectWithSetupFromObject(PyObject*& obj, const char* name) {
        PyObject* _obj = PyObject_GetAttrString(obj, name);
        this->setupObject(_obj);
        return _obj;
    }
}