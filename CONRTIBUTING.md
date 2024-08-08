
**
TO BE UPDATED
**

### Updating external libraries (Phydro, libpspm, and Flare)

Do not edit the code in external/ folder unless you are familiar with git subtrees and know how to contribute it back to the parent repos. 

#### Adding a new external lib from a specific ref

```
git remote add -f <newlib> https://github.com/xxx/newlib
git subtree add --prefix external/<newlib> <newlib> <ref> --squash
```

#### updating an external library from a specific ref

```
git subtree pull --prefix external/<lib> <lib> <ref> --squash
```

#### *Danger Zone:* Updating the parent library with changes made locally

```
git subtree push --prefix external/<lib> <lib> <ref<
```

#### Updating the python library

If you want to expose any new functions or classes in python you will need to edit the `pybinds/pypfate.cpp` file and either add a new class object or extend an existing class:

```
	py::class_<namespace::class_name>(m, "class_name", py::dynamic_attr())
		.def(py::init<argument_types>()) // constructor
		.def("function_name", &namespace::class_name::function_name)
		.def_readwrite("class_variable", &namespace::class_name::class_variable);
```

Full documentation for the python binding can be found here: https://pybind11.readthedocs.io/en/stable/basics.html.