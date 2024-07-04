# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

name = "py-plantFATE"
__version__ = "0.0.1"

package_dir = {name: "pybinds"}
depends = ["./inst/include/plantfate_patch.h",
           "./inst/include/pspm_interface.h",
           "./src/plantfate_patch.cpp",
           "./build/plantfate_patch.o",
           "./src/pspm_interface.cpp",
           "./build/pspm_interface.o"]
libraries = ["pspm", "stdc++fs"]
library_dirs = ["./external/libpspm/lib", "./lib"]
extra_compile_args = ["-O3",
                      "-pg", 
                      "-Wall", 
                      "-Wextra", 
                      "-DPHYDRO_ANALYTICAL_ONLY", 
                      "-Wno-sign-compare", 
                      "-Wno-unused-variable", 
                      "-Wno-unused-but-set-variable", 
                      "-Wno-float-conversion", 
                      "-Wno-unused-parameter", 
                    #   "-L" + root_dir + "libpspm/lib",
                    #   "-L./lib",
                      "-lpspm",
                      "-libpfate"
                      ]
include_dirs = ['./inst/include', 
                './external/phydro/inst/include', 
                "./external/libpspm/include",
                "./external/flare/include", 
                "./src", "./build"]
# The main interface is through Pybind11Extension.
# * You can add cxx_std=11/14/17, and then build_ext can be removed.
# * You can set include_pybind11=false to add the include directory yourself,
#   say from a submodule.
#
# Note:
#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/python_example/pull/53)

ext_modules = [
    Pybind11Extension("pypfate",
        ["pybinds/pypfate.cpp"],
        # Example: passing in the version to the compiled code
        define_macros = [('VERSION_INFO', __version__)],
        libraries = libraries, 
        library_dirs = library_dirs, 
        extra_compile_args = extra_compile_args,
        include_dirs=include_dirs,
        ),
]

setup(
    name=name,
    version=__version__,
    # author="Sylvain Corlay",
    # author_email="sylvain.corlay@gmail.com",
    # url="https://github.com/pybind/python_example",
    # description="A test project using pybind11",
    # long_description="",
    # packages = "plantFATE",
    # package_dir="pybinds",
    ext_modules=ext_modules,
    # Currently, build_ext only provides an optional "highest supported C++
    # level" feature, but in the future it may provide more features.
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7",
)