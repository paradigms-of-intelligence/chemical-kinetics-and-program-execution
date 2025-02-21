#!/bin/bash

# Copyright 2025 Google LLC
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     https://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# Simplistic bare-bones build script.
# More sophisticated settings would use some Makefile or BUILD rule -
# for this work, we want to do without any such additional tool
# dependency and show the basic steps that need to be done.

gambit_so=libgambit.so

if [ ! -r "${gambit_so}" ] ; then
    echo "Library ${gambit_so} not found in current directory. Trying to use 'locate'."
    candidate=$(locate libgambit.so |
                    while read -r f; do
                        if [ ! -L "${f}" ] ; then echo -n "$f" ; break ; fi
                    done)
    read -r -p "Symlink ${candidate} -> ${gambit_so}? (y/N) " yn
    case "${yn}" in
        y|Y)
            ln -s "${candidate}" "${gambit_so}"
            ;;
        *)
            ;;
    esac
fi


echo "*** (re)building tapes_py_interface.so ***"
rm 2>&1 >/dev/null -f tapes_py_interface{.so,.scm~,.c,_.c} || /bin/true
gsc -link tapes_py_interface
gcc -shared -fPIC -O2 -o tapes_py_interface.so \
    tapes_py_interface_c_glue.c tapes_py_interface.c tapes_py_interface_.c \
    -L. -lgambit
