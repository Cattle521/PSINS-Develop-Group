<!-- vim-markdown-toc GFM -->

- [Open PSINS document in HTML](#open-psins-document-in-html)
- [CMake Build Setup:](#cmake-build-setup)
  - [Unix/Linux:](#unixlinux)
  - [windows](#windows)

<!-- vim-markdown-toc -->

# Open PSINS document in HTML

URL: [https://yinflying.github.io/](https://yinflying.github.io/)

# CMake Build Setup:

## Unix/Linux:
Under current direcotry:
```sh
$ mkdir ./build
$ cd ./build
$ cmake -DCMAKE_BUILD_TYPE=Debug ../
$ make
$ make test
```
## windows
NOTE: have tested under: win10, cmake-gui 3.13.1, visual studio 2013

1. download cmake from https://cmake.org/download/  and install.
2. open cmake-gui
3. "Browse Source..." select "<DIR>/PSINS-Develop-Group/Source"
4. Create directory "<DIR>/PSINS-Develop-Group/Source/build"
5. "Browse Build..." select "<DIR>/PINS-Develop-Group/Source/build"
6. Click "Configure" and "Generate" to generate the Project file
7. Check the directory "<DIR>/PINS-Developo-Group/Source/build", 
    and open vs project
8. Set the "demo" as "start project"
9. Start Debug
10. copy <DIR>PSINS-Develop-Group/Source/demo/MTI_gaotie.txt to 
    <DIR>PSINS-Develop-Group/Source/build/bin/Debug/, so that demo.exe 
    could test the file
