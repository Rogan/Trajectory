@echo off



rem Build gesop_model.dll

rem

rem The first calling argument is expected to be GESOP_HOME

rem This is the Version for a mixed pure Fortran model

rem

rem Written by Denis Fischer, 27. June 2001

rem

rem Copyright (c) 2001 by IFR

rem $Header: /examples/GESOP_Examples/Bryson-Ho.aps/C_Example.gtp/model_dir/dll/make_model.bat 3     10/19/01 6:15p Peter $



if not exist lib mkdir lib
if ERRORLEVEL 1 goto error_occurred
cd lib

if "%1"=="" goto no_home_dir



rem For convenience, define environment variable

set GESOP_C=%1\src\Gesop_Model\C
set GESOP_ADA=%1\src\Ada95


rem Compiler options

set COMPILER_OPT=/I"%GESOP_C%" /nologo /I..\src\include /MT /W3 /EHsc /O2 /D"WIN32" /D"NDEBUG" /D"_WINDOWS" /D"_MBCS" /D"_USRDLL" /D"C_GESOP_DLL_INTERFACE_EXPORTS" /D"MAKE_A_DLL" /Fp"c_gesop_dll_interface.pch" /FD /c 



rem Linker options

set LINKER_OPT=kernel32.lib /nologo /dll /incremental:no /pdb:"c_gesop_dll_interface.pdb" /machine:I386 /out:"gesop_model.dll" /implib:"c_gesop_dll_interface.lib" 



rem Build C Library




rem check, if source files have been updated since gesop_model.dll was built

dir /S /O:-D /B /W ..\src\*.* > check_files.txt

%1\bin\Check_Date check_files.txt gesop_model.dll

if not ERRORLEVEL 1 goto end_exit



rem Location of gesop_main.lib

set GESOP_MAIN_LIB=%1\lib


rem Execute compiler

cl %COMPILER_OPT% %GESOP_C%\gesop_model.c

if ERRORLEVEL 1 goto error_occured

cl %COMPILER_OPT% ..\src\*.c ..\src\src\*.c
if ERRORLEVEL 1 goto error_occured

if exist ..\src\*.cpp cl %COMPILER_OPT% ..\src\*.cpp

if ERRORLEVEL 1 goto error_occured



rem Execute linker

link *.obj %GESOP_MAIN_LIB%\gesop_main.lib ..\src\cspice\lib\*.lib %LINKER_OPT%

if ERRORLEVEL 1 goto error_occured



goto end_exit



:pause_exit

	pause

	exit



:end_exit

	exit



:no_home_dir

	echo Gesop homedirectory not specified

	echo No gesop_home > DLL_Error.txt

	pause

	exit



:error_occured

	echo ERROR #%ERRORLEVEL% OCCURED

	echo %ERRORLEVEL% > DLL_Error.txt

	pause

	exit
