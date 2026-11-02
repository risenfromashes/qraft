@echo off
rem Checks an unpacked QRAFT release for Windows: the built-in selftest, then
rem the example dataset against the tree this release infers from it.
setlocal
set "here=%~dp0"
set "scratch=%TEMP%"
if "%scratch%"=="" set "scratch=%here%"
set "out=%scratch%\qraft-install-check.tre"
"%here%qraft.exe" --version || goto failed
"%here%qraft.exe" selftest || goto failed
"%here%qraft.exe" -i "%here%examples\genes-100taxa.tre" -o "%out%" -q || goto failed
fc /b "%out%" "%here%examples\expected-100taxa.tre" >nul || goto differs
del "%out%"
echo example: same tree as expected
echo install check passed
exit /b 0
:differs
echo example: the tree differs from examples\expected-100taxa.tre 1>&2
:failed
echo install check FAILED 1>&2
exit /b 1
