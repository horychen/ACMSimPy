@echo off
setlocal

set "SRC=C:\Users\lenovo\Codes\ACMSimPy\course"
set "DST=C:\Users\lenovo\Codes\horychen.github.io\public\motor-html"

echo Copying lecture, homework, and coding project HTML files...
for %%F in ("%SRC%\lecture*.html" "%SRC%\homework*.html" "%SRC%\codingProject*.html") do (
    echo   %%~nxF
    copy /Y "%%F" "%DST%\"
)

echo.
echo Running html-to-dir.bat...
pushd "%DST%"
call "%DST%\html-to-dir.bat"
popd

echo.
echo Copying published landing page...
copy /Y "%SRC%\index-public.html" "%DST%\index.html"

echo.
echo Committing and pushing to GitHub...
pushd "C:\Users\lenovo\Codes\horychen.github.io"
git add -A
git commit -m "update public html"
git push origin || (echo Push failed! & pause & exit /b 1)
popd

endlocal
