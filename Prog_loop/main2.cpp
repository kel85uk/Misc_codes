
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <stdio.h>
#include <windows.h>

using namespace std;

int main()
{
    PROCESS_INFORMATION pi;
    STARTUPINFO si;
    
    ZeroMemory( &si, sizeof(si) );
    si.cb = sizeof(si);
    ZeroMemory( &pi, sizeof(pi) );
    
    // Call up the 1st .exe file, initial_input.exe, which records 2 values input by the user
    if(!CreateProcess("c:\\Documents and Settings\\R.Dunphy\\My Documents\\NetBeansProjects\\Batch file\\initial_input.exe", NULL, NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi))
    {
        printf( "CreateProcess failed (%d). \n", GetLastError() );
    }
    else
    {
        printf("Process Creation Success");
    }
    // Wait until initial_input.exe exits.
    WaitForSingleObject( pi.hProcess, INFINITE);
    
    // Call up the 2nd .exe file, sum.exe, which adds the 2 values recorded by initial_input.exe
    if(!CreateProcess("c:\\Documents and Settings\\R.Dunphy\\My Documents\\NetBeansProjects\\Batch file\\sum.exe", NULL, NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi))
    {
        printf( "CreateProcess failed (%d). \n", GetLastError() );
    }
    else
    {
        printf("Process Creation Success");
    }
    // Wait until sum.exe exits.
    WaitForSingleObject( pi.hProcess, INFINITE);
    
    // Call up the 3rd .exe file, product.exe, which reads the sum from sum.exe, and multiplies it by 2. 
    if(!CreateProcess("c:\\Documents and Settings\\R.Dunphy\\My Documents\\NetBeansProjects\\Batch file\\product.exe", NULL, NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi))
    {
        printf( "CreateProcess failed (%d). \n", GetLastError() );
    }
    else
    {
        printf("Process Creation Success");
        return 0;
    }
    
    CloseHandle( pi.hProcess );
    CloseHandle( pi.hThread ); 
}

