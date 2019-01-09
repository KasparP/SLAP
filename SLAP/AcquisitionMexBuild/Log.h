
#include <string>
#include <iostream>
#include <fstream>

//#define CONSOLEDBG
#define LOGFILE

#ifdef CONSOLEDBG
   #define STARTCONSOLE AllocConsole(); freopen("CONOUT$", "w", stdout);
   #define STOPCONSOLE FreeConsole();
   #define COUT(x) cout << x
#elif defined LOGFILE
   #define STARTCONSOLE
   #define STOPCONSOLE
   #define COUT(x) writeToLog(x)
#else
   #define STARTCONSOLE
   #define STOPCONSOLE
   #define COUT(x)
#endif

using namespace std;

void writeToLog(const string str);
