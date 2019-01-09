
#include "Log.h"

void writeToLog(const string str)
{
   ofstream f("C:\\acqMexLog.txt", ofstream::app);
   f << str << endl;
   f.close();
}