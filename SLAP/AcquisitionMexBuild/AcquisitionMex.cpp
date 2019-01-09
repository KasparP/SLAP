// AcquisitionMex.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "AcqEngine.h"
#include "Log.h"

#define MAX_COMMAND_LEN 10

cAcqEngine *gAcqEngine = NULL;

void onExit(void)
{
   if (gAsyncMex)
      AsyncMex_destroy(&gAsyncMex);
}

void asyncMexMATLABCallback(LPARAM lParam, void* params)
{
   string r = (gAcqEngine ? " (acq eng valid)" : " (acq eng invalid)");
   COUT("Entering async cb from matlab thread" + r);

   if (gAcqEngine)
      gAcqEngine->asyncCallback(gAcqEngine->m_finishCb);

   COUT("Exiting async cb");
}


void __cdecl mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
	char cmdStr[MAX_COMMAND_LEN];

   if (!gAsyncMex)
   {
      gAsyncMex = AsyncMex_create((AsyncMex_Callback*)asyncMexMATLABCallback, NULL);
      if (gAsyncMex)
         mexAtExit(&onExit);
   }

	if(nrhs < 1) {
		mexErrMsgTxt("MEX: No command specified.");
		returnInt32Val(-2, nlhs, plhs);
		return;
	}

	mxGetString(prhs[0],cmdStr,MAX_COMMAND_LEN);

   string r = (gAcqEngine ? " (acq eng valid)" : " (acq eng invalid)");
   COUT("Entering mex: " + string(cmdStr) + r);

	if(!strcmp(cmdStr, "init"))
	{
		if(!gAcqEngine)
			gAcqEngine = new cAcqEngine(prhs[1]);
		returnInt32Val((int32_t)gAcqEngine->m_numPages, nlhs, plhs);

      COUT("Exiting mex");
		return;
	}

	if(!strcmp(cmdStr, "start"))
	{
		if(gAcqEngine)
			returnInt32Val(gAcqEngine->start(prhs[1]), nlhs, plhs);
		else
			returnInt32Val(-100, nlhs, plhs);

      COUT("Exiting mex");
		return;
	}

	if(!strcmp(cmdStr, "finish"))
	{
		if(gAcqEngine)
			returnInt32Val(gAcqEngine->finish(prhs[1]), nlhs, plhs);
		else
			returnInt32Val(-100, nlhs, plhs);

      COUT("Exiting mex");
		return;
	}

	if(!strcmp(cmdStr, "stop"))
	{
		if(gAcqEngine)
			returnInt32Val(gAcqEngine->stop(), nlhs, plhs);
		else
			returnInt32Val(-100, nlhs, plhs);

      COUT("Exiting mex");
		return;
	}


	if(!strcmp(cmdStr, "exit"))
	{
		if(gAcqEngine)
		{
			delete gAcqEngine;
			gAcqEngine = NULL;
		}
		returnInt32Val(0, nlhs, plhs);

      COUT("Exiting mex");
		return;
	}

	mexErrMsgTxt("MEX: Invalid command.");
	returnInt32Val(-1, nlhs, plhs);
   COUT("Exiting mex (was err)");
}
