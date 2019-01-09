#include "stdafx.h"
#include "AcqEngine.h"
#include "Log.h"

AsyncMex *gAsyncMex;


cAcqEngine::cAcqEngine(const mxArray *sObj)
{
	double val;
	int ctr;
	
	m_fpgaStat = NiFpga_Initialize();
	m_pmtReadingThdRunning = false;
	m_pmtWritingThdRunning = false;
	m_galvoLogging = false;
	m_acqStopReq = false;
   m_acqStopWhenFinishedReq = false;
	m_metaFilename = NULL;
	m_pDataFilename = NULL;
	m_gDataFilename = NULL;
	m_headerString = NULL;
   m_finishCb = NULL;
	
	//Get memory config parameters
	getMtlbObjPropertyDbl(sObj, "memoryPages", &val, 1000);
	m_numPages = (uint32_t)val;
	
	getMtlbObjPropertyDbl(sObj, "memoryPageSize", &val, 250000);
	m_pageSize = (uint32_t)val;

	//try to allocate the memory
	m_ppDataBufs = new uint64_t*[m_numPages];
	m_pPageElements = new size_t[m_numPages];
	m_pPageDone = new bool[m_numPages];
	for(ctr = 0; ctr < m_numPages; ctr++)
	{
		try
		{
			m_ppDataBufs[ctr] = NULL;
			m_ppDataBufs[ctr] = new uint64_t[m_pageSize];
			memset(m_ppDataBufs[ctr],0,m_pageSize*8);			//doing this memset ensures the memory really is allocated
		}
		catch(exception e)
		{
			if(m_ppDataBufs[ctr])
				delete[] m_ppDataBufs[ctr];
			m_ppDataBufs[ctr] = NULL;
			break;
		}
	}
	m_numPages = ctr;
}


cAcqEngine::~cAcqEngine(void)
{
	stop();
	NiFpga_Finalize();
	
	if(m_metaFilename)
		delete[] m_metaFilename;
	if(m_pDataFilename)
		delete[] m_pDataFilename;
	if(m_gDataFilename)
		delete[] m_gDataFilename;
	if(m_headerString)
		delete[] m_headerString;

	for(m_pageSize = 0; m_pageSize < m_numPages; m_pageSize++)
	{
		try
		{
			if(m_ppDataBufs[m_pageSize])
			{
				delete[] m_ppDataBufs[m_pageSize];
				m_ppDataBufs[m_pageSize] = NULL;
			}
		}
		catch(exception e)
		{
			
		}
	}
	if(m_pPageElements)
		delete[] m_pPageElements;
	m_pPageElements = NULL;
	if(m_pPageDone)
		delete[] m_pPageDone;
	m_pPageDone = NULL;
}


int32_t cAcqEngine::start(const mxArray *sObj)
{
	double val;

   cancelCallback();

	if(m_pmtReadingThdRunning || m_pmtWritingThdRunning || m_galvoLogging)
	{
		int32_t err = stop();

		if(err)
			return err;
	}

	//Get fpga session
	if(getMtlbObjPropertyDbl(sObj, "fpgaSessionId", &val, 0))
		m_fpgaSession = (NiFpga_Session)val;
	else
		return -1;

	//Get fifo numbers
	if(getMtlbObjPropertyDbl(sObj, "pmtFifoNum", &val, 0))
		m_pmtFifoNum = (uint32_t)val;
	else
		return -2;
	if(getMtlbObjPropertyDbl(sObj, "galvoFifoNum", &val, 0))
		m_galvoFifoNum = (uint32_t)val;
	else
		return -3;

	//Get the filenames
	if(m_metaFilename)
		delete[] m_metaFilename;
	if(m_pDataFilename)
		delete[] m_pDataFilename;
	if(m_gDataFilename)
		delete[] m_gDataFilename;
	getMtlbObjPropertyStr(sObj, "metaFileName", &m_metaFilename);
	getMtlbObjPropertyStr(sObj, "pmtDataFileName", &m_pDataFilename);
	getMtlbObjPropertyStr(sObj, "galvoDataFileName", &m_gDataFilename);
	getMtlbObjPropertyStr(sObj, "metaHeaderString", &m_headerString);
	
	for(m_writePointer = 0; m_writePointer < m_numPages; m_writePointer++)
	{
		m_pPageElements[m_writePointer] = 0;
		m_pPageDone[m_writePointer] = false;
	}
	m_writePointer = 0;
	m_readPointer = 0;
	m_pmtFileWriteErrors = 0;
	m_galvoFileWriteErrors = 0;

	m_hPmtFifoReadThd = _beginthreadex(NULL, 0, &cAcqEngine::pmtFifoReadThreadFcn, (LPVOID)this, 0, NULL);
	if(m_hPmtFifoReadThd == 0)
		return -11;

	m_hPmtFileWriteThd = _beginthreadex(NULL, 0, &cAcqEngine::pmtFileWriteThreadFcn, (LPVOID)this, 0, NULL);
	if(m_hPmtFileWriteThd == 0)
		return -12;

	m_hGalvoLoggingThd = _beginthreadex(NULL, 0, &cAcqEngine::galvoLoggingThreadFcn, (LPVOID)this, 0, NULL);
	if(m_hGalvoLoggingThd == 0)
		return -13;
	
	return 0;
}

int32_t cAcqEngine::finish(const mxArray *finishCb)
{
   if (m_pmtReadingThdRunning || m_pmtWritingThdRunning || m_galvoLogging)
   {
      if (m_finishCb)
         mxDestroyArray(m_finishCb);
      m_finishCb = NULL;

      // prepare callback
      if(finishCb)
         m_finishCb = mxDuplicateArray(finishCb);

      if (m_finishCb)
      {
         mexMakeArrayPersistent(m_finishCb);
         m_acqStopWhenFinishedReq = true;
      }
   }
   else
   {
      // call callback now
      asyncCallback(finishCb);
   }

   return 0;
}

void cAcqEngine::cancelCallback(void)
{
   if (m_finishCb)
   { // this order of operations will not create a race condition!!
      mxArray *temp = m_finishCb;
      m_finishCb = NULL;
      mxDestroyArray(temp);
   }
   m_acqStopWhenFinishedReq = false;
}

void cAcqEngine::asyncCallback(const mxArray *finishCb)
{
   if (finishCb)
   {
      // call callback now
      mxArray **rhs;
      mxArray* mException = NULL;

      rhs = (mxArray **)&finishCb;
      mException = mexCallMATLABWithTrap(0, NULL, 1, rhs, "feval");

      if (mException)
      {
         char *message = NULL;

         OutputDebugStringA("Matlab Callback Error\n");

         getMtlbObjPropertyStr(mException, "message", &message);

         if (message)
         {
            mexPrintf("Error in AcquisitionMex finish callback. Error message:\n%s\n", message);
            delete[] message;
         }
         else
            mexPrintf("Error in AcquisitionMex finish callback.\n");

         mxDestroyArray(mException);
      }
   }
}


int32_t cAcqEngine::stop(void)
{
   cancelCallback();

	if(m_pmtReadingThdRunning || m_pmtWritingThdRunning || m_galvoLogging)
	{
		DWORD st;

		m_acqStopReq = true;
		st = GetTickCount();

		while(m_pmtReadingThdRunning || m_pmtWritingThdRunning || m_galvoLogging)
		{
			if((GetTickCount() - st) > 5000)
			{
				m_acqStopReq = false;
				return 1;
			}
			Sleep(10);
		}
		
		m_acqStopReq = false;
	}
	return 0;
}


unsigned int WINAPI cAcqEngine::pmtFifoReadThreadFcn(LPVOID userData)
{
	cAcqEngine *pA = (cAcqEngine*)userData;
	NiFpga_Status status;
	size_t elementsRemaining;
	uint64_t *wp = pA->m_ppDataBufs[pA->m_writePointer];
	size_t currentPageFreeElements = pA->m_pageSize;

	pA->m_pmtReadingThdRunning = true;
   COUT("PMT FIFO read thread started");

	while(!pA->m_acqStopReq)
	{
		//make sure we are ready to read more data
		if(pA->m_pPageDone[pA->m_writePointer])
		{
			//the current page is full, meaning we have caught up to the file write thread in the revovling buffer. No choice but to wait.
			Sleep(1);
		}
		else
		{
			status = NiFpga_ReadFifoU64(pA->m_fpgaSession, pA->m_pmtFifoNum, wp, currentPageFreeElements, 0, &elementsRemaining);
			
			if(status == NiFpga_Status_Success)
			{
            // the next block of code use elementsRemaining as an indication of how many elements were read
				elementsRemaining = currentPageFreeElements;
			}
			else if(status == NiFpga_Status_FifoTimeout)
			{
				if(elementsRemaining > 0)
				{
					// the read timed out but there are some elements available. lets read them
					status = NiFpga_ReadFifoU64(pA->m_fpgaSession, pA->m_pmtFifoNum, wp, elementsRemaining, 0, NULL);
				}
				else if (pA->m_pPageElements[pA->m_writePointer])
				{
					//nothing to read. if this page was partially written move on to the next page so that the reading thread will consume it
					wp = pA->advanceWritePointer();
					currentPageFreeElements = pA->m_pageSize;
				}
				else if (pA->m_acqStopWhenFinishedReq)
					break; //if we have a soft quit request, quit now
				else
					Sleep(2); // there was no data in the fifo and nothing to do. give the cpu a little break
			}

			if(status == NiFpga_Status_Success)
			{
				pA->m_pPageElements[pA->m_writePointer] += elementsRemaining;
				wp += elementsRemaining;
				currentPageFreeElements -= elementsRemaining;

				if(!currentPageFreeElements)
				{
					// the current page is full. advance the write pointer to the next one
					wp = pA->advanceWritePointer();
					currentPageFreeElements = pA->m_pageSize;
				}
			}
		}
	}

   COUT("PMT FIFO read thread exiting");
	pA->m_pmtReadingThdRunning = false;
	return 0;
}

uint64_t *cAcqEngine::advanceWritePointer(void)
{
	m_pPageDone[m_writePointer] = true;
	m_writePointer++;
	if(m_writePointer == m_numPages)
		m_writePointer = 0;
	return m_ppDataBufs[m_writePointer];
}


unsigned int WINAPI cAcqEngine::pmtFileWriteThreadFcn(LPVOID userData)
{
	cAcqEngine *pA = (cAcqEngine*)userData;
	uint64_t samplesWritten = 0;
	fstream metaFile;
	fstream dataFile;
   DWORD st;
   DWORD timeout;
	size_t elements;

	pA->m_pmtWritingThdRunning = true;
   COUT("PMT logger thread started");

	// write meta data file
	metaFile.open(pA->m_metaFilename, fstream::out | fstream::app);
	if(metaFile.is_open())
	{
		metaFile << "SLAPMi data log" << endl;
		metaFile << pA->m_headerString;
		metaFile.close();
	}

	// open data file
	dataFile.open(pA->m_pDataFilename, fstream::out | fstream::app | fstream::binary);

	// log data
	if(dataFile.is_open())
	{
		while(!pA->m_acqStopReq)
		{
			if(pA->m_pPageDone[pA->m_readPointer])
			{
            //the page at the current position is ready to be written to disk. write it and advance the read pointer
				elements = pA->m_pPageElements[pA->m_readPointer];
				if(elements)
				{
					samplesWritten += elements;

					try
					{
						dataFile.write((char*)pA->m_ppDataBufs[pA->m_readPointer], elements*8);
					}
					catch(exception e)
					{
						pA->m_pmtFileWriteErrors++;
					}

					pA->m_pPageElements[pA->m_readPointer] = 0;
					pA->m_pPageDone[pA->m_readPointer] = false;
				}

				pA->m_readPointer++;
				if(pA->m_readPointer == pA->m_numPages)
					pA->m_readPointer = 0;
			}
         else
         {
            //the page at the current position is not ready to be written to disk. do nothing

            if (pA->m_acqStopWhenFinishedReq && !pA->m_pmtReadingThdRunning)
               break;
            else
               Sleep(2);
         }
		}

		dataFile.close();
	}


   // if this is a soft stop request, give the galvo logging thread lots of time to finish
   if (!pA->m_acqStopReq)
      timeout = 20000;
   else
      timeout = 5000;

	// wait for galvo data thread to quit
	st = GetTickCount();
	while(pA->m_galvoLogging)
	{
		if((GetTickCount() - st) > timeout)
		{
         COUT("Galvo logger didnt quit?");
			break;
		}
		Sleep(10);
	}
	
	if(samplesWritten || pA->m_gSamplesWritten)
	{
		// write samples written to meta file
		metaFile.open(pA->m_metaFilename, fstream::out | fstream::app);
		if(metaFile.is_open())
		{
			metaFile << "pmtSamplesWritten=" << samplesWritten*2 << endl;
			metaFile << "galvoSamplesWritten=" << pA->m_gSamplesWritten << endl;
			metaFile.close();
		}
	}
	else
	{
		//delete meta and data files
		remove(pA->m_metaFilename);
		remove(pA->m_gDataFilename);
		remove(pA->m_pDataFilename);
	}

   if(pA->m_acqStopWhenFinishedReq)
      AsyncMex_postEventMessage(gAsyncMex, NULL);
	
   COUT("PMT logger thread exiting");
	pA->m_pmtWritingThdRunning = false;
	return 0;
}


unsigned int WINAPI cAcqEngine::galvoLoggingThreadFcn(LPVOID userData)
{
	cAcqEngine *pA = (cAcqEngine*)userData;
	fstream dataFile;
	uint64_t *pData;
	NiFpga_Status status;
	size_t elementsRemaining;

	pA->m_galvoLogging = true;
	pA->m_gSamplesWritten = 0;
   COUT("Galvo logger thread started");

	//open data file
	dataFile.open(pA->m_gDataFilename, fstream::out | fstream::app | fstream::binary);

	//log data
	if(dataFile.is_open())
	{
		pData = new uint64_t[25000];

		while(!pA->m_acqStopReq)
		{
			status = NiFpga_ReadFifoU64(pA->m_fpgaSession, pA->m_galvoFifoNum, pData, 25000, 0, &elementsRemaining);
			
			if(status == NiFpga_Status_Success)
			{
				elementsRemaining = 25000;
			}
			else if(status == NiFpga_Status_FifoTimeout)
			{
				if (elementsRemaining > 0)
					status = NiFpga_ReadFifoU64(pA->m_fpgaSession, pA->m_galvoFifoNum, pData, elementsRemaining, 0, NULL);
				else if (pA->m_acqStopWhenFinishedReq)
					break; // handle soft stop request
				else
					Sleep(2); // nothing to read. give the cpu a break
			}

			if(status == NiFpga_Status_Success)
			{
				try
				{
					dataFile.write((char*)pData, elementsRemaining*8);
				}
				catch(exception e)
				{
					pA->m_galvoFileWriteErrors++;
				}
			}
		}

		delete[] pData;
		dataFile.close();
	}


   COUT("Galvo logger thread exiting");
	pA->m_galvoLogging = false;
	return 0;
}




