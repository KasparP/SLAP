#pragma once
#include "vidrioMex.h"
#include "AsyncMex.h"

extern AsyncMex *gAsyncMex;

class cAcqEngine
{
public:
	NiFpga_Status m_fpgaStat;
	bool m_acqStopReq;
   bool m_acqStopWhenFinishedReq;

	bool m_pmtReadingThdRunning;
	bool m_pmtWritingThdRunning;
	bool m_galvoLogging;
	uintptr_t m_hPmtFifoReadThd;
	uintptr_t m_hPmtFileWriteThd;
	uintptr_t m_hGalvoLoggingThd;

	uint64_t m_gSamplesWritten;

	char *m_metaFilename;
	char *m_pDataFilename;
	char *m_gDataFilename;
	char *m_headerString;
	
	size_t m_numPages;
	size_t m_numPagesAllocated;
	size_t m_pageSize;
	size_t m_writePointer;
	size_t m_readPointer;
	uint64_t **m_ppDataBufs;
	size_t *m_pPageElements;
	bool *m_pPageDone;
	
	size_t m_pmtFileWriteErrors;
	size_t m_galvoFileWriteErrors;

	NiFpga_Session m_fpgaSession;
	uint32_t m_pmtFifoNum;
	uint32_t m_galvoFifoNum;

   mxArray *m_finishCb;
   AsyncMex *m_AsyncMex;

	cAcqEngine(const mxArray *sObj);
	~cAcqEngine();

	int32_t start(const mxArray *sObj);
   int32_t finish(const mxArray *finishCb);
   void cancelCallback(void);
   void asyncCallback(const mxArray *finishCb);
	int32_t stop();

	uint64_t *advanceWritePointer(void);
	
	static unsigned int WINAPI pmtFifoReadThreadFcn(LPVOID userData);
	static unsigned int WINAPI pmtFileWriteThreadFcn(LPVOID userData);
	static unsigned int WINAPI galvoLoggingThreadFcn(LPVOID userData);
};

