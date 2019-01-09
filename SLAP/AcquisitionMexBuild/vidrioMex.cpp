#include "vidrioMex.h"

void returnInt32Val(int32_t r, int nlhs, mxArray *plhs[])
{
	if(nlhs)
	{
		int *datPtr;
		mwSize sz = 1;

		plhs[0] = mxCreateNumericArray(1, &sz, mxINT32_CLASS, mxREAL);
		datPtr = (int*)mxGetData(plhs[0]);

		if(datPtr)
			*datPtr = (int)r;
	}
}

void returnLogicalVal(bool r, int nlhs, mxArray *plhs[])
{
	if(nlhs)
	{
		uint8_t *datPtr;
		mwSize sz = 1;

		plhs[0] = mxCreateNumericArray(1, &sz, mxLOGICAL_CLASS, mxREAL);
		datPtr = (uint8_t*)mxGetData(plhs[0]);

		if(datPtr)
			*datPtr = (int)r;
	}
}

bool getLogical(const mxArray *arr, bool defaultVal)
{
	void *dat = mxGetData(arr);

	if(!dat)
		return defaultVal;

	switch(mxGetClassID(arr))
	{
	case mxLOGICAL_CLASS :
	case mxINT8_CLASS :
	case mxUINT8_CLASS :
		return (*((uint8_t*)dat)) != 0;

	case mxINT16_CLASS :
	case mxUINT16_CLASS :
		return (*((uint16_t*)dat)) != 0;

	case mxINT32_CLASS :
	case mxUINT32_CLASS :
		return (*((uint32_t*)dat)) != 0;

	case mxINT64_CLASS :
	case mxUINT64_CLASS :
		return (*((uint64_t*)dat)) != 0;

	case mxSINGLE_CLASS :
		return (*((float*)dat)) != 0;

	case mxDOUBLE_CLASS :
		return (*((double*)dat)) != 0;

	default :
		return defaultVal;
	}
}


double getDouble(const mxArray *arr, double defaultVal)
{
	void *dat = mxGetData(arr);

	if(!dat)
		return defaultVal;

	switch(mxGetClassID(arr))
	{
	case mxLOGICAL_CLASS :
	case mxUINT8_CLASS :
		return (*((uint8_t*)dat));

	case mxINT8_CLASS :
		return (*((int8_t*)dat));

	case mxUINT16_CLASS :
		return (*((uint16_t*)dat));

	case mxINT16_CLASS :
		return (*((int16_t*)dat));

	case mxUINT32_CLASS :
		return (*((uint32_t*)dat));

	case mxINT32_CLASS :
		return (*((int32_t*)dat));

	case mxUINT64_CLASS :
		return (*((uint64_t*)dat));

	case mxINT64_CLASS :
		return (*((int64_t*)dat));

	case mxSINGLE_CLASS :
		return (*((float*)dat));

	case mxDOUBLE_CLASS :
		return (*((double*)dat));

	default :
		return defaultVal;
	}
}


int getInteger(const mxArray *arr, int defaultVal)
{
	void *dat = mxGetData(arr);

	if(!dat)
		return defaultVal;

	switch(mxGetClassID(arr))
	{
	case mxLOGICAL_CLASS :
	case mxUINT8_CLASS :
		return (*((uint8_t*)dat));

	case mxINT8_CLASS :
		return (*((int8_t*)dat));

	case mxUINT16_CLASS :
		return (*((uint16_t*)dat));

	case mxINT16_CLASS :
		return (*((int16_t*)dat));

	case mxUINT32_CLASS :
		return (*((uint32_t*)dat));

	case mxINT32_CLASS :
		return (*((int32_t*)dat));

	case mxUINT64_CLASS :
		return (*((uint64_t*)dat));

	case mxINT64_CLASS :
		return (*((int64_t*)dat));

	case mxSINGLE_CLASS :
		return (*((float*)dat));

	case mxDOUBLE_CLASS :
		return (*((double*)dat));

	default :
		return defaultVal;
	}
}

bool getMtlbObjPropertyDbl(const mxArray *mtlbObj, const char *propname, double * const val, double defaultVal)
{
	mxArray* propVal = mxGetProperty(mtlbObj,0,propname);
	if(propVal)
	{
		*val = (double)mxGetScalar(propVal);
		mxDestroyArray(propVal);
		return true;
	}
	else
	{
		*val = defaultVal;
		return false;
	}
}

bool getMtlbObjPropertyStr(const mxArray *mtlbObj, const char *propname, char **val)
{
	bool success = false;
	mxArray* propVal = mxGetProperty(mtlbObj,0,propname);
	*val = NULL;

	if(propVal)
	{
		size_t buflen = mxGetNumberOfElements(propVal) + 1;
		if(buflen > 1)
		{
			char *buf = new char[buflen];
			if(buf)
			{
				*val = buf;
				success = (mxGetString(propVal,buf,(mwSize)buflen) == 0);
			}
		}
		mxDestroyArray(propVal);
	}

	return success;
}