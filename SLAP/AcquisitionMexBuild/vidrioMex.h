#pragma once
#include "stdafx.h"
#include "mex.h"

void returnInt32Val(int32_t r, int nlhs, mxArray *plhs[]);

void returnLogicalVal(bool r, int nlhs, mxArray *plhs[]);

bool getLogical(const mxArray *arr, bool defaultVal);

double getDouble(const mxArray *arr, double defaultVal);

int getInteger(const mxArray *arr, int defaultVal);

bool getMtlbObjPropertyDbl(const mxArray *mtlbObj, const char *propname, double * const val, double defaultVal);

bool getMtlbObjPropertyStr(const mxArray *mtlbObj, const char *propname, char **val);
