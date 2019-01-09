/**
 * @file AsyncMex.h
 *
 * @brief Interface to a generic module for passing asynchronous events back onto the Matlab thread.
 *
 * @author Timothy O'Connor
 * @date 4/24/10
 *
 * <em><b>Copyright</b> - Northwestern University/Howard Hughes Medical Institute 2010</em>
 *
 */

//#pragma once

#pragma message("Start of AsyncMex\n")

#ifndef __AsyncMex_h_ //Multiple include protection.
#pragma message("First time finding AsyncMex\n")
#define __AsyncMex_h_

#ifdef __cplusplus
extern "C" {
#endif

/*********************************
*         OS DEFINES            *
*********************************/
/**
 * @brief This is required to use Windows functions, which are the core of the event passing.
 *
 * This declares that the minumum system used must be NT4.0.
 * @see http://msdn.microsoft.com/en-us/library/ms686352.aspx
 */
//#define _WIN32_WINNT 0x0400

#include <limits.h>
#include <windows.h>

/**
 * @brief A version number for this object.
 *
 * This number should get incremented when changes are made to the code.
 */
#define ASYNCMEX_VERSION 0.1
#ifdef _AsyncMex_c
///@brief The string used to generate the unique identifier (determined at runtime) for AsyncMex events.
LPCTSTR ASYNCMEX_WINDOWMESSAGE_NAME = "AsyncSLAPMiMex_Event";
///@brief The unique identifier (determined at runtime) for AsyncMex events.
UINT ASYNCMEX_WINDOWMESSAGE_ID;
///@brief The ID of the (one and only) hook function.
HHOOK ASYNCMEX_MESSAGE_PUMP_HOOK_ID;
///@brief The string used to generate the unique identifier (determined at runtime) for AsyncMex events.
LPCTSTR ASYNCMEX_WINDOWCLASS_NAME = "AsyncSLAPMiMex_Window_Class";
#endif
///@brief The prototype for the callback function, which must be implemented by the user code.
typedef void (AsyncMex_Callback)(LPARAM lParam, void* userData);

/**
 * @brief The representation of an <tt>AsyncMex</tt> object.
 */
typedef struct
{
  ///@brief The handle to the Matlab thread.
  HANDLE matlabThread;
  ///@brief The ID of the Matlab thread.
  DWORD matlabThreadID;
  ///@brief The message ID used to filter messages.
  int messageID;
  ///@brief The ID of the hook function, necessary for later removal.
  HHOOK messagePumpHookID;
  ///@brief Generic user data that may be passed to the callback.
  void* userData;
  ///@brief The pointer to the user implemented callback function.
  AsyncMex_Callback* callback;
  ///@brief Handle to the window used to recieve messages.
  HWND hwnd;
  ///@brief The class used when creating the message processing client window.
  WNDCLASSEX wndClass;

} AsyncMex;

/**
 * @brief Prints the Windows error message corresponding to <tt>GetLastError()</tt>.
 * @see http://msdn.microsoft.com/en-us/library/ms679351(VS.85).aspx
 */
void AsyncMex_printWindowsErrorMessage(int lineNum);

/**
 * @brief Renders a <tt>AsyncMex</tt> object into a string.
 * @arg <tt>asyncM</tt> - The <tt>AsyncMex</tt> to be rendered.
 * @arg <tt>buff</tt> - The buffer into which to render the string.
 * @arg <tt>buffSize</tt> - The size of the buffer, in bytes.
 * @arg <tt>prefix</tt> - An arbitrary prefix, to control formating, may not be NULL but may be an empty string (ie "").
 * @return Number of characters written.
 */
int AsyncMex_structToString(AsyncMex* asyncM, char* buff, size_t buffSize, const char* prefix);

/**
 * @brief Creates a new <tt>AsyncMex</tt> object.
 * @arg <tt>callback</tt> - The callback used to handle events on the Matlab thread.
 * @arg <tt>userData</tt> - Arbitrary user data that will get passed to the callback at runtime.
 * @return A newly allocated <tt>AsyncMex</tt> object, NULL if unsuccessful.
 */
AsyncMex* AsyncMex_create(AsyncMex_Callback* callback, void* userData);

/**
 * @brief Destructor for an <tt>AsyncMex</tt> object.
 * The pointer passed in is set to NULL, to help prevent illegal memory access later.
 * @arg <tt>asyncM</tt> - The <tt>AsyncMex</tt> to be destroyed.
 */
void AsyncMex_destroy(AsyncMex** asyncM);

/**
 * @brief Posts an event to the Matlab thread, via an <tt>AsyncMex</tt> object.
 * @arg <tt>asyncM</tt> - The <tt>AsyncMex</tt> to be used to pass an event.
 * @arg <tt>lParam</tt> - Event-specific information.
 * @see PostThreadMessage ( http://msdn.microsoft.com/en-us/library/ms644946(VS.85).aspx )
 */
int AsyncMex_postEventMessage(AsyncMex* asyncM, LPARAM lParam);

/**
 * @brief Creates a window, for programs that need to process messages.
 * @arg <tt>asyncM</tt> - The <tt>AsyncMex</tt> for which to create a window.
 * @arg <tt>windowProc</tt> - The function to handle the message processing.
 * @return 0 if successful, non-zero otherwise.
 */
//int AsyncMex_createClientWindow(AsyncMex* asyncM, AsyncMex_WinProcCallback windowProc);

/**
 * @brief Tear down the client window used for Windows messaging.
 * @arg <tt>asyncM</tt> - he <tt>AsyncMex</tt> for which to destroy a window.
 * @return 0 if successful, non-zero otherwise.
 */
//int MCT_destroyClientWindow(AsyncMex* asyncM);
#ifdef __cplusplus
}
#endif

#endif