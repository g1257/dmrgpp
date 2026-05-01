#pragma once

#include <H5Cpp.h>

/* Class to disable HDF5 exception printing.  It is re-enabled when the
 * class is destroyed.  Can comment out dontPrint() call for debugging.
 */
class HDF5DisableExceptionPrinting {
public:

	HDF5DisableExceptionPrinting()
	{
		H5::Exception::getAutoPrint(_func, &_clientData);
		H5::Exception::dontPrint();
	}
	~HDF5DisableExceptionPrinting() { H5::Exception::setAutoPrint(_func, _clientData); }

private:

	H5E_auto2_t _func;
	void*       _clientData;
};
