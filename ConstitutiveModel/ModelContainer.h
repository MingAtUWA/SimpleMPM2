#ifndef __Model_Container_H__
#define __Model_Container_H__

#include "ItemBuffer.hpp"

#include "LinearElasticity.h"
#include "ModifiedCamClay.h"

#define __Model_Container_Add_Model__(md_name)           \
protected:                                               \
MemoryUtilities::ItemBuffer<md_name> md_name ## _buffer; \
public:                                                  \
	md_name *add_## md_name ##(size_t num)               \
	{                                                    \
		md_name *res = md_name ## _buffer.alloc(num);    \
		for (size_t m_id = 0; m_id < num; ++m_id)        \
			new (res + m_id) (md_name);                  \
		return res;                                      \
	}


class ModelContainer
{
	__Model_Container_Add_Model__(LinearElasticity);
	__Model_Container_Add_Model__(ModifiedCamClay);
};

#undef __Model_Container_Add_Model__

#endif