#ifndef __HASH_TABLE_EDGE_H__
#define __HASH_TABLE_EDGE_H__

#include "ItemBuffer.hpp"
#include "TriangleMesh.h"

class HashTableEdge
{
protected:
	struct EdgeInfo
	{
		size_t n1, n2;
		size_t elem_num;
		EdgeInfo *next;
	};
	size_t table_size;
	EdgeInfo **table;
	MemoryUtilities::ItemBuffer<EdgeInfo> mem;

public:
	HashTableEdge(size_t node_num) :
		table_size(node_num/2+1), table(nullptr)
	{
		table = new EdgeInfo*[table_size];
		memset(table, 0, table_size * sizeof(EdgeInfo*));
		mem.set_page_size(table_size * 3);
	}

	~HashTableEdge()
	{
		if (table) delete[] table;
		table = nullptr;
		mem.clear();
	}

	bool add_edge(size_t n1, size_t n2)
	{
		if (n1 > n2)
		{
			size_t n_tmp = n1;
			n1 = n2;
			n2 = n_tmp;
		}

		EdgeInfo *tmp;
		if (tmp = find_edge(n1, n2))
		{
			++tmp->elem_num;
		}
		else
		{
			EdgeInfo *&top = table[n1%table_size];
			tmp = mem.alloc();
			tmp->next = top;
			top = tmp;
			tmp->n1 = n1;
			tmp->n2 = n2;
			tmp->elem_num = 1;
			return true;
		}

		return false;
	}
	
	size_t get_edge_at_boundary(TriangleMesh::Edge *&edges)
	{
		size_t edge_num = 0;
		for (size_t i = 0; i < table_size; ++i)
		{
			for (EdgeInfo *pe = table[i]; pe; pe = pe->next)
			{
				if (pe->elem_num == 1)
					++edge_num;
			}
		}

		if (edge_num == 0)
		{
			edges = nullptr;
			return 0;
		}

		edges = new TriangleMesh::Edge[edge_num];
		size_t k = 0;
		for (size_t i = 0; i < table_size; ++i)
		{
			for (EdgeInfo *pe = table[i]; pe; pe = pe->next)
			{
				if (pe->elem_num == 1)
				{
					edges[k].n1 = pe->n1;
					edges[k].n2 = pe->n2;
					++k;
				}	
			}
		}

		return edge_num;
	}

protected:
	EdgeInfo *find_edge(size_t n1, size_t n2)
	{
		if (n1 == 23)
			int egrgr = 0;
		for (EdgeInfo *cur = table[n1%table_size]; cur; cur = cur->next)
		{
			if (cur->n1 == n1 && cur->n2 == n2)
				return cur;
		}
		return nullptr;
	}
};

#endif