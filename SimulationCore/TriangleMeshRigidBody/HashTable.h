#ifndef __HASH_TABLE_H__
#define __HASH_TABLE_H__

#include <cstring>

class HashTable
{
protected:
	struct Item { size_t key, value; Item *next; };
	size_t table_size;
	Item **table;
	size_t item_num, item_id;
	Item *content;

public:
	HashTable(size_t t_size, size_t i_num = 0) : 
		table_size(t_size), table(nullptr),
		item_num(i_num), content(nullptr), item_id(0)
	{
		table = new Item*[table_size];
		memset(table, 0, table_size*sizeof(Item *));
		if (!item_num) item_num = table_size;
		content = new Item[item_num];
	}
	~HashTable()
	{
		if (table) delete[] table;
		if (content) delete[] content;
	}
	// do not check whether the key is already in table
	bool add_pair(size_t k, size_t v)
	{
		if (item_id < item_num)
		{
			Item *&entry = table[k%table_size];
			content[item_id].next = entry;
			entry = content + item_id;
			content[item_id].key = k;
			content[item_id].value = v;
			++item_id;
			return true;
		}
		return false;
	}
	bool get_pair(size_t k, size_t &v)
	{
		for (Item *cur = table[k%table_size]; cur; cur = cur->next)
		{
			if (cur->key == k)
			{
				v = cur->value;
				return true;
			}
		}
		return false;
	}
};

#endif