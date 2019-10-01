#ifndef __MEMORY_UTILITIES_ITEM_BUFFER_MULTI_ALLOC_HPP__
#define __MEMORY_UTILITIES_ITEM_BUFFER_MULTI_ALLOC_HPP__

#include <string>

// must be power of two
#define MEMORY_ALIGNMENT sizeof(void *)

// padding to round up the cloest power of two
#define MEMORY_ALIGNMENT_PADDING(address) \
	((MEMORY_ALIGNMENT - ((address) & (MEMORY_ALIGNMENT - 1))) & (MEMORY_ALIGNMENT - 1))

// round up the cloest power of two
#define MEMORY_ALIGNED_ADDRESS(address) \
	((address) + MEMORY_ALIGNMENT_PADDING(address))

namespace MemoryUtilities
{
	/*===============================================
	Class ItemBufferMultiAlloc
	-------------------------------------------------
	1. Can alloc multiple items;
	2. Can only be reset, does not support del().
	================================================*/
	template <class Item, size_t fold = 1, size_t pre_alloc_size = 0>
	class ItemBufferMultiAlloc
	{
	protected:
		struct MemPageHeader
		{
			MemPageHeader *next;
			Item *start;
			size_t size;
		};
		size_t base_page_size, page_size; // default size of new memory pool
		MemPageHeader *top_page; // memory pools managed as stack
		// first in-stack memory page
		union
		{
			MemPageHeader first_page;
			struct
			{
				MemPageHeader *first_page_next;
				Item *first_page_start;
				size_t first_page_size;
			};
			char first_page_mem[sizeof(MemPageHeader) + MEMORY_ALIGNMENT + sizeof(Item)*pre_alloc_size];
		};

		MemPageHeader *cur_page;
		Item *cur, *end;

	public:
		ItemBufferMultiAlloc(size_t init_page_size = pre_alloc_size) :
			top_page(&first_page),
			base_page_size(init_page_size ? init_page_size : 1), page_size(base_page_size),
			first_page_next(nullptr),
			first_page_start((Item *)MEMORY_ALIGNED_ADDRESS(size_t(first_page_mem) + sizeof(MemPageHeader))),
			first_page_size(pre_alloc_size),
			cur_page(&first_page), cur(first_page_start), end(first_page_start + first_page_size) {}
		~ItemBufferMultiAlloc() { clear(); }
		inline void set_page_size(size_t init_page_size) noexcept
		{
			base_page_size = init_page_size ? init_page_size : 1;
			page_size = base_page_size;
		}
		inline size_t get_page_size(void) const noexcept { return page_size; }
		inline Item *alloc(size_t num)
		{
			Item *res = cur;
			cur += num;
			if (cur > end) // not enough memory in this slot
			{
				while (next_pool)
				{
					if (size_t(next_pool->end) < num)
					{
						next_pool = next_pool->next;
						continue;
					}
					cur = next_pool->mem;
					end = cur + next_pool->size;
					next_pool = next_pool->next;
					goto memory_found;
				}
				new_mem_pool(num);
				next_pool = nullptr;
				cur = top_pool->mem;
				end = cur + top_pool->size;
			memory_found:	
				res = cur;
				cur += num;
			}
			return res;
		}

		// clear items in buffer
		// preserve allocated memory
		inline void reset(void)
		{
			// reset and delete all items
			if (top_pool)
			{
				next_pool = top_pool->next;
				cur = top_pool->mem;
				end = cur + top_pool->size;
			}
			else
			{
				next_pool = nullptr;
				cur = nullptr;
				end = nullptr;
			}
		}
		// Optimize performance by reallocating
		// the whole chunk of memory
		inline void reset_optimize(void)
		{
			// optimize performance
			if (pool_num > 1)
			{
				size_t merged_pool_size = total_size;
				clear();
				new_mem_pool(merged_pool_size);
			}
			// reset and delete all items
			if (top_pool)
			{
				next_pool = top_pool->next;
				cur = top_pool->mem;
				end = cur + top_pool->size;
			}
			else
			{
				next_pool = nullptr;
				cur = nullptr;
				end = nullptr;
			}
		}

		// clear allocated memory
		void clear(void)
		{
			char *tmp;
			while (top_pool)
			{
				tmp = reinterpret_cast<char *>(top_pool);
				top_pool = top_pool->next;
				delete[] tmp;
			}
			pool_num = 0;
			total_size = 0;
			next_pool = nullptr;
			cur = nullptr;
			end = nullptr;
		}

	protected:
		// allocate new memeory pool with size >= "size"
		void new_mem_pool(size_t size)
		{
			if (size == 0) return;
			if (size < pool_size)
				size = pool_size;
			size_t mem_size = sizeof(MemPageHeader) + size * sizeof(Item) + MEMORY_ALIGNMENT;
			char *mem = new char[mem_size];
			MemPageHeader *pool_header = reinterpret_cast<MemPageHeader *>(mem);
			pool_header->next = top_pool;
			top_pool = pool_header;
			pool_header->size = size;
			// memory alignment
			mem += sizeof(MemPageHeader);
			mem += MEMORY_ALIGNMENT_PADDING(size_t(mem));
			pool_header->mem = reinterpret_cast<Item *>(mem);
			++pool_num;
			total_size += size;
		}
	};
};

#undef MEMORY_ALIGNMENT
#undef MEMORY_ALIGNMENT_PADDING
#undef MEMORY_ALIGNED_ADDRESS

#endif