#ifndef _MEMORY_UTILITIES_ITEM_BUFFER_HPP_
#define _MEMORY_UTILITIES_ITEM_BUFFER_HPP_

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
	Class ItemBuffer
	-------------------------------------------------
	1. This buffer just store item without sequence;
	2. Use linked memory pool;
	3. Set reasonable large page_size for efficiency;
	4. No check to ensure the item is in this buffer
	   when using del().
	================================================*/
	template<typename Item, size_t fold = 1, size_t pre_alloc_size = 0>
	class ItemBuffer
	{
	protected:
		struct MemPageHeader
		{
			MemPageHeader *next;
			Item *start, *end;
		};
		// first in-stack memory page
		union
		{
			char first_page_mem[sizeof(MemPageHeader) + MEMORY_ALIGNMENT + sizeof(Item)*pre_alloc_size];
			MemPageHeader first_page;
			struct
			{
				MemPageHeader *first_page_next;
				Item *first_page_start, *first_page_end;
			};
		};
		MemPageHeader *last_page;
		size_t base_page_size, page_size;
		bool need_optimized;
		// current state
		MemPageHeader *cur_page;
		Item *cur, *start, *end;
		
		struct EmptySlot { EmptySlot *next; };
		union Slot { Item item;	EmptySlot *next; };
		EmptySlot *empty; // head of allocated item

	public:
		ItemBuffer(size_t init_page_size = pre_alloc_size) :
			last_page(&first_page),
			base_page_size(init_page_size ? init_page_size : 1), page_size(base_page_size),
			need_optimized(false),
			first_page_next(nullptr),
			first_page_start((Item *)MEMORY_ALIGNED_ADDRESS(size_t(first_page_mem) + sizeof(MemPageHeader))),
			first_page_end(first_page_start + pre_alloc_size),
			cur_page(&first_page), cur(first_page_start),
			start(first_page_start), end(first_page_end),
			empty(nullptr) {}
		~ItemBuffer() { clear(); }
		inline void set_page_size(size_t init_page_size) noexcept
		{
			base_page_size = init_page_size ? init_page_size : 1;
			page_size = base_page_size;
		}
		inline size_t get_page_size(void) const noexcept { return page_size; }
		// allocate memory for new item
		inline Item *alloc(void)
		{
			Item *tmp;
			if (empty)
			{
				tmp = (Item *)empty;
				empty = empty->next;
			}
			else if (cur < end)
			{
				tmp = (Item *)cur;
				++cur;
			}
			else
			{
				tmp = alloc_from_next_page();
			}
			return tmp;
		}
		// no safety check if pitem is in the buffer.
		inline void del(Item *pitem)
		{
			EmptySlot *tmp = (EmptySlot *)(pitem);
			// add to empty list
			tmp->next = empty;
			empty = tmp;
		}
		inline void reset(void)
		{
			cur_page = &first_page;
			start = first_page_start;
			cur = start;
			end = first_page_end;
		}
		// contract allocated memory pages into one
		void reset_optimize(size_t additional_stack_size = 0)
		{
			if (need_optimized)
			{
				// cal total allocated size
				size_t total_size = 0;
				for (MemPageHeader *pg_iter = first_page.next; pg_iter; pg_iter = pg_iter->next)
					total_size += pg_iter->end - pg_iter->start;
				total_size += additional_stack_size;
				clear(); // clear old pages
				// alloc new page
				union { char *mem; MemPageHeader *mem_page; };
				total_size += pre_alloc_size;
				// take memory of the first page as "preallocated space"
				size_t char_num = sizeof(MemPageHeader) + MEMORY_ALIGNMENT + sizeof(Item)*(total_size+1);
				mem = new char[char_num];
				// abandon original in-stack space
				first_page_start = (Item *)MEMORY_ALIGNED_ADDRESS(size_t(mem) + sizeof(MemPageHeader));
				first_page_end = first_page_start + total_size;
				mem_page->start = first_page_end;
				mem_page->end = mem_page->start + 1;
				// add page to list
				mem_page->next = nullptr;
				last_page = mem_page;
				first_page.next = mem_page;
				// clear flag
				need_optimized = false;
			}
			// reset state
			cur_page = &first_page;
			start = first_page_start;
			end = first_page_end;
			cur = start;
		}
		void clear(void)
		{
			// reset page size
			page_size = base_page_size;
			// clear list
			MemPageHeader *&top_page = first_page.next;
			MemPageHeader *tmp_page = top_page;
			while (tmp_page)
			{
				top_page = top_page->next;
				delete[] (char*)tmp_page;
				tmp_page = top_page;
			}
			first_page_start = ((Item *)MEMORY_ALIGNED_ADDRESS(size_t(first_page_mem) + sizeof(MemPageHeader)));
			first_page_end = first_page_start + pre_alloc_size;
			last_page = &first_page;
			first_page.next = nullptr;
			empty = nullptr;
			need_optimized = false;
			// reset state
			cur_page = &first_page;
			start = first_page_start;
			cur = start;
			end = first_page_end;
		}
	protected:
		// allocate memory from other memory pool
		Item *alloc_from_next_page(void)
		{
			// alloc new page (if necessary)
			if (cur_page->next == nullptr)
			{
				// alloc new page
				union { char *mem; MemPageHeader *mem_page; };
				size_t char_num = sizeof(MemPageHeader) + MEMORY_ALIGNMENT + sizeof(Item)*page_size;
				mem = new char[char_num];
				// init new page
				mem_page->start = (Item *)MEMORY_ALIGNED_ADDRESS(size_t(mem) + sizeof(MemPageHeader));
				mem_page->end = mem_page->start + page_size;
				mem_page->next = nullptr;
				// add to list
				last_page->next = mem_page;
				last_page = mem_page;
				need_optimized = true;
				// update page_size
				page_size *= fold;
			}
			// move the next page
			cur_page = cur_page->next;
			start = cur_page->start;
			cur = start;
			end = cur_page->end;
			return cur++;
		}
	};
}

#undef MEMORY_ALIGNMENT
#undef MEMORY_ALIGNMENT_PADDING
#undef MEMORY_ALIGNED_ADDRESS

#endif