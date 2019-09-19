#ifndef __CONTACT_STATE_H__
#define __CONTACT_STATE_H__

#include <cstddef>

struct ContactState
{
public:
	double ftx_con, fty_con; // tangential force
	bool prev_in_contact; // previously in contact = in ContactStateList
	bool cur_in_contact; // currently in contact
	ContactState *_prev, *_next;
public:
	inline void init(void) noexcept
	{
		ftx_con = 0.0;
		fty_con = 0.0;
		prev_in_contact = false;
		cur_in_contact = false;
		_prev = nullptr;
		_next = nullptr;
	}
	inline ContactState *prev(void) noexcept { return _prev; }
	inline ContactState *next(void) noexcept { return _next; }
	template <typename ParticleType>
	inline ParticleType &get_pcl(void)
	{
		return *(ParticleType *)((char *)this - offsetof(ParticleType, contact_state));
	}
};

// double list of contact state
struct ContactStateList
{
	ContactState head;
	ContactStateList() { head._prev = &head; head._next = &head; }
	inline ContactState *first(void) noexcept { return head._next; }
	inline ContactState *last(void) noexcept { return head._prev; }
	inline ContactState *end(void) noexcept { return &head; }
	inline void add(ContactState &cs) noexcept
	{
		cs.ftx_con = 0.0;
		cs.fty_con = 0.0;
		cs.prev_in_contact = true;
		cs.cur_in_contact = true;
		cs._prev = head._prev;
		cs._next = &head;
		head._prev->_next = &cs;
		head._prev = &cs;
	}
	inline void del(ContactState &cs) noexcept
	{
		cs.ftx_con = 0.0;
		cs.fty_con = 0.0;
		cs.prev_in_contact = false;
		cs.cur_in_contact = false;
		cs._prev->_next = cs._next;
		cs._next->_prev = cs._prev;
		cs._prev = nullptr;
		cs._next = nullptr;
	}
	inline void reset_all_contact_state(void)
	{
		ContactState *pcs = head._next;
		ContactState *pcs_end = &head;
		while (pcs != pcs_end)
		{
			pcs->cur_in_contact = false;
			pcs = pcs->_next;
		}
	}
	inline void clear_not_in_contact(void)
	{
		ContactState *pcs = head._next;
		ContactState *pcs_end = &head;
		ContactState *pcs_tmp;
		while (pcs != pcs_end)
		{
			if (pcs->cur_in_contact)
				pcs = pcs->_next;
			else
			{
				pcs_tmp = pcs;
				pcs = pcs->_next;
				del(*pcs_tmp);
			}
		}
	}
	// Helper functions
	// Assumption: ParticleType has member "ContactState contact_state"
	template <typename ParticleType>
	inline void add_pcl(ParticleType &pcl) noexcept { add(pcl.contact_state); }
	template <typename ParticleType>
	inline void del_pcl(ParticleType &pcl) noexcept { del(pcl.contact_state); }
};

#endif