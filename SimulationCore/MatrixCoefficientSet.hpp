#ifndef __Matrix_Coefficient_Set_hpp__
#define __Matrix_Coefficient_Set_hpp__

#include <iostream>
#include "ItemBuffer.hpp"

template <size_t dim_num, typename Scalar = double, typename Index = size_t>
class MatrixCoefficientSet
{
public:
	struct Coefficient
	{
		friend MatrixCoefficientSet;
	public:
		Index m_row;
		Index m_col;
		Scalar m_value;
		inline Scalar row(void) const { return m_row; }
		inline Scalar col(void) const { return m_col; }
		inline Scalar value(void) const { return m_value; }
	protected:
		Coefficient *prev_col, *next_col;
		Coefficient *prev_row, *next_row;
	};
protected:
	struct RowHeader
	{
		Coefficient *head;
		Coefficient *first_col;
		Coefficient *last_col;
		RowHeader() :
			head((Coefficient *)((char *)this + offsetof(RowHeader, first_col) - offsetof(Coefficient, next_col))),
			first_col(nullptr), last_col(head) {}
		inline void add_coef(Coefficient *coef)
		{
			coef->next_col = nullptr;
			coef->prev_col = last_col;
			last_col->next_col = coef;
			last_col = coef;
		}
		inline void del_coef(Coefficient *coef)
		{
			if (coef->next_col) // not the last item
			{
				coef->prev_col->next_col = coef->next_col;
				coef->next_col->prev_col = coef->prev_col;
			}
			else // the last item
			{
				last_col->prev_col->next_col = nullptr;
				last_col = coef->prev_col;
			}
		}
	};
	struct ColHeader
	{
		Coefficient *head;
		Coefficient *first_row;
		Coefficient *last_row;
		ColHeader() :
			head((Coefficient *)((char *)this + offsetof(ColHeader, first_row) - offsetof(Coefficient, next_row))),
			first_row(nullptr), last_row(head) {}
		inline void add_coef(Coefficient *coef)
		{
			coef->next_row = nullptr;
			coef->prev_row = last_row;
			last_row->next_row = coef;
			last_row = coef;
		}
		inline void del_coef(Coefficient *coef)
		{
			if (coef->next_row) // not the last item
			{
				coef->prev_row->next_row = coef->next_row;
				coef->next_row->prev_row = coef->prev_row;
			}
			else // the last item
			{
				last_row->prev_row->next_row = nullptr;
				last_row = coef->prev_row;
			}
		}
	};
	RowHeader row_header[dim_num];
	ColHeader col_header[dim_num];
	MemoryUtilities::ItemBuffer<Coefficient, 2> coef_buffer;

public:
	MatrixCoefficientSet()
	{
		coef_buffer.set_page_size(dim_num);
	}
	~MatrixCoefficientSet()
	{
		coef_buffer.clear();
	}
	void reset(void)
	{
		for (size_t i = 0; i < dim_num; ++i)
		{
			RowHeader &rh = row_header[i];
			rh.first_col = nullptr;
			rh.last_col = rh.head;
			ColHeader &ch = col_header[i];
			ch.first_row = nullptr;
			ch.last_row = ch.head;
		}
		coef_buffer.reset_optimzie();
	}
	inline void add_coefficient(Index row_id, Index col_id, Scalar value)
	{
		Coefficient *coef = coef_buffer.alloc();
		coef->m_row = row_id;
		coef->m_col = col_id;
		coef->m_value = value;
		row_header[row_id].add_coef(coef);
		col_header[col_id].add_coef(coef);
	}
	inline void del_coefficient(Coefficient *coef)
	{
		row_header[coef->m_row].del_coef(coef);
		col_header[coef->m_col].del_coef(coef);
		coef_buffer.del(coef);
	}
	void print(void)
	{
		for (size_t i = 0; i < dim_num; i++)
		{
			if (row_header[i].first_col)
			{
				for (Coefficient *iter = row_header[i].first_col;
					iter; iter = iter->next_col)
				{
					std::cout << "(" << iter->row() << ", " << iter->col() << ", " << iter->value() << ") ";
				}
				std::cout << "\n";
			}
		}
		std::cout << "\n";
	}
	void print_with_iter(void)
	{
		size_t row_id = begin()->row();
		for (MatrixCoefficientSet<5>::CoefficientIterator iter(begin()); iter != end(); ++iter)
		{
			if (row_id != iter->row())
			{
				row_id = iter->row();
				std::cout << "\n";
			}
			std::cout << "(" << iter->row() << ", " << iter->col() << ", " << iter->value() << ") ";
		}
		std::cout << "\n";
	}
public:
	// delete all coefficients in the idth row and col except the diagonal term
	// return the diagonal term
	inline Scalar del_col_and_row(Index id, double col_coefs[dim_num])
	{
		Coefficient *iter, *tmp;
		for (iter = row_header[id].first_col; iter;)
		{
			if (iter->m_col == id)
			{
				iter = iter->next_col;
			}
			else
			{
				tmp = iter;
				iter = iter->next_col;
				del_coefficient(tmp);
			}
		}
		memset(col_coefs, 0, sizeof(double) * dim_num);
		for (iter = col_header[id].first_row; iter;)
		{
			col_coefs[iter->m_row] += iter->m_value;
			if (iter->m_row == id)
			{
				iter = iter->next_row;
			}
			else
			{
				tmp = iter;
				iter = iter->next_row;
				del_coefficient(tmp);
			}
		}
		return col_coefs[id];
	}
	inline Scalar del_row(Index id)
	{
		Coefficient *iter, *tmp;
		double dig_term = 0.0;
		for (iter = row_header[id].first_col; iter;)
		{
			if (iter->m_col == id)
			{
				dig_term += iter->m_value;
				iter = iter->next_col;
			}
			else
			{
				tmp = iter;
				iter = iter->next_col;
				del_coefficient(tmp);
			}
		}
		return dig_term;
	}

public: // iteractor
	struct CoefficientIterator
	{
		friend MatrixCoefficientSet;
	protected:
		MatrixCoefficientSet *parent;
		size_t cur_row;
		Coefficient *cur_coef;
	public:
		CoefficientIterator() {}
		CoefficientIterator(const CoefficientIterator &another)
		{
			parent = another.parent;
			cur_row = another.cur_row;
			cur_coef = another.cur_coef;
		}
		inline bool operator != (const CoefficientIterator &another) const
		{
			return cur_row != another.cur_row || cur_coef != another.cur_coef;
		}
		inline bool operator == (const CoefficientIterator &another) const
		{
			return cur_row == another.cur_row && cur_coef == another.cur_coef;
		}
		inline void operator++(void)
		{
			if (cur_coef)
			{
				cur_coef = cur_coef->next_col;
				if (cur_coef)
					return;
			}
			for (++cur_row; cur_row < dim_num; ++cur_row)
			{
				cur_coef = parent->row_header[cur_row].first_col;
				if (cur_coef)
					return;
			}
		}
		Coefficient *operator->(void) {	return cur_coef; }
	};
	CoefficientIterator begin(void)
	{
		CoefficientIterator iter;
		iter.parent = this;
		for (size_t row_id = 0; row_id < dim_num; ++row_id)
		{
			RowHeader &rh = row_header[row_id];
			iter.cur_coef = rh.first_col;
			if (iter.cur_coef)
			{
				iter.cur_row = row_id;
				return iter;
			}
		}
		iter.cur_row = dim_num;
		return iter;
	}
	CoefficientIterator end(void)
	{
		CoefficientIterator iter;
		iter.parent = this;
		iter.cur_row = dim_num;
		iter.cur_coef = nullptr;
		return iter;
	}
};

#endif