#ifndef __TIME_HISTORY_OUTPUT_CONSOLE_PROGRESS_BAR_H__
#define __TIME_HISTORY_OUTPUT_CONSOLE_PROGRESS_BAR_H__

#include <chrono>

#include "TimeHistoryOutput.h"

int time_history_output_func_console_progress_bar(TimeHistoryOutput &_self);

/* ===========================================================
Class TimeHistoryOutput_ConsoleProgressBar
=========================================================== */
class TimeHistoryOutput_ConsoleProgressBar : public TimeHistoryOutput
{
protected:
	int width; // width must <= 200
	size_t prev_pos, cur_pos;
	float width_div_100;
	
	std::chrono::system_clock::time_point start_time;
	std::chrono::system_clock::time_point end_time;
	
public:
	TimeHistoryOutput_ConsoleProgressBar();
	~TimeHistoryOutput_ConsoleProgressBar();

	// Initialize each steps
	int init_per_step(void) override;
	// output function
	friend int time_history_output_func_console_progress_bar(TimeHistoryOutput &_self);
	// Finalize each steps
	void finalize_per_step(void) override;

	inline void set_width(int wd) noexcept
	{
		if (wd > 0)
		{
			width = wd < 200 ? wd : 200;
			width_div_100 = (float)width / 100.0f;
		}
	}

protected:
	void print_progress(void);
};


#endif