#ifndef _TIMEHISTORY_CONSOLEPROGRESSBAR_H_
#define _TIMEHISTORY_CONSOLEPROGRESSBAR_H_

#include <chrono>
#include "TimeHistory.h"

class TimeHistory_ConsoleProgressBar : public TimeHistory
{
protected:
	int width; // width must <= 200
	size_t prev_pos, cur_pos;
	float width_div_100;
	
	std::chrono::system_clock::time_point start_time;
	std::chrono::system_clock::time_point end_time;
	
public:
	TimeHistory_ConsoleProgressBar();
	~TimeHistory_ConsoleProgressBar();

	// Initialize each steps
	int init_per_step(void);
	// Finalize each steps
	void finalize_per_step(void);
	// Output funtion
	int output(void);

	inline void set_width(int wd)
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