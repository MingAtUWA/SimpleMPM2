#ifndef __RESULT_FILE_H__
#define __RESULT_FILE_H__

class ResultFile
{
public:
	ResultFile() {}
	~ResultFile() {}
};

class TimeHistory;
typedef int (*OutputFunc) (TimeHistory &_th, ResultFile &_rf);

#endif