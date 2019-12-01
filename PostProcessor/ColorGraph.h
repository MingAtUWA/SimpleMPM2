#ifndef __Color_Graph_H__
#define __Color_Graph_H__

class ColorGraph
{
public:
	struct Color { float r, g, b; };
	struct ValueColorPair { double va; float r, g, b; };

	ColorGraph() : pairs(nullptr), color_num(0) {}
	~ColorGraph() { clear(); }

	int init(ValueColorPair *vcps, size_t num)
	{
		clear();
		if (!vcps || num == 0) return -1;
		pairs = new ValueColorPair[num];
		color_num = num;
		for (size_t i = 0; i < num; ++i)
		{
			pairs[i].va = vcps[i].va;
			pairs[i].r = vcps[i].r;
			pairs[i].g = vcps[i].g;
			pairs[i].b = vcps[i].b;
		}
		return 0;
	}

	inline Color get_color(double va)
	{
		Color res;
		if (va < pairs[0].va)
		{
			// white
			res.r = 1.0;
			res.g = 1.0;
			res.b = 1.0;
			return res;
		}
		else if (va > pairs[color_num - 1].va)
		{
			// grep
			res.r = 0.5;
			res.g = 0.5;
			res.b = 0.5;
			return res;
		}

		// find va is in which internval
		size_t low_id = 0, up_id = color_num - 1, mid_id;
		mid_id = (low_id + up_id) / 2;
		do
		{
			if (pairs[mid_id].va > va)
				up_id  = mid_id;
			else
				low_id = mid_id;
			mid_id = (low_id + up_id) / 2;
		} while (low_id != mid_id);
		
		// interpolate color (can use nonlinear interpolation in the future)
		res.r = pairs[low_id].r + (pairs[up_id].r - pairs[low_id].r) / (pairs[up_id].va - pairs[low_id].va) * (va - pairs[low_id].va);
		res.g = pairs[low_id].g + (pairs[up_id].g - pairs[low_id].g) / (pairs[up_id].va - pairs[low_id].va) * (va - pairs[low_id].va);
		res.b = pairs[low_id].b + (pairs[up_id].b - pairs[low_id].b) / (pairs[up_id].va - pairs[low_id].va) * (va - pairs[low_id].va);
		return res;
	}

protected:
	size_t color_num;
	ValueColorPair *pairs;

	void clear(void)
	{
		if (pairs)
		{
			delete[] pairs;
			pairs = nullptr;
		}
		color_num = 0;
	}
};


#endif