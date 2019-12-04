#include "TestsWithGL_pcp.h"

#include "ColorGraph.h"
#include "test_post_processor.h"

namespace
{
	void print_color(ColorGraph::Colorf &c)
	{
		std::cout << c.r << ", " << c.g << ", " << c.b << "\n";
	}
}

void test_color_graph(void)
{
	ColorGraph cg;
	ColorGraph::ValueColorPair vc_pairs[] = {
		{ 0.00, 0.2f, 0.3f, 0.4f },
		{ 0.20, 0.3f, 0.3f, 0.5f },
		{ 0.30, 0.4f, 0.5f, 0.6f },
		{ 0.45, 0.6f, 0.6f, 0.8f }
	};
	cg.init(vc_pairs, sizeof(vc_pairs) / sizeof(ColorGraph::ValueColorPair));

	ColorGraph::Colorf c;
	c = cg.get_color(-1.0);
	print_color(c);
	c = cg.get_color(1.0);
	print_color(c);

	c = cg.get_color(0.1);
	print_color(c);
	c = cg.get_color(0.25);
	print_color(c);
	c = cg.get_color(0.4);
	print_color(c);

	system("pause");
}
